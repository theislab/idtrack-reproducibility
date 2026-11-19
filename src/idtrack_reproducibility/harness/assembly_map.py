"""Gene places mapped between assemblies through Ensembl's assembly chain files (`truth/<organism>/chain/`).

Imports nothing of idtrack. A chain file (UCSC chain format, Ensembl's sequence names) aligns blocks of one assembly
(the chain's target, "from") to another (its query, "to"). A place (sequence, start, end, strand; GTF coordinates,
1-based and inclusive) is mapped by taking every aligned block that overlaps it, grouping the mapped pieces by the
sequence and strand they land on, and keeping the group with the most aligned bases: the mapped place spans that
group's pieces, on the gene's strand (flipped where the chain aligns to the reverse strand). Where no direct chain
exists between two assemblies but one exists through a third (mouse NCBIM37 <-> GRCm39 through GRCm38), the place
is mapped twice. A place is unmapped when no chain path exists, or when fewer than `MIN_ALIGNED` of its bases are
aligned in the kept group (the reason says which).
"""
import bisect
import collections
import gzip
import os

from idtrack_reproducibility.harness import harness_lib as H

MIN_ALIGNED = 0.5  # share of the gene's bases that the kept group must align (checked in `check_handoffs`)


class Chain:
    """The aligned blocks of one chain file, per sequence of the 'from' assembly, for overlap queries."""

    def __init__(self, path):
        blocks = collections.defaultdict(list)  # tname -> [(tstart, tend, qname, qstrand, qsize, qstart)]
        with gzip.open(path, "rt") as fh:
            head = None
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if f[0] == "chain":
                    tname, tstart = f[2], int(f[5])
                    qname, qsize, qstrand, qstart = f[7], int(f[8]), f[9], int(f[10])
                    head = [tname, tstart, qname, qsize, qstrand, qstart]
                    continue
                size = int(f[0])
                tname, t, qname, qsize, qstrand, qq = head
                blocks[tname].append((t, t + size, qname, qstrand, qsize, qq))
                if len(f) == 3:
                    head[1] = t + size + int(f[1])
                    head[5] = qq + size + int(f[2])
        self.index = {}
        for tname, bl in blocks.items():
            bl.sort()
            maxend, m = [], -1
            for b in bl:
                m = max(m, b[1])
                maxend.append(m)
            self.index[tname] = (bl, maxend)

    def map(self, seqname, start, end, strand):
        """Map a 1-based inclusive place; returns (place or None, aligned share, reason)."""
        if seqname not in self.index:
            return None, 0.0, f"sequence {seqname} not in the chain"
        bl, maxend = self.index[seqname]
        s0, e0 = start - 1, end
        i = bisect.bisect_right(maxend, s0)
        groups = collections.defaultdict(lambda: [0, None, None])
        while i < len(bl) and bl[i][0] < e0:
            t1, t2, qname, qstrand, qsize, q1 = bl[i]
            a, b = max(t1, s0), min(t2, e0)
            if a < b:
                qa, qb = q1 + (a - t1), q1 + (b - t1)
                if qstrand == "-":  # reverse-strand coordinates of the query sequence -> forward coordinates
                    qa, qb = qsize - qb, qsize - qa
                g = groups[(qname, qstrand)]
                g[0] += b - a
                g[1] = qa if g[1] is None else min(g[1], qa)
                g[2] = qb if g[2] is None else max(g[2], qb)
            i += 1
        if not groups:
            return None, 0.0, "no aligned block"
        (qname, qstrand), (n, qa, qb) = max(groups.items(), key=lambda kv: kv[1][0])
        share = n / (e0 - s0)
        if share < MIN_ALIGNED:
            return None, share, f"{share:.0%} of the place aligned"
        new_strand = strand if qstrand == "+" else {"+": "-", "-": "+"}.get(strand, strand)
        return dict(seqname=qname, start=qa + 1, end=qb, strand=new_strand), share, ""


class Mapper:
    """All chain paths of one organism between the assemblies of its GTF releases."""

    def __init__(self, org):
        self.org = org
        d = os.path.join(H.truth_dir(org), "chain")
        self.files = {}
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.endswith(".chain.gz") and "_to_" in f:
                    a, b = f[: -len(".chain.gz")].split("_to_")
                    self.files[(a, b)] = os.path.join(d, f)
        self.loaded = {}

    def chain(self, a, b):
        if (a, b) not in self.loaded:
            self.loaded[(a, b)] = Chain(self.files[(a, b)])
        return self.loaded[(a, b)]

    def path(self, a, b):
        """The assemblies to pass through from a to b (direct, else through one other), or None."""
        if (a, b) in self.files:
            return [a, b]
        for (x, y) in self.files:
            if x == a and (y, b) in self.files:
                return [a, y, b]
        return None

    def map(self, place, to_assembly):
        """Map a place dict (assembly, seqname, start, end, strand) to another assembly.

        Returns a dict: the mapped place with `assembly`, `via` (the chain path) and `aligned` (the smallest aligned
        share along the path), or `unmapped` with the reason.
        """
        a = place["assembly"]
        p = self.path(a, to_assembly)
        if p is None:
            return dict(unmapped=f"no assembly chain from {a} to {to_assembly} on Ensembl's FTP")
        cur, low = dict(seqname=place["seqname"], start=int(place["start"]), end=int(place["end"]),
                        strand=place["strand"]), 1.0
        for x, y in zip(p, p[1:]):
            new, share, why = self.chain(x, y).map(cur["seqname"], cur["start"], cur["end"], cur["strand"])
            if new is None:
                return dict(unmapped=f"{x} -> {y}: {why}")
            cur, low = new, min(low, share)
        return dict(assembly=to_assembly, via=" -> ".join(p), aligned=round(low, 3), **cur)


class Exceptions:
    """Ensembl's assembly exceptions of each release (`truth/<organism>/exceptions/<release>.tsv`, fetch_tables.py
    exceptions): a haplotype, patch or alternative locus is a full-length copy of its reference chromosome whose
    segment [start, end] replaces the reference segment [ref_start, ref_end]; a pseudoautosomal region (PAR) of Y is
    the X segment it equals. From human 110 on (scaffolds with coordinates of their own) the rows are the GRC's
    placements of the scaffold segment on its parent chromosome. `project` puts a place on such a sequence onto the
    reference sequence: a PAR place by its offset (the segments are equal base for base); a HAP/PATCH place outside
    the segment by the offset of the nearer segment end, inside it through the GRC's alignment of the scaffold to
    its chromosome (`grc_alignments/<assembly>.tsv`, `aligned_by`). An end inside the segment that no alignment
    places (a part of the scaffold aligned to nothing, or a scaffold without alignment: NCBI36's haplotypes, five
    GRCh38 scaffolds) makes the place `unplaced` (the reason): the judge then cannot compare places (undetermined).
    Two earlier versions placed such ends by linear scaling, then by the span of both segment ends; both were wrong
    where the copy's sequence is rearranged against the reference (NCBI36 c6_COX genes 30 kb away)."""

    def __init__(self, org):
        self.org, self.cache, self.aln = org, {}, {}
        self.asm = H.gtf_assemblies(org)

    def alignments(self, asm):
        """The GRC's alignments of the alt scaffolds of an assembly (`grc_alignments/<assembly>.tsv`, fetch_tables.py
        grc_alignments): {scaffold name: [rows]} and {(parent, parent_start, length): name} (Ensembl's GRCh37 and
        pre-110 names differ from the GRC's: `HSCHR6_MHC_COX` for `HSCHR6_MHC_COX_CTG1`, `CHR_HG1012_PATCH`)."""
        if asm not in self.aln:
            p = os.path.join(H.truth_dir(self.org), "grc_alignments", f"{asm}.tsv")
            by, key = collections.defaultdict(list), {}
            if os.path.exists(p):
                for r in H.read_tsv(p):
                    by[r["alt_name"]].append(r)
                    key[(r["parent"], r["parent_start"], r["alt_length"])] = r["alt_name"]
            self.aln[asm] = (dict(by), key)
        return self.aln[asm]

    @staticmethod
    def _map_pos(rows, sp):
        """A scaffold position through the GRC's match lines (GFF3 Gap: M both advance, I the scaffold, D the
        chromosome) to (chromosome position, target strand), or None outside every aligned match."""
        for m in rows:
            t1, t2 = int(m["tgt_start"]), int(m["tgt_end"])
            if not t1 <= sp <= t2:
                continue
            plus = m["tgt_strand"] != "-"
            ref, tpos = int(m["ref_start"]), (t1 if plus else t2)
            for op in m["gap"].split():
                kind, n = op[0], int(op[1:])
                if kind in "MI":
                    hit = (tpos <= sp < tpos + n) if plus else (tpos - n < sp <= tpos)
                    if hit:
                        return (ref + (sp - tpos if plus else tpos - sp) if kind == "M" else ref), m["tgt_strand"]
                    tpos = tpos + n if plus else tpos - n
                    if kind == "M":
                        ref += n
                elif kind == "D":
                    ref += n
        return None

    def _aligned(self, place, r, rel):
        """The place on its chromosome through the GRC's alignment of its scaffold: an end inside the exception
        segment through the alignment, an end outside it by the offset of the segment end. None when the GRC has no
        alignment of the scaffold; a place with `unplaced` (the reason) when an end inside the segment lies in a part
        of the scaffold that the GRC aligns to nothing (novel sequence, or rearranged beyond the alignment)."""
        asm = self.asm.get(rel)
        by, key = self.alignments(asm) if asm else ({}, {})
        grc_coords = not str(r["from_release"]).isdigit()  # rows from the GRC's placement: scaffold coordinates
        s1, s2, r1, r2 = int(r["start"]), int(r["end"]), int(r["ref_start"]), int(r["ref_end"])
        name = r["seqname"] if r["seqname"] in by else r["seqname"][4:] if r["seqname"].startswith("CHR_") and \
            r["seqname"][4:] in by else key.get((r["ref_seqname"], r["ref_start"], str(s2 - s1 + 1)))
        if name not in by:
            name = None
        pos, strands = [], []
        for x in (int(place["start"]), int(place["end"])):
            if x < s1:
                pos.append(x - s1 + r1)
            elif x > s2:
                pos.append(x - s2 + r2)
            else:
                got = self._map_pos(by[name], x if grc_coords else x - s1 + 1) if name else None
                if got is None:
                    why = (f"{place['seqname']}:{x} lies in a part of scaffold {name} that the GRC aligns to nothing "
                           f"on {r['ref_seqname']}") if name else (f"{place['seqname']}:{x} lies inside the exception "
                                                                   f"segment and the GRC has no alignment of it ({asm})")
                    return dict(place, projected_from=place["seqname"], exception=r["type"], unplaced=why)
                pos.append(got[0])
                strands.append(got[1])
        strand = place["strand"]
        if strands and strands[0] == "-":
            strand = "-" if strand == "+" else "+"
        lo, hi = sorted(pos)
        return dict(place, seqname=r["ref_seqname"], start=lo, end=hi, strand=strand, projected_from=place["seqname"],
                    exception=r["type"], aligned_by="GRC alignment")

    def rows(self, rel):
        if rel not in self.cache:
            p = os.path.join(H.truth_dir(self.org), "exceptions", f"{rel}.tsv")
            by = collections.defaultdict(list)
            if os.path.exists(p):
                for r in H.read_tsv(p):
                    by[r["seqname"]].append(r)
            self.cache[rel] = dict(by)
        return self.cache[rel]

    def has(self, rel):
        return bool(self.rows(rel))

    def project(self, place, rel):
        """The place (dict with seqname, start, end, strand; other keys kept) on the reference sequence, with
        `projected_from` (the sequence it was on) and `exception` (its type); the place itself when its sequence is no
        exception at that release (or a Y place outside the PARs)."""
        rows = self.rows(rel).get(place["seqname"]) if place else None
        if not rows:
            return place
        a, b = int(place["start"]), int(place["end"])
        par = [r for r in rows if r["type"] == "PAR" and int(r["start"]) <= b and int(r["end"]) >= a]
        hp = [r for r in rows if r["type"] != "PAR"]
        if par:
            r = par[0]
        elif hp:
            # the segment the place overlaps most (a sequence has one HAP/PATCH segment in Ensembl's tables)
            r = max(hp, key=lambda r: min(b, int(r["end"])) - max(a, int(r["start"])))
        else:
            return place
        if r["type"] != "PAR":
            return self._aligned(place, r, rel)
        # a PAR of Y is the X segment it equals, base for base: the offset
        s1, r1 = int(r["start"]), int(r["ref_start"])
        lo, hi = a - s1 + r1, b - s1 + r1
        strand = place["strand"]
        return dict(place, seqname=r["ref_seqname"], start=lo, end=hi, strand=strand, projected_from=place["seqname"],
                    exception=r["type"])


def overlaps(place, g):
    """A GTF gene row overlaps a place on the same sequence and strand."""
    return bool(g and g["seqname"] == place["seqname"] and g["strand"] == place["strand"]
                and int(g["start"]) <= int(place["end"]) and int(g["end"]) >= int(place["start"]))


def check_handoffs(org):
    """The mapping checked on genes whose id lives across each assembly handoff: the place at the last release of
    one assembly, mapped to the next, must overlap the place the GTF gives at the first release of the next (same
    sequence and strand). Both directions. Returns rows (handoff, direction, genes, mapped, overlapping, shares)."""
    asm = H.gtf_assemblies(org)
    rels = sorted(r for r in H.gtf_releases(org) if r <= H.ORGANISMS[org][2])
    m = Mapper(org)
    rows = []
    for r1, r2 in zip(rels, rels[1:]):
        a1, a2 = asm.get(r1), asm.get(r2)
        if a1 == a2:
            continue
        g1 = H.gtf_genes(org, r1, "primary")
        g2 = H.gtf_genes(org, r2, "primary")
        both = sorted(set(g1) & set(g2))
        for (ra, aa, ga), (rb, ab, gb) in (((r1, a1, g1), (r2, a2, g2)), ((r2, a2, g2), (r1, a1, g1))):
            n = mapped = hit = 0
            why = collections.Counter()
            for g in both:
                n += 1
                p = dict(assembly=aa, seqname=ga[g]["seqname"], start=int(ga[g]["start"]), end=int(ga[g]["end"]),
                         strand=ga[g]["strand"])
                res = m.map(p, ab)
                if "unmapped" in res:
                    why[res["unmapped"].split(":")[-1].strip() if "%" not in res["unmapped"] else "under half aligned"] += 1
                    continue
                mapped += 1
                hit += overlaps(res, gb[g])
            rows.append(dict(handoff=f"{aa} {ra} -> {ab} {rb}", genes=n, mapped=mapped, overlapping=hit,
                             share_mapped=round(mapped / n, 4) if n else None,
                             share_overlapping_of_mapped=round(hit / mapped, 4) if mapped else None,
                             unmapped=dict(why.most_common(4))))
    return rows


if __name__ == "__main__":
    import json
    import sys
    for o in sys.argv[1:]:
        for row in check_handoffs(o):
            print(json.dumps(dict(organism=o, **row)), flush=True)
