"""Fetch Ensembl's core history tables, xref dumps and BioMart tables over HTTPS (no idtrack involved).

Usage:
  fetch_tables.py --organism <name> core [--release 115]        stable_id_event, mapping_session, gene, meta of the core db
  fetch_tables.py --organism <name> biomart --releases 115,105  whole-genome gene <-> external tables of archive BioMarts
  fetch_tables.py --organism <name> vgnc                         VGNC's gene set of the organism (VGNC id -> symbol)
  fetch_tables.py --organism <name> vega --releases 75,80        Vega gene names of the core xref tables of a release
  fetch_tables.py --organism <name> chains                       Ensembl's assembly chain files between the organism's assemblies
  fetch_tables.py --organism <name> sequences --releases eras     the sequence regions of each assembly and their kind (karyotype)

Writes into `<work>/truth/<organism>/core/`, `biomart/<release>/<group>.tsv`, `vgnc/`, `vega/` and `chain/`, one line
per file in SOURCES.tsv.
Resumable: existing files are skipped. A BioMart group whose attributes the archive lacks gets a `.unavailable` file;
an archive that does not answer gets `biomart/<release>/_archive_unavailable`.
"""
import argparse
import bz2
import datetime
import gzip
import hashlib
import os
import re
import sys
import time
import urllib.parse

from idtrack_reproducibility.harness import harness_lib as H

FTP = "https://ftp.ensembl.org/pub"
# Ensembl archive sites (https://www.ensembl.org/info/website/archives/index.html)
ARCHIVE = {116: "jun2026", 115: "sep2025", 114: "may2025", 113: "oct2024", 112: "may2024", 111: "jan2024",
           110: "jul2023", 109: "feb2023", 108: "oct2022", 107: "jul2022", 106: "apr2022", 105: "dec2021",
           104: "may2021", 103: "feb2021", 102: "nov2020", 101: "aug2020", 100: "apr2020", 99: "jan2020",
           98: "sep2019", 97: "jul2019", 96: "apr2019", 95: "jan2019", 94: "oct2018", 93: "jul2018", 92: "apr2018",
           91: "dec2017", 90: "aug2017", 89: "may2017", 88: "mar2017", 87: "dec2016", 86: "oct2016", 85: "jul2016",
           84: "mar2016", 83: "dec2015", 82: "sep2015", 81: "jul2015", 80: "may2015", 79: "mar2015", 78: "dec2014",
           77: "oct2014", 76: "aug2014", 75: "feb2014", 74: "dec2013", 67: "may2012", 54: "may2009"}
GROUPS = {
    "genes": [["ensembl_gene_id"], ["ensembl_gene_id_version"], ["external_gene_name", "external_gene_id"],
              ["external_gene_source", "external_gene_db"], ["gene_biotype"]],
    "entrez": [["ensembl_gene_id"], ["entrezgene_id", "entrezgene"], ["entrezgene_accession"]],
    "hgnc": [["ensembl_gene_id"], ["hgnc_symbol"], ["hgnc_id"]],
    "mgi": [["ensembl_gene_id"], ["mgi_symbol"], ["mgi_id"]],
    "vgnc": [["ensembl_gene_id"], ["vgnc_symbol", "vgnc"], ["vgnc_id"]],
    "uniprot_gn": [["ensembl_gene_id"], ["uniprot_gn_symbol", "uniprot_genename"], ["uniprot_gn_id"]],
    "swissprot": [["ensembl_gene_id"], ["uniprotswissprot", "uniprot_swissprot", "uniprot_swissprot_accession"]],
    # the corrected judge: BioMart 75's `uniprot_swissprot` is the entry name (display label); the accession is its own
    # attribute there (only archives up to 77 have it; fetched with --only swissprot_acc --releases 75)
    "swissprot_acc": [["ensembl_gene_id"], ["uniprot_swissprot_accession"]],
    "sptrembl": [["ensembl_gene_id"], ["uniprotsptrembl", "uniprot_sptrembl"]],
    "synonym": [["ensembl_gene_id"], ["external_synonym"]],
    "refseq_mrna": [["ensembl_gene_id"], ["refseq_mrna"]],
    "refseq_mrna_predicted": [["ensembl_gene_id"], ["refseq_mrna_predicted"]],
    "refseq_ncrna": [["ensembl_gene_id"], ["refseq_ncrna"]],
    "refseq_ncrna_predicted": [["ensembl_gene_id"], ["refseq_ncrna_predicted"]],
    "refseq_peptide": [["ensembl_gene_id"], ["refseq_peptide"]],
    "refseq_peptide_predicted": [["ensembl_gene_id"], ["refseq_peptide_predicted"]],
    "rfam": [["ensembl_gene_id"], ["rfam"]],
    "wikigene": [["ensembl_gene_id"], ["wikigene_name"], ["wikigene_id"]],
    "transcripts": [["ensembl_gene_id"], ["ensembl_transcript_id"], ["ensembl_peptide_id"]],
}


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def save(org, rel_path, data, url):
    dst = os.path.join(H.truth_dir(org), rel_path)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst + ".tmp", "wb") as fo:
        fo.write(data)
    os.replace(dst + ".tmp", dst)
    H.add_source(org, rel_path, url, len(data), hashlib.md5(data).hexdigest(), now())


def core(org, rel):
    html = H.http_get(f"{FTP}/release-{rel}/mysql/").decode()
    dbs = sorted(set(re.findall(rf'href="({org}_core_{rel}_[^"/]+)/?"', html)))
    if len(dbs) != 1:
        raise SystemExit(f"core databases of {org} at {rel}: {dbs}")
    for t in ["stable_id_event", "mapping_session", "gene", "meta"]:
        rp = f"core/{dbs[0]}.{t}.txt.gz"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        url = f"{FTP}/release-{rel}/mysql/{dbs[0]}/{t}.txt.gz"
        save(org, rp, H.http_get(url), url)
        print("saved", rp, flush=True)


def biomart(org, rels, only=None):
    ds = H.ORGANISMS[org][0]
    for rel in rels:
        base = f"https://{ARCHIVE[rel]}.archive.ensembl.org/biomart/martservice"
        d = os.path.join(H.truth_dir(org), "biomart", str(rel))
        os.makedirs(d, exist_ok=True)
        apath = os.path.join(d, "_attributes.txt")
        if not os.path.exists(apath):
            try:
                txt = H.http_get(f"{base}?type=attributes&dataset={ds}", tries=3, timeout=300)
            except Exception as e:  # noqa: BLE001
                open(os.path.join(d, "_archive_unavailable"), "w").write(f"{base}: {e!r} at {now()}\n")
                print(rel, "archive unavailable", e, flush=True)
                continue
            if b"\t" not in txt[:2000]:
                open(os.path.join(d, "_archive_unavailable"), "w").write(f"{base}: no attribute list for {ds}\n")
                print(rel, "no dataset", flush=True)
                continue
            open(apath, "wb").write(txt)
        avail = {line.split("\t")[0] for line in open(apath) if line.strip()}
        for group, slots in GROUPS.items():
            if only and group not in only:
                continue
            dst = os.path.join(d, f"{group}.tsv")
            if os.path.exists(dst) or os.path.exists(dst + ".unavailable"):
                continue
            chosen = [next(a for a in slot if a in avail) for slot in slots if any(a in avail for a in slot)]
            if len(chosen) < 2 or chosen[0] != "ensembl_gene_id":
                open(dst + ".unavailable", "w").write(f"attributes not in archive: {slots}\n")
                continue
            attrs = "".join(f'<Attribute name="{a}"/>' for a in chosen)
            xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" '
                   'formatter="TSV" header="0" uniqueRows="1" datasetConfigVersion="0.6" completionStamp="1">'
                   f'<Dataset name="{ds}" interface="default">{attrs}</Dataset></Query>')
            ok = False
            lines = ["(no answer)"]
            for attempt in range(3):
                try:
                    txt = H.http_get(base, data=urllib.parse.urlencode({"query": xml}).encode(), timeout=1800)
                except Exception as e:  # noqa: BLE001  a failing table is skipped; a rerun tries it again
                    lines = [repr(e)]
                    time.sleep(30)
                    continue
                lines = txt.decode("utf-8", "replace").rstrip("\n").split("\n")
                if lines and lines[-1].strip() == "[success]":
                    ok = True
                    break
                time.sleep(30)
            if not ok:
                print(rel, group, "INCOMPLETE", lines[-1][:200], flush=True)
                continue
            body = ("\t".join(chosen) + "\n" + "\n".join(lines[:-1]) + "\n").encode()
            save(org, f"biomart/{rel}/{group}.tsv", body, f"{base} (dataset {ds}; attributes {','.join(chosen)})")
            print(rel, group, len(lines) - 1, flush=True)


VGNC = "https://storage.googleapis.com/public-download-files/vgnc/tsv"
VGNC_SPECIES = {"sus_scrofa": "pig", "macaca_mulatta": "macaque", "macaca_fascicularis": "crab-eating_macaque"}


def vgnc(org):
    """VGNC's gene set of the organism (vgnc_id, symbol, prev_symbol, ...): BioMart's VGNC column is the id only."""
    sp = VGNC_SPECIES.get(org)
    if sp is None:
        print(org, "has no VGNC gene set", flush=True)
        return
    rp = f"vgnc/{sp}_vgnc_gene_set_All.tsv"
    if os.path.exists(os.path.join(H.truth_dir(org), rp)):
        return
    url = f"{VGNC}/{sp}/{sp}_vgnc_gene_set_All.tsv"
    save(org, rp, H.http_get(url), url)
    print("saved", rp, flush=True)


def _columns(sql, table):
    m = re.search(rf"CREATE TABLE `{table}` \((.*?)\n\)", sql, re.S)
    return [c for c in re.findall(r"^\s+`(\w+)`", m.group(1), re.M)] if m else None


def vega(org, rels):
    """Vega gene names of a release: core xrefs of `external_db.db_name` `Vega_gene` on genes (display label).

    The dumps of `external_db`, `xref`, `object_xref` and `gene` (and `gene_stable_id` where the gene table has no
    stable id) are read from the FTP; each is listed in SOURCES.tsv (URL, size, md5, time) and not kept; the
    extract `vega/<release>.tsv` (gene stable id, Vega name, OTTG accession) is listed too.
    """
    for rel in rels:
        rp = f"vega/{rel}.tsv"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        html = H.http_get(f"{FTP}/release-{rel}/mysql/").decode()
        dbs = sorted(set(re.findall(rf'href="({org}_core_{rel}_[^"/]+)/?"', html)))
        if len(dbs) != 1:
            raise SystemExit(f"core databases of {org} at {rel}: {dbs}")
        base = f"{FTP}/release-{rel}/mysql/{dbs[0]}"
        sql = gzip.decompress(H.http_get(f"{base}/{dbs[0]}.sql.gz")).decode()

        def table(name):
            url = f"{base}/{name}.txt.gz"
            data = H.http_get(url)
            H.add_source(org, f"{rp} (read from {name}.txt.gz; not kept)", url, len(data), hashlib.md5(data).hexdigest(), now())
            cols = _columns(sql, name)
            held = ""
            for line in gzip.decompress(data).decode("utf-8", "replace").split("\n"):
                if line.endswith("\\") and not line.endswith("\\\\"):  # an escaped newline inside a text field
                    held += line[:-1] + "\n"
                    continue
                line, held = held + line, ""
                fields = line.split("\t")
                if line and len(fields) == len(cols):
                    yield dict(zip(cols, fields))

        ids = {r["external_db_id"] for r in table("external_db") if r["db_name"] == "Vega_gene"}
        rows = []
        if ids:
            xrefs = {r["xref_id"]: (r["display_label"], r["dbprimary_acc"]) for r in table("xref") if r["external_db_id"] in ids}
            links = [(r["ensembl_id"], r["xref_id"]) for r in table("object_xref")
                     if r["ensembl_object_type"] == "Gene" and r["xref_id"] in xrefs]
            if "stable_id" in (_columns(sql, "gene") or []):
                stable = {r["gene_id"]: r["stable_id"] for r in table("gene")}
            else:
                stable = {r["gene_id"]: r["stable_id"] for r in table("gene_stable_id")}
            rows = sorted({(stable[g], *xrefs[x]) for g, x in links if g in stable})
        body = ("gene_id\tvega_name\tottg\n" + "".join("\t".join(r) + "\n" for r in rows)).encode()
        save(org, rp, body, f"{base} (external_db Vega_gene, object_xref Gene, xref display_label)")
        print(rel, "Vega gene names", len(rows), flush=True)


EXC_HEADER = ["seqname", "start", "end", "type", "ref_seqname", "ref_start", "ref_end", "ori", "seq_length", "ref_length",
              "from_release"]


# The GRC's own placement of alternative loci and patches on the chromosomes, for the releases whose Ensembl tables
# hold none (human GRCh38 from release 110: the dumps of assembly_exception are empty, and the patch and
# alternative-locus sequences are scaffolds with coordinates of their own; their GRC names are Ensembl's).
GRC_PLACEMENTS = {
    "homo_sapiens": {"GRCh38": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/"
        "GCA_000001405.29_GRCh38.p14_assembly_structure/all_alt_scaffold_placement.txt",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/"
        "GCA_000001405.29_GRCh38.p14_assembly_structure/PATCHES/alt_scaffolds/patch_type",
        "GRC GCA_000001405.29 (GRCh38.p14)")},
}


# The GRC's alignments of every alternative locus and patch scaffold to its chromosome (GFF3 `match` lines with
# `Target` and `Gap`), per assembly; `exceptions` places a copy inside its segment through them (assembly_map).
GRC_STRUCTURE = {
    "homo_sapiens": {
        "GRCh37": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/"
                  "GCA_000001405.14_GRCh37.p13_assembly_structure",
        "GRCh38": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/"
                  "GCA_000001405.29_GRCh38.p14_assembly_structure"},
    "mus_musculus": {
        "GRCm38": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_GRCm38.p6/"
                  "GCA_000001635.8_GRCm38.p6_assembly_structure"},
}
ALN_HEADER = ["alt_name", "alt_acc", "alt_length", "parent", "parent_start", "parent_stop", "ref_start", "ref_end",
              "tgt_start", "tgt_end", "tgt_strand", "gap"]


def grc_alignments(org):
    """`grc_alignments/<assembly>.tsv`: one row per aligned `match` of each alt scaffold of the GRC's placement file
    (ALT_REF_LOCI_* and PATCHES; mouse strain loci are no Ensembl exceptions and are skipped). Every file read is
    listed in SOURCES.tsv; the GFFs are not kept."""
    for asm, base in GRC_STRUCTURE.get(org, {}).items():
        rp = f"grc_alignments/{asm}.tsv"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        place_url = f"{base}/all_alt_scaffold_placement.txt"
        data = H.http_get(place_url)
        H.add_source(org, f"{rp} (read from all_alt_scaffold_placement.txt; not kept)", place_url, len(data),
                     hashlib.md5(data).hexdigest(), now())
        rows = []
        for line in data.decode().split("\n"):
            f = line.split("\t")
            if len(f) < 13 or line.startswith("#") or f[4] != "CHROMOSOME":
                continue
            if not (f[0].startswith("ALT_REF_LOCI") or f[0] == "PATCHES"):
                continue
            url = f"{base}/{f[0]}/alt_scaffolds/alignments/{f[3]}_{f[6]}.gff"
            try:
                gff = H.http_get(url, tries=3)
            except Exception as e:  # noqa: BLE001
                print(org, asm, f[2], "no alignment file", url, repr(e)[:80], flush=True)
                continue
            H.add_source(org, f"{rp} (read from {f[3]}_{f[6]}.gff; not kept)", url, len(gff),
                         hashlib.md5(gff).hexdigest(), now())
            for g in gff.decode().split("\n"):
                c = g.split("\t")
                if len(c) < 9 or g.startswith("#") or c[2] != "match":
                    continue
                attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
                tgt = attrs.get("Target", "").split()
                if len(tgt) < 4:
                    continue
                rows.append([f[2], f[3], f[10], f[5], f[11], f[12], c[3], c[4], tgt[1], tgt[2], tgt[3],
                             attrs.get("Gap", f"M{int(c[4]) - int(c[3]) + 1}")])
        body = ("\t".join(ALN_HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows)).encode()
        save(org, rp, body, f"{base} (alt_scaffolds/alignments/*.gff, match lines)")
        print(org, asm, "alignment rows", len(rows), flush=True)


def core_places(org, rels):
    """`core_places/<release>.tsv`: the place (sequence, start, end, strand) of every gene of the release's core
    database that its GTF lacks (human 48-51: genes on unplaced scaffolds; cynomolgus 112: 47 genes), from the core's
    `gene` (with `gene_stable_id` where the gene table has no stable id) and `seq_region` dumps (listed in SOURCES.tsv,
    not kept). The judge places such genes with it (score.Places.core, lines revised by `core_places`)."""
    for rel in rels:
        rp = f"core_places/{rel}.tsv"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        html = H.http_get(f"{FTP}/release-{rel}/mysql/").decode()
        dbs = sorted({d for d in re.findall(rf'href="({org}_core_{rel}_[^"/]+?)/?"', html) if not d.endswith(".gz")})
        db = max(dbs, key=lambda d: int(re.search(r"_(\d+)[a-z]*$", d).group(1)))
        base = f"{FTP}/release-{rel}/mysql/{db}"
        listing = H.http_get(base + "/").decode()

        def table(name):
            url = f"{base}/{name}.txt.gz"
            data = H.http_get(url)
            H.add_source(org, f"{rp} (read from {name}.txt.gz; not kept)", url, len(data), hashlib.md5(data).hexdigest(), now())
            return [line.split("\t") for line in gzip.decompress(data).decode("utf-8", "replace").split("\n") if line]

        gene = table("gene")
        if 'href="gene_stable_id.txt.gz"' in listing:
            stable = {r[0]: (r[1], r[2]) for r in table("gene_stable_id")}
            genes = [(stable[r[0]][0], stable[r[0]][1], r[3], r[4], r[5], r[6]) for r in gene if r[0] in stable]
        else:
            genes = [(r[12], r[13], r[3], r[4], r[5], r[6]) for r in gene]
        sr = {r[0]: r[1] for r in table("seq_region")}
        gtf = H.gtf_genes(org, rel) or {}
        rows = sorted({(g, v, sr.get(s, ""), a, b, "+" if st == "1" else "-") for g, v, s, a, b, st in genes if g not in gtf})
        body = ("gene_id\tgene_version\tseqname\tstart\tend\tstrand\n" + "".join("\t".join(r) + "\n" for r in rows)).encode()
        save(org, rp, body, f"{base} (gene, seq_region: genes the GTF lacks)")
        print(org, rel, "core-only genes with places", len(rows), flush=True)


def grc_placements(org, rp, grc, lengths):
    """Rows (EXC_HEADER) from the GRC's alt_scaffold_placement: an alternative locus is HAP, a patch PATCH_FIX or
    PATCH_NOVEL (patch_type); the scaffold segment [alt_scaf_start, alt_scaf_stop] stands for the parent chromosome's
    [parent_start, parent_stop] (ori + or -)."""
    url, purl, label = grc
    data, pdata = H.http_get(url), H.http_get(purl)
    for u, d in ((url, data), (purl, pdata)):
        H.add_source(org, f"{rp} (read from {u.rsplit('/', 1)[1]}; not kept)", u, len(d), hashlib.md5(d).hexdigest(), now())
    ptype = {}
    for line in pdata.decode().split("\n"):
        f = line.split("\t")
        if len(f) >= 3 and not line.startswith("#"):
            ptype[f[0]] = f[2]
    rows = []
    for line in data.decode().split("\n"):
        f = line.split("\t")
        if len(f) < 13 or line.startswith("#") or f[4] != "CHROMOSOME":
            continue
        kind = "HAP" if f[0].startswith("ALT_REF_LOCI") else f"PATCH_{ptype.get(f[2], 'FIX')}"
        rows.append([f[2], f[9], f[10], kind, f[5], f[11], f[12], "1" if f[8] == "+" else "-1",
                     lengths.get(f[2], ""), lengths.get(f[5], ""), label])
    return rows


def core_reader(org, rel, rp):
    """(the core database's folder on the FTP, a reader of its table dumps) of an organism at a release: `table(name)`
    returns the rows of `<name>.txt.gz` (or `.txt.gz.bz2`), each source listed in SOURCES.tsv under `rp` (not kept)."""
    html = H.http_get(f"{FTP}/release-{rel}/mysql/").decode()
    dbs = sorted(set(re.findall(rf'href="({org}_core_{rel}_[^"/]+?)/?"', html)) - {''})
    dbs = [d for d in dbs if not d.endswith(".gz")]
    if len(dbs) > 1:  # human 112 on also serves the GRCh37 archive core (`_37`): the release's own is the newest
        dbs = [max(dbs, key=lambda d: int(re.search(r"_(\d+)[a-z]*$", d).group(1)))]
    if len(dbs) != 1:
        raise SystemExit(f"core databases of {org} at {rel}: {dbs}")
    base = f"{FTP}/release-{rel}/mysql/{dbs[0]}"
    listing = H.http_get(base + "/").decode()

    def table(name):
        # most releases dump `<table>.txt.gz`; some (human and mouse 89) `<table>.txt.gz.bz2`
        fn = f"{name}.txt.gz" if f'href="{name}.txt.gz"' in listing else f"{name}.txt.gz.bz2"
        url = f"{base}/{fn}"
        data = H.http_get(url)
        H.add_source(org, f"{rp} (read from {fn}; not kept)", url, len(data), hashlib.md5(data).hexdigest(), now())
        raw = bz2.decompress(data) if fn.endswith(".bz2") else data
        raw = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
        return [line.split("\t") for line in raw.decode().split("\n") if line]
    return base, table


SEQ_HEADER = ["seqname", "coord_system", "coord_system_rank", "length", "karyotype_rank", "codon_table"]


def sequences(org, rels):
    """The sequence regions of the default coordinate systems at a release, with what Ensembl's core says of each:
    the coordinate system (`coord_system`, its `rank`), the length, the `karyotype_rank` attribute (a sequence of the
    organism's karyotype: a chromosome, or the mitochondrial genome) and the `codon_table` attribute (another genetic
    code than the standard one: the mitochondrial genome). Written as `sequences/<release>.tsv` (SEQ_HEADER); the
    dumps (`coord_system`, `seq_region`, `seq_region_attrib`, `attrib_type`) are listed in SOURCES.tsv and not kept.
    One release per assembly is enough (a sequence's kind is the assembly's): `--releases eras` takes the last GTF
    release of each assembly of the window."""
    for rel in rels:
        rp = f"sequences/{rel}.tsv"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        base, table = core_reader(org, rel, rp)
        cs = {c[0]: (c[2], c[4]) for c in table("coord_system") if "default_version" in c[-1]}
        codes = {r[0]: r[1] for r in table("attrib_type")}
        want = {k for k, v in codes.items() if v in ("karyotype_rank", "codon_table")}
        sr = {r[0]: [r[1], cs[r[2]][0], cs[r[2]][1], r[3], "", ""] for r in table("seq_region") if r[2] in cs}
        for r in table("seq_region_attrib"):
            if r[1] in want and r[0] in sr:
                sr[r[0]][4 if codes[r[1]] == "karyotype_rank" else 5] = r[2]
        rows = sorted(sr.values())
        body = ("\t".join(SEQ_HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows)).encode()
        save(org, rp, body, f"{base} (coord_system, seq_region, seq_region_attrib, attrib_type; see fetch_tables.sequences)")
        print(org, rel, "sequences", len(rows), "karyotype", sum(1 for r in rows if r[4]), flush=True)


def exceptions(org, rels):
    """Assembly exceptions of a release: the core database's `assembly_exception` (haplotypes HAP, patches
    PATCH_FIX/PATCH_NOVEL, pseudoautosomal regions PAR: a segment of one sequence region and the segment of the
    reference sequence it stands for), with sequence names from `seq_region`. Written as
    `exceptions/<release>.tsv` (EXC_HEADER; `from_release` = the release whose table the row comes from). The dumps
    are listed in SOURCES.tsv and not kept.

    Ensembl's FTP dumps of `assembly_exception` are empty for human 110 to 116 although those releases have patches
    and alternative loci (their GTFs place genes on them). For a release whose dump is empty, the rows of the latest
    earlier release on the same assembly with rows are taken where both sequences exist in this release's
    `seq_region` (default coordinate systems) with the same lengths (`from_release` says which: human 110 on keeps
    only the PARs, since from 110 the patches and alternative loci are scaffolds of their own, named as the GRC
    names them, e.g. `HSCHR6_MHC_COX_CTG1`, `HG1012_PATCH`); patch sequences of the GTF still unplaced are placed by
    the GRC's own `alt_scaffold_placement` (GRC_PLACEMENTS; `from_release` names it). A release with no earlier rows
    (mouse GRCm39, the organisms without exceptions) is written empty.
    """
    asm = H.gtf_assemblies(org)
    for rel in rels:
        rp = f"exceptions/{rel}.tsv"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        base, table = core_reader(org, rel, rp)
        ae = table("assembly_exception")
        rows = []
        if ae:
            sr = {r[0]: (r[1], r[3]) for r in table("seq_region")}
            for r in ae:
                (n1, l1), (n2, l2) = sr[r[1]], sr[r[5]]
                rows.append([n1, r[2], r[3], r[4], n2, r[6], r[7], r[8], l1, l2, str(rel)])
        else:
            earlier = [x for x in range(rel) if asm.get(x) == asm.get(rel)
                       and os.path.exists(os.path.join(H.truth_dir(org), f"exceptions/{x}.tsv"))]
            prev = []
            for x in reversed(earlier):
                prev = [r for r in H.read_tsv(os.path.join(H.truth_dir(org), f"exceptions/{x}.tsv"))
                        if r["from_release"] == str(x)]
                if prev:
                    break
            if prev:
                # sequences of this release's default coordinate systems (the archive's older assemblies are kept
                # in other coordinate systems under the same names, e.g. GRCh37's chromosome Y)
                default = {c[0] for c in table("coord_system") if "default_version" in c[-1]}
                lengths = {r[1]: r[3] for r in table("seq_region") if r[2] in default}
                kept = [r for r in prev if lengths.get(r["seqname"]) == r["seq_length"]
                        and lengths.get(r["ref_seqname"]) == r["ref_length"]]
                rows = [[r[k] for k in EXC_HEADER] for r in kept]
                print(rel, f"dump empty: {len(kept)} of {len(prev)} rows of release {prev[0]['from_release']} "
                      f"(sequences of the same names and lengths)", flush=True)
                grc = GRC_PLACEMENTS.get(org, {}).get(asm.get(rel))
                placed = {r[0] for r in rows}
                gtf = {g["seqname"] for g in (H.gtf_genes(org, rel) or {}).values() if g["region"] == "patch"}
                if grc and gtf - placed:
                    add = grc_placements(org, rp, grc, lengths)
                    rows += [r for r in add if r[0] not in placed]
                    print(rel, f"{len(gtf - placed)} patch sequences of the GTF unplaced by Ensembl's tables: "
                          f"{len(add)} placements of {grc[2]}", flush=True)
                left = gtf - {r[0] for r in rows}
                if left:
                    print(rel, f"patch sequences of the GTF with no placement: {len(left)}, e.g. {sorted(left)[:5]}",
                          flush=True)
        body = ("\t".join(EXC_HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows)).encode()
        save(org, rp, body, f"{base} (assembly_exception with seq_region names; see fetch_tables.exceptions)")
        print(org, rel, "exceptions", len(rows), flush=True)


# Ensembl's assembly chain files (UCSC chain format, Ensembl sequence names). The newest dump of each organism:
# human and mouse in release-115/assembly_chain; pig (Sscrofa9 <-> Sscrofa10.2) last in release-89; none for
# Sscrofa10.2 -> Sscrofa11.1, rhesus or cynomolgus (checked in every release folder 76 to 116 and in assembly_mapping/).
CHAINS = {
    "homo_sapiens": ("release-115/assembly_chain/homo_sapiens", ["NCBI36_to_GRCh37", "GRCh37_to_NCBI36", "GRCh37_to_GRCh38",
                                                                  "GRCh38_to_GRCh37", "NCBI36_to_GRCh38", "GRCh38_to_NCBI36"]),
    "mus_musculus": ("release-115/assembly_chain/mus_musculus", ["NCBIM37_to_GRCm38", "GRCm38_to_NCBIM37", "GRCm38_to_GRCm39",
                                                                 "GRCm39_to_GRCm38"]),
    "sus_scrofa": ("release-89/assembly_chain/sus_scrofa", ["Sscrofa9_to_Sscrofa10.2", "Sscrofa10.2_to_Sscrofa9"]),
}


def chains(org):
    if org not in CHAINS:
        print(org, "has no assembly chain file on Ensembl's FTP", flush=True)
        return
    folder, names = CHAINS[org]
    for n in names:
        rp = f"chain/{n}.chain.gz"
        if os.path.exists(os.path.join(H.truth_dir(org), rp)):
            continue
        url = f"{FTP}/{folder}/{n}.chain.gz"
        save(org, rp, H.http_get(url), url)
        print("saved", rp, flush=True)


def eras_of(org):
    """[(assembly, first release, last release)] of the GTF releases of the window."""
    out = []
    asm = H.gtf_assemblies(org)
    for r in H.gtf_releases(org):
        if r > H.ORGANISMS[org][2]:
            continue
        if out and out[-1][0] == asm.get(r):
            out[-1][2] = r
        else:
            out.append([asm.get(r), r, r])
    return [tuple(x) for x in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("what", choices=["core", "biomart", "vgnc", "vega", "chains", "exceptions", "grc_alignments",
                                     "core_places", "sequences"])
    ap.add_argument("--release", type=int, default=None)
    ap.add_argument("--releases", default="115")
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    if a.what == "core":
        core(a.organism, a.release or H.ORGANISMS[a.organism][2])
    elif a.what == "vgnc":
        vgnc(a.organism)
    elif a.what == "vega":
        vega(a.organism, [int(x) for x in a.releases.split(",")])
    elif a.what == "exceptions":
        rr = H.gtf_releases(a.organism) if a.releases == "all" else [int(x) for x in a.releases.split(",")]
        exceptions(a.organism, sorted(rr))
    elif a.what == "sequences":
        if a.releases == "eras":
            rr = [hi for _a, _lo, hi in eras_of(a.organism)]
        else:
            rr = [int(x) for x in a.releases.split(",")]
        sequences(a.organism, sorted(rr))
    elif a.what == "core_places":
        core_places(a.organism, [int(x) for x in a.releases.split(",")])
    elif a.what == "grc_alignments":
        grc_alignments(a.organism)
    elif a.what == "chains":
        chains(a.organism)
    else:
        biomart(a.organism, [int(x) for x in a.releases.split(",")], a.only.split(",") if a.only else None)
    print("done", flush=True)


if __name__ == "__main__":
    main()
