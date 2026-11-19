"""The weakness hunts' reproducers: one table of calls, one command (docs/hunts.md, "Reproducing").

Usage:
  python -m idtrack_reproducibility.hunts.reproduce <angle> <id> [--tree d4c2dda] [--organism o] [--json file]
  python -m idtrack_reproducibility.hunts.reproduce <angle|all> all ...    every reproducer of the angle (or of all)
  python -m idtrack_reproducibility.hunts.reproduce --list                 the table, without a graph
  python -m idtrack_reproducibility.hunts.reproduce <angle|all> <id|all> --submit [--tree d4c2dda,71690ab] [--out d]
                                                   one job per organism (jobs.submit), each tree in turn

Each finding of the hunts (angles `inputs`, `history`, `regions`) has here the calls the hunt made, the answer the
truth expects with the truth's source, and the answer the hunt recorded. A call is

- `convert`: one conversion asked through the families' adapter (`families.adapter`, as widened by family E's
  `families.e.adapter_e`), so that `d4c2dda`'s api (`convert_identifier`) and `71690ab`'s (`convert`) are read
  alike: the ids of the view asked (`best` or `all`), the status, the match, the copies reported apart, the error;
- `infer`: source inference on a list read from the truth store (the gene names of Ensembl's GTF of a release, one
  per gene of the primary assembly, as the hunt read them), `infer_identifier_source` on `d4c2dda`, `infer_source`
  on `71690ab`;
- `truth`: a fact of the truth store the finding rests on (a gene present or absent in Ensembl's GTF or in the core
  database's gene table of a release); no graph.

The recorded answer is a set of fields, each compared with the answer got (`compare`): `outcome`, `status`, `match`,
`error_type` (equal), `got` (the same ids), `got_base` (the same stable ids, versions aside, where the hunt named no
version), `got_has` (these ids among those got), `n_got` (as many ids), `apart_any` (copies reported apart, or not),
`present` (a truth fact). A finding holds when every call answers as recorded.

A run loads the configured tree's graph of each organism its calls need (`config`: `[trees.<tree>]` source, python
and graphs; the process restarts itself in the tree's Python with `PYTHONHASHSEED=0`; one graph at a time), makes
the calls, prints expected, recorded and got, and exits 0 when every finding run still holds as recorded, 1 when an
answer changed (the finding was fixed, or changed: read the lines marked CHANGED), 2 when it cannot run.
"""
import argparse
import gc
import json
import os
import re
import shutil
import sys
import tempfile
import time

from idtrack_reproducibility import config as CFG

ANGLES = ("inputs", "history", "regions")
HS, MM, SS, MMU, MFA = "homo_sapiens", "mus_musculus", "sus_scrofa", "macaca_mulatta", "macaca_fascicularis"
FIELDS = ("outcome", "status", "match", "error_type", "got", "got_base", "got_has", "n_got", "apart_any", "present")
ENS_VERSIONED = re.compile(r"^(ENS[A-Z]*[GTPE]\d{11}|LRG_\d+)\.\d+$")


def conv(org, ident, rec, expect, source, view="best", **kw):
    """A conversion: `convert(ident, **kw)` read in view `view`."""
    return dict(kind="convert", organism=org, id=ident, view=view, kwargs=kw, recorded=rec, expected=expect,
                source=source)


def infer(org, release, mode, rec, expect, source):
    """Source inference on the gene names of Ensembl's GTF of `release` (primary assembly, one per gene)."""
    return dict(kind="infer", organism=org, release=release, mode=mode, recorded=rec, expected=expect, source=source)


def fact(org, table, release, gene, present, expect, source):
    """A fact of the truth store: `gene` is (or is not) in `table` (`gtf` or `core`) of `release`."""
    return dict(kind="truth", organism=org, table=table, release=release, gene=gene, recorded=dict(present=present),
                expected=expect, source=source)


GTF115 = "Ensembl's GTF of release 115"
MISSING = "pandas' and R's spellings of a missing value"
TP53 = "TP53's gene ENSG00000141510.20"

# ---------------- the table: every finding with a reproducer (docs/hunts.md) ----------------

FINDINGS = [
    # ---- inputs ----
    dict(angle="inputs", id="IN-1", title="missing values become genes", calls=[
        conv(HS, "nan", dict(status="converted", match="synonym:case", got=["ENSG00000168356.13"]),
             "no gene (`nan` is pandas' text of a missing value)", MISSING, to_release=115),
        conv(HS, "NA", dict(status="converted", match="synonym:exact", got=["ENSG00000047597.7"]),
             "no gene (`NA` is R's missing value)", MISSING, to_release=115),
        conv(HS, "NaT", dict(status="converted", match="synonym:case", got=["ENSG00000204256.16"]),
             "no gene (`NaT` is pandas' missing time)", MISSING, to_release=115),
        conv(HS, "unknown", dict(status="converted", match="case", got=["ENSG00000119684.17"]),
             "no gene (UniProt's placeholder gene name UNKNOWN, releases 63 to 78)",
             "UniProtKB Gene Name of releases 63 to 78 (a placeholder)", to_release=115),
        conv(MM, "-", dict(status="converted", match="synonym:separator", got=["ENSMUSG00000031354.17"]),
             "no gene (a placeholder without a letter or digit)", "a string without a letter or digit",
             to_release=115),
        conv(MM, "nan", dict(status="tie", match="synonym:case",
                             got=["ENSMUSG00000034115.11", "ENSMUSG00000054191.10"]),
             "no gene", MISSING, to_release=115),
        conv(MM, "unknown", dict(status="tie", match="synonym:exact",
                                 got=["ENSMUSG00000026005.16", "ENSMUSG00000034751.18", "ENSMUSG00000052310.10",
                                      "ENSMUSG00000061313.11"]),
             "no gene", "UniProtKB Gene Name (a placeholder)", to_release=115),
        conv(SS, "N/A", dict(status="converted", match="synonym:exact", got=["ENSSSCG00000001527.5"]),
             "no gene", MISSING, to_release=115),
    ]),
    dict(angle="inputs", id="IN-2", title="mouse symbols in another letter case reach other genes", calls=[
        conv(MM, "AR", dict(status="tie", match="synonym:exact",
                            got=["ENSMUSG00000001642.19", "ENSMUSG00000029378.6", "ENSMUSG00000032274.10"]),
             "Ar ENSMUSG00000046532.9", "MGI symbol Ar; " + GTF115, to_release=115),
        conv(MM, "CCK", dict(status="converted", match="synonym:exact", got=["ENSMUSG00000041135.16"]),
             "Cck ENSMUSG00000032532.8", "MGI symbol Cck; " + GTF115, to_release=115),
        conv(MM, "AGT", dict(status="converted", match="synonym:exact", got=["ENSMUSG00000054612.5"]),
             "Agt ENSMUSG00000031980.11", "MGI symbol Agt; " + GTF115, to_release=115),
        conv(MM, "GRIN1", dict(status="converted", match="synonym:exact", got=["ENSMUSG00000069227.5"]),
             "Grin1 ENSMUSG00000026959.14", "MGI symbol Grin1; " + GTF115, to_release=115),
        conv(MM, "nos2", dict(status="converted", match="synonym:exact", got=["ENSMUSG00000051965.9"]),
             "Nos2 ENSMUSG00000020826.10", "MGI symbol Nos2; " + GTF115, to_release=115),
        conv(SS, "afp", dict(status="converted", match="synonym:exact", got=["ENSSSCG00000001233.7"]),
             "AFP ENSSSCG00000008949.5", GTF115, to_release=115),
        conv(SS, "Has2", dict(status="converted", match="synonym:exact", got=["ENSSSCG00000057079.1"]),
             "HAS2 ENSSSCG00000005992.6", GTF115, to_release=115),
    ]),
    dict(angle="inputs", id="IN-3", title="names made unique by single-cell tools reach other genes", calls=[
        conv(HS, "SRF.1", dict(status="converted", match="synonym:separator", got=["ENSG00000132122.13"]),
             "SRF ENSG00000112658.9 (a copy of SRF), or nothing", GTF115, to_release=115),
        conv(HS, "MB.1", dict(status="converted", match="synonym:separator", got=["ENSG00000105369.10"]),
             "MB ENSG00000198125.14 (a copy of MB), or nothing", GTF115, to_release=115),
        conv(HS, "AIP_1", dict(status="converted", match="synonym:separator", got=["ENSG00000138443.18"]),
             "AIP ENSG00000110711.12 (a copy of AIP), or nothing", GTF115, to_release=115),
    ]),
    dict(angle="inputs", id="IN-4", title="invisible characters and look-alikes are not read", calls=[
        conv(HS, "﻿TP53", dict(status="not_found", got=[]), TP53, GTF115, to_release=115),
        conv(HS, "TP53​", dict(status="not_found", got=[]), TP53, GTF115, to_release=115),
        conv(HS, '"TP53"', dict(status="not_found", got=[]), TP53, GTF115, to_release=115),
        conv(HS, "TP５３", dict(status="not_found", got=[]), TP53, GTF115, to_release=115),
        conv(HS, "HLA–A", dict(status="not_found", got=[]), "HLA-A's gene ENSG00000206503.16", GTF115,
             to_release=115),
    ]),
    dict(angle="inputs", id="IN-5", title="unversioned transcript ids in another letter case are not found", calls=[
        conv(HS, "enst00000016171", dict(status="not_found", got=[]),
             "ENSG00000014919.14 (the gene of ENST00000016171)", GTF115 + " (transcripts)", to_release=115),
    ]),
    dict(angle="inputs", id="IN-6", title="the graph's internal names are read as input", calls=[
        conv(HS, "synonym_id::TP53", dict(status="converted", match="exact", got=["ENSG00000067369.15"]),
             "nothing, or an input error (no database writes this string)", "the docs of idtrack's input forms",
             to_release=115),
        conv(HS, "ENSG00000141510.Void", dict(status="converted", match="exact", got=["ENSG00000141510.20"]),
             "nothing, or an input error (`Void` is no version)", GTF115 + " (versions)", to_release=115),
        conv(HS, "ENSG00000141510.void", dict(status="converted", match="case", got=["ENSG00000141510.20"]),
             "nothing, or an input error (`void` is no version)", GTF115 + " (versions)", to_release=115),
    ]),
    dict(angle="inputs", id="IN-7", title="lists of names are placed at another release or assembly", calls=[
        infer(HS, 55, "release", dict(got=["88"]), "release 55", "Ensembl's GTF of release 55 (the list)"),
        infer(HS, 75, "release", dict(got=["115"]), "release 75", "Ensembl's GTF of release 75 (the list)"),
        infer(HS, 80, "release", dict(got=["84"]), "release 80", "Ensembl's GTF of release 80 (the list)"),
        infer(HS, 105, "release", dict(got=["115"]), "release 105", "Ensembl's GTF of release 105 (the list)"),
        infer(HS, 110, "release", dict(got=["115"]), "release 110", "Ensembl's GTF of release 110 (the list)"),
        infer(HS, 75, "assembly", dict(got=["38"]), "assembly 37 (GRCh37)", "Ensembl's GTF of release 75 (the list)"),
        infer(MM, 60, "release", dict(got=["91"]), "release 60", "Ensembl's GTF of release 60 (the list)"),
        infer(MM, 75, "release", dict(got=["115"]), "release 75", "Ensembl's GTF of release 75 (the list)"),
        infer(MM, 110, "release", dict(got=["112"]), "release 110", "Ensembl's GTF of release 110 (the list)"),
        infer(MM, 67, "assembly", dict(got=["38"]), "assembly 37 (NCBIM37)", "Ensembl's GTF of release 67 (the list)"),
    ]),
    dict(angle="inputs", id="IN-8", title="a version of more than 4,300 digits raises a bare ValueError", calls=[
        conv(HS, "ENSG00000141510." + "9" * 5000, dict(outcome="error", error_type="ValueError"),
             "nothing, or an input error naming the argument", "the docs' error classes", to_release=115),
    ]),
    dict(angle="inputs", id="IN-9", title="bare numbers are read as old numeric ids of other databases", calls=[
        conv(HS, "56386", dict(status="converted", match="exact", got=["ENSG00000198000.13"]),
             "nothing (56386 is a mouse NCBI gene id, no human one)", "NCBI gene_info", to_release=115),
    ]),
    dict(angle="inputs", id="IN-10", title="names equal up to separators: one is taken silently", calls=[
        conv(MM, "AY172335_1", dict(status="converted", match="separator", got=["ENSMUSG00000064336.1"]),
             "a tie of ENSMUSG00000064336.1 (AY172335-1) and ENSMUSG00000064359.1 (AY172335.1), or nothing",
             "Ensembl's cross-references of release 115 (two names, two genes)", to_release=115),
    ]),
    dict(angle="inputs", id="IN-11", title="names with stray punctuation from Ensembl's tables", calls=[
        conv(HS, "ADORA3;", dict(status="converted", match="exact", got=["ENSG00000121933.19"]),
             "ADORA3 ENSG00000282608.2, or nothing", GTF115, to_release=115),
        conv(MM, "'Agt'", dict(status="converted", match="synonym:case", got=["ENSMUSG00000054612.5"]),
             "Agt ENSMUSG00000031980.11, or nothing", "MGI symbol Agt; " + GTF115, to_release=115),
    ]),
    # ---- history ----
    dict(angle="history", id="hh-1", title="a branch that left an id while it lived on outranks its own line", calls=[
        conv(SS, "ENSSSCG00000014532", dict(status="converted", got=["ENSSSCG00000042508.1"]),
             "ENSSSCG00000042580, ENSSSCG00000045532, ENSSSCG00000047590 (the line the id lives on)",
             "Ensembl's stable_id_event (strict walk); REST archive/id possible_replacement",
             from_release=56, to_release=98),
    ]),
    dict(angle="history", id="hh-2", title="an id that Ensembl moved to another chromosome, answered unmarked", calls=[
        conv(HS, "ENSG00000115934", dict(status="converted", got=["ENSG00000115934.11"]),
             "ENSG00000115934.11, marked as moved (chr19 protein-coding gene at 61, chr12 transcript from 62)",
             "Ensembl's stable_id_event (115934.10 -> .11 at 61 -> 62, score 0.046); GTFs of 61 and 62",
             from_release=61, to_release=115),
    ]),
    dict(angle="history", id="hh-3", title="a returned id that reuses an old version answers two ways", calls=[
        conv(MM, "ENSMUSG00000118483", dict(status="converted", got=["ENSMUSG00000118430.3"]),
             "the same answer as from release 109 (one gene throughout)",
             "Ensembl's stable_id_event; GTFs of 97 to 115", from_release=102, to_release=108),
        conv(MM, "ENSMUSG00000118483", dict(status="no_conversion", got=[]),
             "the same answer as from release 102 (one gene throughout)",
             "Ensembl's stable_id_event; GTFs of 97 to 115", from_release=109, to_release=108),
    ]),
    dict(angle="history", id="hh-4", title="bridges where Ensembl's history has no successor", calls=[
        conv(HS, "ENSG00000198477", dict(status="converted", got=["ENSG00000275004.1"]),
             "nothing (Ensembl retires the id without successor)", "Ensembl's stable_id_event; the locus rule",
             from_release=48, to_release=76),
        conv(HS, "ENSG00000198477", dict(got_has=["ENSG00000168214.17"]),
             "nothing (Ensembl retires the id without successor)", "Ensembl's stable_id_event; the locus rule",
             view="all", from_release=48, to_release=76),
    ]),
    dict(angle="history", id="hh-5", title="the truth's GTFs lack genes the core database has", calls=[
        fact(MFA, "gtf", 112, "ENSMFAG00000000631", False, "absent (the GTF of 112 lacks 47 genes)",
             "Ensembl's GTF of release 112"),
        fact(MFA, "gtf", 111, "ENSMFAG00000000631", True, "present", "Ensembl's GTF of release 111"),
        fact(MFA, "core", 112, "ENSMFAG00000000631", True, "present",
             "Ensembl's core database of release 112 (gene table)"),
        conv(MFA, "ENSMFAG00000000631", dict(status="converted", got=["ENSMFAG00000000631.2"]),
             "ENSMFAG00000000631.2 (CYP4X1, in the core database at 112)",
             "Ensembl's core database of release 112 (gene table)", from_release=106, to_release=112),
    ]),
    dict(angle="history", id="hh-6", title="the judge's place for ids Ensembl later moved; copies before 82", calls=[
        conv(HS, "ENSG00000187929", dict(status="converted", got=["ENSG00000229259.1"]),
             "ENSG00000229259.1 is the continuation of the model the id was moved to at 49 -> 50 (chr1); the "
             "locus rule at the release of origin (chr17) calls it wrong", "GTFs of 48 to 56; Havana gene "
             "OTTHUMG00000038009", from_release=48, to_release=56),
    ]),
    # ---- regions ----
    dict(angle="regions", id="R1", title="mitochondrial gene ids reused for nuclear genes at an assembly change",
         calls=[
             conv(MMU, "MT-ND1", dict(status="converted", got=["ENSMMUG00000028699.3"]),
                  "ND1 ENSMMUG00000065372.1 (the mitochondrial gene of the name)", GTF115, from_release=97),
             conv(MMU, "MT-CO1", dict(status="converted", got_base=["ENSMMUG00000028689"]),
                  "COX1 ENSMMUG00000065358.1", GTF115, from_release=97),
             conv(MMU, "MT-ND2", dict(status="converted", got_base=["ENSMMUG00000028695"]),
                  "ND2 ENSMMUG00000065356.1", GTF115, from_release=97),
             conv(MMU, "MT-ND6", dict(status="converted", got_base=["ENSMMUG00000028672"]),
                  "ND6 ENSMMUG00000065359.1", GTF115, from_release=97),
             conv(MFA, "COX1", dict(status="converted", got=["ENSMFAG00000046417.2"]),
                  "COX1 ENSMFAG00000060073.1", GTF115, from_release=102),
         ]),
    dict(angle="regions", id="R2", title="a bridge through a shared synonym beats the gene's own name or NCBI id",
         calls=[
             conv(MMU, "MT-CO2", dict(status="converted", got=["ENSMMUG00000016857.3"]),
                  "COX2 ENSMMUG00000065380.1 (the mitochondrial gene; PTGS2 shares the synonym COX2)", GTF115,
                  from_release=97),
             conv(HS, "ENSG00000212875", dict(status="converted", got=["ENSG00000095303.18"]),
                  "MT-CO1 ENSG00000198804 (its HGNC symbol), not PTGS1 (synonym COX1)", "HGNC; " + GTF115,
                  from_release=49),
             conv(HS, "ENSG00000280256", dict(status="converted", got=["ENSG00000125363.11"]),
                  "AMELY ENSG00000099721 (its own NCBI gene id 266), not AMELX", "NCBI gene_info; GTF of 76",
                  from_release=77, to_release=76),
             conv(HS, "ENSG00000279950", dict(status="converted", got=["ENSG00000100764.15"]),
                  "RPS4Y1 ENSG00000129824 (its own NCBI gene id 6192), not PSMC1 (synonym S4)",
                  "NCBI gene_info; " + GTF115, to_release=115),
             conv(HS, "ENSG00000292345", dict(status="converted", got=["ENSG00000215788.11"]),
                  "ZBED1 ENSG00000214717 (its own NCBI gene id 9189), not TNFRSF25 (synonym TRAMP)",
                  "NCBI gene_info; GTF of 109", to_release=109),
             conv(MM, "ENSMUSG00000091705", dict(status="converted", got=["ENSMUSG00000061232.5"]),
                  "H2-Q2's own line (its NCBI gene id 15013), not H2-K1 (UniProt synonym H2-K)",
                  "NCBI gene_info; GTF of 48", from_release=101, to_release=48),
         ]),
    dict(angle="regions", id="R3", title="GENCODE's _PAR_Y gene ids are not found", calls=[
        conv(HS, "ENSG00000182378.14_PAR_Y", dict(status="not_found", got=[]),
             "the Y copy of PLCXD1, ENSG00000292344 (its own id from release 110)",
             "GENCODE's id form up to release 109; " + GTF115),
        conv(HS, "ENSG00000182378_PAR_Y", dict(status="not_found", got=[]),
             "the Y copy of PLCXD1, ENSG00000292344", "GENCODE's id form up to release 109; " + GTF115),
    ]),
    dict(angle="regions", id="R4", title="names the X and Y copies of a PAR gene share tie two genes", calls=[
        conv(HS, "PLCXD1", dict(status="tie", got=["ENSG00000182378.16", "ENSG00000292344.2"]),
             "both carry the name (a tie is right by the truth); the decided rule answers the X gene with the Y "
             "copy apart", GTF115),
    ]),
    dict(angle="regions", id="R5", title="UniProt entry names, display names at 48 and 49, are not found", calls=[
        conv(HS, "COX1_HUMAN", dict(status="not_found", got=[]),
             "ENSG00000198804 (MT-CO1; COX1_HUMAN is its display name at 49)", "Ensembl's GTF of release 49",
             from_release=49),
    ]),
    dict(angle="regions", id="R6", title="a patch or haplotype copy without successor goes to a paralog", calls=[
        conv(HS, "ENSG00000268214", dict(status="converted", got_base=["ENSG00000215764"]),
             "KIR2DL4 (the reference gene of the haplotype copy), not KIR2DL2",
             "Ensembl's alt_allele groups; GTF of 55 (same name)", from_release=70, to_release=55),
        conv(HS, "ENSG00000288366", dict(status="converted", got=["ENSG00000141854.4"]),
             "the MIR1199 gene ENSG00000277805 at the copy's coordinates, not the host gene MISP3",
             "GTFs of 76 and 102 (coordinates)", from_release=102, to_release=76),
        conv(HS, "ENSG00000284979", dict(status="converted", got_base=["ENSG00000186572"]),
             "DEFB107B ENSG00000198129, not DEFB107A", "GTFs of 76 and 100 (same name)", from_release=100,
             to_release=76),
    ]),
    dict(angle="regions", id="R7", title="an RFAM family label names every member of the family", calls=[
        conv(HS, "U6", dict(status="tie", n_got=1279),
             "the 16 genes whose display name is U6 at 75, walked to 82", "Ensembl's GTF of release 75",
             from_release=75, to_release=82),
    ]),
    dict(angle="regions", id="R8", title="another gene's live synonym is read before a former name", calls=[
        conv(HS, "RNR1", dict(status="converted", match="synonym:exact", got=["ENSG00000153234.17"]),
             "MT-RNR1 ENSG00000211459 (RNR1 is NCBI Gene 4549's symbol)",
             "NCBI gene_info; Ensembl's NCBI gene names of releases 86 to 95"),
    ]),
    dict(angle="regions", id="R9", title="an LRG id answers itself, not the gene it describes", calls=[
        conv(HS, "LRG_1", dict(status="converted", got=["LRG_1.1"]),
             "COL1A1 ENSG00000108821 (the gene of the LRG record)", "the LRG record LRG_1; Ensembl's alt_allele"),
    ]),
    dict(angle="regions", id="J1", title="judge: copies on patches judged by the place of their own sequence", calls=[
        conv(MM, "ENSMUSG00000116201", dict(status="converted", got=["ENSMUSG00000073791.5"]),
             "Efcab7 ENSMUSG00000073791 (same coordinates as the patch copy): right",
             "Ensembl's assembly_exception; GRC alignments of the patch", from_release=100, to_release=68),
    ]),
    dict(angle="regions", id="J2", title="judge: the GTF of cynomolgus 112 lacks 47 genes", calls=[
        fact(MFA, "gtf", 112, "ENSMFAG00000000631", False, "absent", "Ensembl's GTF of release 112"),
        fact(MFA, "core", 112, "ENSMFAG00000000631", True, "present",
             "Ensembl's core database of release 112 (gene table)"),
        conv(MFA, "ENSMFAG00000000631", dict(status="converted", got=["ENSMFAG00000000631.2"]),
             "ENSMFAG00000000631.2: right (the core database has the gene at 112)",
             "Ensembl's core database of release 112 (gene table)", from_release=106, to_release=112),
    ]),
    dict(angle="regions", id="J3", title="judge: copies before 82 labelled primary", calls=[
        conv(HS, "KIR2DL1", dict(status="converted", got_base=["ENSG00000125498"], apart_any=True),
             "the reference gene, its copies reported apart: right",
             "Ensembl's assembly_exception of release 75 (GTFs before 82 have no patch file)", view="all",
             from_release=75, to_release=75),
    ]),
]


# ---------------- selection and description ----------------

def select(angle, fid, organism=None):
    """The findings of `angle` (or `all`) with id `fid` (or `all`), their calls cut to `organism` when given."""
    if angle not in ANGLES + ("all",):
        raise SystemExit(f"unknown angle {angle!r}: one of {', '.join(ANGLES)}, all")
    out = []
    for f in FINDINGS:
        if angle != "all" and f["angle"] != angle:
            continue
        if fid != "all" and f["id"] != fid:
            continue
        calls = [c for c in f["calls"] if organism is None or c["organism"] == organism]
        if calls:
            out.append(dict(f, calls=calls))
    if not out and fid != "all":
        known = sorted(f["id"] for f in FINDINGS if angle in ("all", f["angle"]))
        raise SystemExit(f"no reproducer {angle} {fid}" + (f" for {organism}" if organism else "")
                         + f"; known: {', '.join(known)}")
    return out


def literal(s):
    """A string as Python source, ASCII only; a long run of one character at its end as `"x" * n`."""
    n = len(s) - len(s.rstrip(s[-1])) if s else 0
    if len(s) > 60 and n > 40:
        return f"{json.dumps(s[:-n], ensure_ascii=True)} + {json.dumps(s[-1], ensure_ascii=True)} * {n}"
    return json.dumps(s, ensure_ascii=True)


def describe(call):
    """The call as a line of Python (`convert("nan", to_release=115)`)."""
    if call["kind"] == "convert":
        args = [literal(call["id"])] + [f"{k}={v!r}" for k, v in call["kwargs"].items()]
        if call["view"] != "best":
            args.append(f'strategy="{call["view"]}"')
        return f"convert({', '.join(args)})"
    if call["kind"] == "infer":
        return f'infer_source(<gene names of GTF {call["release"]}>, mode="{call["mode"]}")'
    return f'truth: {call["gene"]} in the {call["table"]} genes of release {call["release"]}'


def base_id(x):
    m = ENS_VERSIONED.match(str(x))
    return m.group(1) if m else str(x)


def compare(recorded, observed):
    """The recorded fields the observed answer does not meet (empty: it answers as recorded)."""
    out = []
    got = list(observed.get("got") or [])
    for k, want in recorded.items():
        if k not in FIELDS:
            raise ValueError(f"unknown recorded field {k!r}")
        if k == "got":
            ok = sorted(got) == sorted(want)
            have = sorted(got)
        elif k == "got_base":
            have = sorted({base_id(x) for x in got})
            ok = have == sorted(want)
        elif k == "got_has":
            have = sorted(set(want) - set(got))
            ok = not have
            have = f"missing {have}"
        elif k == "n_got":
            have = len(got)
            ok = have == want
        elif k == "apart_any":
            have = bool(observed.get("apart"))
            ok = have == want
        elif k == "outcome":
            have = observed.get("outcome") or "answer"
            ok = have == want
        else:
            have = observed.get(k)
            ok = have == want
        if not ok:
            out.append(f"{k}: recorded {want!r}, got {have!r}")
    if observed.get("outcome") not in (None, "answer") and "outcome" not in recorded and "error_type" not in recorded:
        out.append(f"outcome: {observed.get('outcome')} ({observed.get('error')})")
    return out


def evaluate(selected, answer):
    """Ask every call (`answer(finding, index, call)` returns the observed answer) and compare; one result per finding."""
    results = []
    for f in selected:
        rows = []
        for i, c in enumerate(f["calls"]):
            obs = answer(f, i, c)
            rows.append(dict(call=c, observed=obs, differences=compare(c["recorded"], obs)))
        results.append(dict(angle=f["angle"], id=f["id"], title=f["title"], calls=rows,
                            holds=all(not r["differences"] for r in rows)))
    return results


def exit_code(results):
    return 0 if all(r["holds"] for r in results) else 1


def show_answer(obs):
    if obs is None:
        return "-"
    if "present" in obs:
        return f"present={obs['present']}"
    parts = []
    if obs.get("outcome") not in (None, "answer"):
        parts.append(f"outcome={obs['outcome']} error={obs.get('error_type')}")
    for k in ("status", "match"):
        if obs.get(k):
            parts.append(f"{k}={obs[k]}")
    got = list(obs.get("got") or [])
    parts.append(f"got={got if len(got) <= 6 else got[:6] + [f'... {len(got)} in all']}")
    if obs.get("apart"):
        parts.append(f"apart={len(obs['apart'])}")
    return " ".join(parts)


def report(results, tree, out=sys.stdout):
    for r in results:
        for row in r["calls"]:
            c = row["call"]
            out.write(f"{r['id']}  {c['organism']}  {describe(c)}\n")
            out.write(f"  expected: {c['expected']} (truth: {c['source']})\n")
            out.write(f"  recorded: {c['recorded']}\n")
            out.write(f"  got:      {show_answer(row['observed'])}\n")
            out.write("  -> as recorded\n" if not row["differences"] else
                      "  -> CHANGED: " + "; ".join(row["differences"]) + "\n")
        out.write(f"{r['angle']} {r['id']} ({r['title']}): "
                  f"{'holds as recorded' if r['holds'] else 'CHANGED'} on {tree}\n\n")


# ---------------- asking a tree ----------------

def truth_fact(call):
    if call["table"] == "gtf":
        from idtrack_reproducibility.harness import harness_lib as H
        genes = H.gtf_genes(call["organism"], call["release"])
    else:
        from idtrack_reproducibility.families.astar import truth as T
        genes = T.core_gene_ids(call["organism"], call["release"])
    if genes is None:
        return dict(outcome="error", error_type="MissingTruth",
                    error=f"no {call['table']} genes of release {call['release']} in the truth store")
    return dict(present=call["gene"] in genes)


def gtf_names(org, release):
    """The gene names of Ensembl's GTF of `release`, primary assembly, one per gene (the hunt's list)."""
    from idtrack_reproducibility.harness import harness_lib as H
    genes = H.gtf_genes(org, release, "primary")
    if genes is None:
        raise SystemExit(f"no GTF genes of {org} release {release} in the truth store")
    return [r["gene_name"] for r in genes.values() if r.get("gene_name")]


class Session:
    """The configured tree's adapter, one organism at a time (a graph is let go before the next one is loaded)."""

    def __init__(self, tree, timeout, workdir):
        cfg = CFG.get()
        self.tree = cfg.tree(tree)
        if not self.tree["source"]:
            raise SystemExit(f"tree {tree}: no `source` (a checkout of idtrack) in the configuration")
        if self.tree["source"] not in sys.path:
            sys.path.insert(0, self.tree["source"])  # the tree's idtrack before any other
        self.timeout, self.workdir = timeout, workdir
        self.organism, self.ad, self.info = None, None, {}

    def adapter(self, org):
        if org != self.organism:
            self.close()
            from idtrack_reproducibility.families import common as C
            from idtrack_reproducibility.families.e import adapter_e
            graph = CFG.get().graph(self.tree["name"], org)
            if not os.path.exists(graph):
                raise SystemExit(f"tree {self.tree['name']}: graph of {org} not found: {graph}")
            self.ad = adapter_e.AdapterE(self.tree["source"], graph, org, C.release_of_graph(graph),
                                         os.path.join(self.workdir, org), timeout=self.timeout)
            self.organism = org
            self.info[org] = dict(self.ad.info)
        return self.ad

    def close(self):
        if self.ad is not None:
            self.ad = None
            gc.unfreeze()
            gc.collect()
            self.organism = None

    def answer(self, finding, i, call):
        if call["kind"] == "truth":
            return truth_fact(call)
        ad = self.adapter(call["organism"])
        qid = f"hunt|{finding['angle']}|{finding['id']}|{i}"
        if call["kind"] == "infer":
            names = gtf_names(call["organism"], call["release"])
            line = dict(qid=qid, id=f"names of GTF {call['release']}", call="infer_source", batch=names,
                        kwargs=dict(mode=call["mode"]), strategy="best")
        else:
            kw = call["kwargs"]
            line = dict(qid=qid, id=call["id"], strategy=call["view"], from_release=kw.get("from_release"),
                        to_release=kw.get("to_release"), final_database=kw.get("final_database"),
                        from_database=kw.get("from_database"))
        a = ad.ask(line)[0]
        obs = dict(outcome=a.get("outcome"), status=a.get("status"), match=a.get("match"), got=list(a.get("got") or []),
                   apart=list(a.get("apart") or []), error_type=a.get("error_type"), error=a.get("error"),
                   seconds=a.get("seconds"))
        if call["kind"] == "infer":
            obs["n_items"] = len(line["batch"])
        return obs


def restart_in_tree_python(tree, argv):
    """Run this command (`argv`, its arguments) again in the tree's Python with PYTHONHASHSEED=0 and this package on
    PYTHONPATH."""
    if os.environ.get("IDTRACK_REPRO_HUNTS_RESTARTED") == "1":
        return
    want = CFG.get().tree(tree)["python"]
    exe = shutil.which(want) or want
    if os.path.realpath(exe) == os.path.realpath(sys.executable) and os.environ.get("PYTHONHASHSEED") == "0":
        return
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1", IDTRACK_REPRO_HUNTS_RESTARTED="1")
    env["PYTHONPATH"] = os.pathsep.join(x for x in (CFG.SRC, os.environ.get("PYTHONPATH")) if x)
    sys.stdout.flush()
    os.execve(exe, [exe, "-m", "idtrack_reproducibility.hunts.reproduce"] + list(argv), env)


def run(a, argv):
    selected = select(a.angle, a.id, a.organism)
    tree = a.tree
    if any(c["kind"] != "truth" for f in selected for c in f["calls"]):
        restart_in_tree_python(tree, argv)
    work = tempfile.mkdtemp(prefix="idtrack_hunts_")
    try:
        sess = Session(tree, a.timeout, work)
        t0 = time.time()
        results = evaluate(selected, sess.answer)
        sess.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)  # only links to the graphs live there
    report(results, tree)
    code = exit_code(results)
    print(f"{sum(r['holds'] for r in results)} of {len(results)} findings hold as recorded on {tree} "
          f"({time.time() - t0:.0f} s)")
    if a.json:
        cfg = CFG.get()
        os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
        doc = dict(tree=tree, commit=cfg.tree(tree)["commit"], graphs=sess.info, when=time.strftime("%Y-%m-%dT%H:%M:%S"),
                   python=sys.version.split()[0], pythonhashseed=os.environ.get("PYTHONHASHSEED"), exit=code,
                   results=[dict(r, calls=[dict(call=dict(row["call"], call=describe(row["call"])), observed=row["observed"],
                                                differences=row["differences"]) for row in r["calls"]])
                            for r in results])
        with open(a.json + ".tmp", "w") as fh:
            json.dump(doc, fh, indent=1, default=str)
        os.replace(a.json + ".tmp", a.json)
    return code


def submit(a):
    """One job per organism of the selected calls; in each, the reproducers of that organism on every tree asked."""
    from idtrack_reproducibility import jobs
    cfg = CFG.get()
    selected = select(a.angle, a.id, a.organism)
    trees = a.tree.split(",")
    for t in trees:
        cfg.tree(t)
    out = os.path.abspath(a.out or os.path.join(cfg.hunts, "reproduce"))
    orgs = sorted({c["organism"] for f in selected for c in f["calls"]})
    ids = []
    for org in orgs:
        cmds = []
        for t in trees:
            js = os.path.join(out, t, f"{a.angle}_{a.id}_{org}.json")
            cmds.append(jobs.module(cfg.python, "idtrack_reproducibility.hunts.reproduce",
                                    f"{a.angle} {a.id} --organism {org} --tree {t} --timeout {a.timeout} --json {js}"))
        what = re.sub(r"[^a-z0-9_]", "_", f"hunts_{a.angle}_{a.id}_{org}".lower())
        ids.append(jobs.submit(what, " ; ".join(cmds), os.path.join(out, "logs"), cpus=2,
                               mem_gb=jobs.memory_for(org), time=a.time))
        print(org, ids[-1])
    return ids


def listing():
    for f in FINDINGS:
        print(f"{f['angle']:8} {f['id']:6} {f['title']}")
        for c in f["calls"]:
            print(f"    {c['organism']:20} {describe(c)}")
            print(f"    {'':20} expected: {c['expected']} (truth: {c['source']}); recorded: {c['recorded']}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.hunts.reproduce",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("angle", nargs="?", help="inputs, history, regions, or all")
    ap.add_argument("id", nargs="?", help="a finding's id (IN-1, hh-1, R1, J1, ...) or all")
    ap.add_argument("--tree", default="d4c2dda", help="a configured tree (with --submit: a comma-separated list)")
    ap.add_argument("--organism", default=None, help="only the calls of this organism")
    ap.add_argument("--json", default=None, help="write the calls, answers and differences to this file")
    ap.add_argument("--timeout", type=float, default=300.0, help="seconds per call (the adapter's budget)")
    ap.add_argument("--list", action="store_true", help="print the table and exit")
    ap.add_argument("--submit", action="store_true", help="run as jobs, one per organism (jobs.submit)")
    ap.add_argument("--out", default=None, help="with --submit: folder of the answers and logs (default <hunts>/reproduce)")
    ap.add_argument("--time", default="04:00:00", help="with --submit: the jobs' time limit")
    argv = sys.argv[1:] if argv is None else list(argv)
    a = ap.parse_args(argv)
    if a.list:
        listing()
        return 0
    if not a.angle or not a.id:
        ap.error("give an angle and a finding id (or all), or --list")
    if a.submit:
        submit(a)
        return 0
    if "," in a.tree:
        ap.error("one tree per run (several only with --submit)")
    try:
        return run(a, argv)
    except SystemExit as e:  # a message: the run could not be made (no tree, graph or truth file)
        if isinstance(e.code, str):
            print(e.code, file=sys.stderr)
            return 2
        raise


if __name__ == "__main__":
    sys.exit(main())
