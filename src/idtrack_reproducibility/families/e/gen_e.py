"""Family E's query sets E1 to E5, E8's static relations and the input hunt's reproducers, per organism.

Usage: gen_e.py --organism <o> [--sets e1,e2,...] [--dry]      (idtrack_dev_env; reads only the truth store and
                                                               the harness's published sets)

Every set is published twice over (docs/families/f0.md, "Query sets"): the base set judges by **71690ab's
docs** (judge `now`: the docstrings of `API.convert`, `convert_many`, `infer_source` and the error hierarchy, as
built); its revision set `<set>_dec` holds the lines whose expectation changes under **the decided rules not built yet**
(judge `decided`: rules `missing_values`, `no_synonym_target`, `inference_per_release`, `single_case_symbols`, `made_unique_copies`, `stray_punctuation`, `unicode_normalised`, `bounded_arguments`, `bare_numbers` and `one_case_schemes` of docs/rules.md, marked decided): `--sets <set>` for the first,
`--sets <set>,<set>_dec` for the second (score_e.py does both).

Strata (docs/methods.md, E; the input hunt's findings IN-1 to IN-11 as strata, docs/hunts.md):

- e1 missing values (IN-1, rule `missing_values`): every spelling of rule `missing_values` in several letter cases, with white
  space, and the hunt's further spellings, alone, with `from_database` named, and in batches at every position.
- e2 malformed strings (IN-3, -4, -5, -6, -8, -10, -11; rules `made_unique_copies`, `stray_punctuation`, `unicode_normalised`, `bounded_arguments` and `one_case_schemes`): forms of real identifiers (names,
  gene, transcript and translation ids, NCBI, RefSeq, UniProt; the harness's published lines give the expected
  gene) with invisible and control characters, look-alikes, quotes and stray punctuation, make-unique suffixes,
  versions of every length, many separators, internal names; very long strings, injection strings, numbers.
- e2_case (IN-2, rule `single_case_symbols`): every MGI symbol of mouse in capitals and in lower case.
- e3 types (rule `typed_api`): non-strings and every container of a batch.
- e4 another organism: every current gene id of every other organism (versioned), samples of their unversioned
  ids, transcripts, translations, symbols and NCBI gene ids (IN-9), asked in this organism's graph.
- e5 arguments at and beyond their bounds: releases, databases, budgets, switches, modes, the other api calls.
- e8 static relations: identity (every gene id at every release it exists at, to that release), the versioned and
  unversioned forms of a current id, a unique name and its case variants, white space, and `infer_source` under
  shuffling and duplicates (the rule `relation`).
- hunt_inputs: the reproducers of findings IN-1 to IN-11 (and IN-7's source inference of name lists).

Expectations come from Ensembl's tables (truth_e.py; the harness's published lines) and from the docs; never
from idtrack. A form that is itself an identifier the truth knows (any table, any release, in any case and
separator form) is `open` (the truth cannot say which the user meant). Line fields of family E: adapter_e.py;
`e_rule` (`contains`: the expected values are among the answer's; `or_error`: an error of the classes in
`or_error` is right as well as the line's rule) is applied by score_e.py after F0's scorer.
"""
import argparse
import collections
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.e import truth_e as T

FAMILY = "e"
SEED = 20261118
NAME_DB = {"homo_sapiens": "HGNC Symbol", "mus_musculus": "MGI Symbol", "sus_scrofa": "VGNC Symbol",
           "macaca_mulatta": "VGNC Symbol", "macaca_fascicularis": "VGNC Symbol"}
WINDOW = {"homo_sapiens": (48, 115), "mus_musculus": (48, 115), "sus_scrofa": (56, 115), "macaca_mulatta": (48, 115),
          "macaca_fascicularis": (91, 115)}
SNAP = 115

# error classes: 71690ab's class | d4c2dda's built-in (adapter_e tags d4c2dda's classes `d4c2dda:<class>`)
E_TYPE = "ArgumentTypeError|d4c2dda:TypeError"
E_MISSING = "MissingValueError|d4c2dda:ValueError"
E_RELEASE = "UnknownReleaseError|d4c2dda:ValueError"
E_DB = "UnknownDatabaseError|d4c2dda:ValueError"
E_INVALID = "InvalidArgumentError|d4c2dda:ValueError"
E_UNKNOWN_ID = "UnknownIdentifierError|d4c2dda:ValueError"
E_INPUT = "InputError"


def enc(v):
    """An encoded value for adapter_e.decode (plain JSON where it is one)."""
    if isinstance(v, float) and (v != v or v in (float("inf"), float("-inf"))):
        return {"$": "float", "v": str(v)}
    return v


class Gen:
    def __init__(self, org):
        self.org = org
        self.rng = random.Random(f"{SEED}:{org}")
        self.names = T.Names(org)
        self.cur = T.Current(org)
        self.sets = collections.defaultdict(list)   # set -> lines (judge now)
        self.dec = collections.defaultdict(dict)     # set -> qid -> changes (judge decided)
        self.seen = set()
        self._harness = None

    # ---------------- populations from the harness's published lines (truth; never idtrack) ----------------

    def harness(self, kind, k, pred=None):
        if self._harness is None:
            right = self.reference_right()
            lines, _, _ = L.load("harness", self.org)
            lines = [q for q in lines if q["qid"] in right]
            by = collections.defaultdict(list)
            for q in lines:
                if q.get("from_release") is None and q.get("final_database") is None and q.get("to_release") == SNAP \
                        and (q.get("strategy") or "best") == "best" and q.get("rule") == "set" and q.get("expect") \
                        and not q.get("open") and not q.get("apart") and not q.get("broad"):
                    by[q["kind"]].append(q)
            for v in by.values():
                v.sort(key=lambda q: q["qid"])
            self._harness = by
        pop = [q for q in self._harness.get(kind, []) if pred is None or pred(q)]
        r = random.Random(f"{SEED}:{self.org}:{kind}")
        return r.sample(pop, min(k, len(pop)))

    def reference_right(self):
        """qids of the harness's lines that the reference d4c2dda answers right (F0's stored check run): the bases of
        the forms, so that a wrong answer of a form is the form's own, not its base's (docs/families/e.md)."""
        import glob
        import pyarrow.parquet as pq
        d = os.path.join(C.BENCH, "f0", "runs", "check_d4c2dda", self.org, "outcomes")
        ok = set()
        for p in glob.glob(os.path.join(d, "*.parquet")):
            t = pq.read_table(p, columns=["qid", "view", "outcome"]).to_pandas()
            ok |= set(t.loc[(t["outcome"] == "right") & (t["view"] == "best"), "qid"])
        if not ok:
            raise SystemExit(f"no stored reference outcomes in {d}")
        return ok

    # ---------------- lines ----------------

    def add(self, set_, stratum, key, ident, expect=(), rule="set", dec=None, unit=None, **kw):
        qid = f"e|{stratum}|{key}"
        if qid in self.seen:
            return None
        self.seen.add(qid)
        line = dict(qid=qid, kind=set_, stratum=stratum, unit=unit or C.base_id(str(ident))[:80], id=ident,
                    to_release=kw.pop("to_release", SNAP), from_release=kw.pop("from_release", None),
                    final_database=kw.pop("final_database", None), strategy=kw.pop("strategy", "best"),
                    rule=rule, cmp=kw.pop("cmp", "versioned"), expect=list(expect), broad=None, locus=None)
        line.update({k: v for k, v in kw.items() if v is not None})
        self.sets[set_].append(line)
        if dec:
            self.dec[set_][qid] = dec
        return line

    def form_line(self, set_, stratum, base, form, how, now, dec=None, **kw):
        """A form of a base line: judged by `now` ("set": found; "safe": the gene or nothing; "nothing") unless the
        form is itself an identifier the truth knows (then open)."""
        E = list(base["expect"])
        key = f"{how}|{base['id']}"
        known = self.names.known(form) and not (self.names.owners_folded(form) <= {C.base_id(e) for e in E})

        def rule_of(j):
            if j == "set":
                return dict(rule="set", expect=E)
            if j == "safe":
                return dict(rule="safe", expect=E)
            if j == "nothing":
                return dict(rule="set", expect=[])
            if isinstance(j, dict):
                return j
            raise ValueError(j)
        r = rule_of(now)
        extra = {}
        if known:
            extra["open"] = "the form is also an identifier of another gene in Ensembl's tables"
        d = None
        if dec is not None and not known:
            d = rule_of(dec)
            if d == r:
                d = None
        return self.add(set_, stratum, key, form, unit=C.base_id(base["id"]), dec=d, note=how, base=base["id"],
                        **r, **extra, **kw)

    # ---------------- E1 missing values ----------------

    def e1(self):
        s = "e1"
        core = ["nan", "NaN", "NAN", "Nan", "na", "NA", "Na", "nA", "n/a", "N/A", "N/a", "none", "None", "NONE",
                "null", "Null", "NULL"]
        more = ["NaT", "nat", "NAT", "<NA>", "#N/A", "#NA", "N.A.", "-", ".", "?", "--", "unknown", "Unknown",
                "UNKNOWN", "nil", "NIL", "missing", "undefined", "NA_character_", "None.", "#VALUE!", "#REF!",
                "#DIV/0!", "inf", "-inf", "Infinity"]
        blanks = ["", " ", "\t", "\n", "\r\n", "  \t ", " ", "　", "  "]
        # rule `missing_values`'s spellings and rule `missing_values`'s additions (NaT, the placeholders); the rest are open both ways
        d52 = {x.casefold() for x in ["nan", "na", "n/a", "none", "null"]}
        d81 = {"nat", "unknown"}

        def judged(tok):
            t = tok.strip().casefold()
            if not tok.strip():
                return "blank"
            if t in d52:
                return "d52"
            if t in d81:
                return "d81"
            return "other"
        toks = []
        for t in core + more:
            toks += [t, f" {t}", f"{t}\t", f" {t}\n"] if t in core[:4] + ["NA", "N/A", "None", "null", "NaT", "unknown"] else [t]
        toks += blanks
        names = self.harness("name_unique_to_gene", 2)
        ids = self.harness("gene_versioned_to_self", 1)
        if not names or not ids:
            return
        v1, v2 = names[0], ids[0]
        for tok in dict.fromkeys(toks):
            j = judged(tok)
            if j == "blank":
                now = dict(rule="set", expect=[], expect_error=E_MISSING)
                dec = None
            elif j in ("d52", "d81"):
                now = dict(rule="set", expect=[], open="71690ab reads it as a name (decision %s not built)" %
                                                         ("52" if j == "d52" else "81"))
                dec = dict(open=None, expect_error=E_MISSING)
            else:
                now = dict(rule="set", expect=[], open="a spelling of a missing value outside rule `missing_values`")
                dec = None
            tag = json.dumps(tok)
            self.add(s, f"e1_single_{j}", f"single|{tag}", tok, unit=tag, dec=dec, **now)
            # with the source database named: looked up (rule `missing_values`); no missing-value error either way
            for db in ("UniProtKB Gene Name", NAME_DB[self.org]):
                if j == "blank":
                    self.add(s, "e1_fromdb_blank", f"fromdb|{db}|{tag}", tok, unit=tag, from_database=db,
                             **dict(rule="set", expect=[], expect_error=E_MISSING))
                else:
                    self.add(s, f"e1_fromdb_{j}", f"fromdb|{db}|{tag}", tok, unit=tag, from_database=db,
                             rule="set", expect=[], open="looked up in the named database (rule `missing_values`)")
            # in a batch at every position
            for pos in range(3):
                batch = [v1["id"], v2["id"]]
                batch.insert(pos, tok)
                if j == "blank":
                    now = dict(rule="set", expect=[], expect_error=E_MISSING)
                    dec = None
                elif j in ("d52", "d81"):
                    now = dict(rule="set", expect=[], open="71690ab reads it as a name (decision %s not built)" %
                                                             ("52" if j == "d52" else "81"))
                    dec = dict(open=None, expect_error=E_MISSING)
                else:
                    now = dict(rule="set", expect=[], open="a spelling of a missing value outside rule `missing_values`")
                    dec = None
                self.add(s, f"e1_batch_{j}", f"batch|{pos}|{tag}", tok, unit=tag, call="convert_many", batch=batch,
                         pos=pos, dec=dec, **now)
                # the other items of that batch keep their answer (E8: an added item changes no other)
                if j != "blank":
                    for k, v in ((0 if pos else 1, None), (2 if pos != 2 else 1, None)):
                        it = batch[k]
                        exp = v1["expect"] if it == v1["id"] else v2["expect"]
                        self.add(s, "e1_batch_neighbour", f"batchnb|{pos}|{k}|{tag}", it, unit=tag,
                                 call="convert_many", batch=batch, pos=k, rule="set", expect=exp)

    # ---------------- E2 malformed strings ----------------

    def e2(self):
        s = "e2"
        bases = dict(name=self.harness("name_unique_to_gene", 400), gv=self.harness("gene_versioned_to_self", 150),
                     gu=self.harness("gene_unversioned_to_self", 150), tx=self.harness("transcript_to_gene", 150),
                     tl=self.harness("translation_to_gene", 150), ncbi=self.harness("ncbi_to_gene", 200),
                     refseq=self.harness("refseq_to_gene", 100), uniprot=self.harness("uniprot_to_gene", 100),
                     syn=self.harness("synonym_to_gene", 150))
        CYR = dict(zip("ABCEHKMOPTXaceopxy", "АВСЕНКМОРТХасеорху"))
        GRK = dict(zip("ABEHIKMNOPTXYZ", "ΑΒΕΗΙΚΜΝΟΡΤΧΥΖ"))

        def full(x):
            return "".join(chr(ord(c) + 0xFEE0) if 0x21 <= ord(c) <= 0x7E else c for c in x)

        def homoglyph(x, table):
            for i, c in enumerate(x):
                if c in table:
                    return x[:i] + table[c] + x[i + 1:]
            return None
        for kind, pop in bases.items():
            for b in pop:
                x = b["id"]
                mid = max(1, len(x) // 2)
                # white space: stripped (71690ab strips str.isspace): found both ways
                for how, f in (("ws_lead", " " + x), ("ws_trail", x + " "), ("ws_tabnl", "\t" + x + "\r\n"),
                               ("ws_nbsp", " " + x + " "), ("ws_thin", " " + x), ("ws_ideo", "　" + x)):
                    self.form_line(s, f"e2_ws_{kind}", b, f, how, "set")
                # zero-width characters and BOM: rule `unicode_normalised` removes them
                for how, f in (("zwsp_trail", x + "​"), ("zwsp_lead", "​" + x), ("zwsp_inner", x[:mid] + "​" + x[mid:]),
                               ("bom_lead", "﻿" + x), ("wj_trail", x + "⁠"), ("zwj_trail", x + "‍"),
                               ("zwnj_inner", x[:mid] + "‌" + x[mid:])):
                    self.form_line(s, f"e2_zw_{kind}", b, f, how, "safe", dec="set")
                # other format characters (bidi marks, soft hyphen): no decision reads them
                for how, f in (("lrm_lead", "‎" + x), ("rlm_trail", x + "‏"), ("rlo", "‮" + x + "‬"),
                               ("lri", "⁦" + x + "⁩"), ("shy_inner", x[:mid] + "­" + x[mid:])):
                    self.form_line(s, f"e2_bidi_{kind}", b, f, how, "safe")
                # control characters
                for how, f in (("nul_trail", x + "\x00"), ("nul_lead", "\x00" + x), ("nul_inner", x[:mid] + "\x00" + x[mid:]),
                               ("soh", x + "\x01"), ("ansi", "\x1b[31m" + x + "\x1b[0m"), ("del", x + "\x7f"),
                               ("bel", "\x07" + x)):
                    self.form_line(s, f"e2_ctrl_{kind}", b, f, how, "safe")
                # full-width letters and digits: rule `unicode_normalised` (NFKC)
                self.form_line(s, f"e2_fullwidth_{kind}", b, full(x), "fullwidth", "safe", dec="set")
                # Unicode dashes for '-': rule `unicode_normalised`
                if "-" in x:
                    for cp in (0x2010, 0x2011, 0x2012, 0x2013, 0x2014, 0x2212, 0xFE63, 0xFF0D):
                        self.form_line(s, f"e2_dash_{kind}", b, x.replace("-", chr(cp)), f"dash_{cp:04x}", "safe", dec="set")
                # look-alike letters of other scripts and emoji: never read
                for how, tab in (("cyrillic", CYR), ("greek", GRK)):
                    f = homoglyph(x, tab)
                    if f:
                        self.form_line(s, f"e2_lookalike_{kind}", b, f, how, "safe")
                for how, f in (("emoji_trail", x + "\U0001f9ec"), ("emoji_lead", "\U0001f600" + x)):
                    self.form_line(s, f"e2_emoji_{kind}", b, f, how, "safe")
                # quotes and stray punctuation (IN-11): rule `stray_punctuation` reads them away
                for how, f in (("q_single", f"'{x}'"), ("q_double", f'"{x}"'), ("semicolon", x + ";"), ("comma", x + ","),
                               ("q_curly", f"“{x}”"), ("q_curly1", f"‘{x}’")):
                    self.form_line(s, f"e2_punct_{kind}", b, f, how, "safe", dec="set")
                if kind not in ("gv", "refseq", "uniprot"):  # `X.` of a versioned id or an accession is a version form
                    self.form_line(s, f"e2_punct_{kind}", b, x + ".", "dot", "safe", dec="set")
                for how, f in (("paren", f"({x})"), ("bracket", f"[{x}]"), ("angle", f"<{x}>"), ("backtick", f"`{x}`")):
                    self.form_line(s, f"e2_brackets_{kind}", b, f, how, "safe")
                # inner space; prefixes of databases (IN-4)
                self.form_line(s, f"e2_innerspace_{kind}", b, x[:mid] + " " + x[mid:], "inner_space", "safe")
                if kind in ("name", "syn"):
                    # make-unique suffixes (IN-3, rule `made_unique_copies`): never another gene
                    for how, f in (("mu_dash1", x + "-1"), ("mu_dot1", x + ".1"), ("mu_us1", x + "_1"), ("mu_dash2", x + "-2"),
                                   ("mu_dot2", x + ".2")):
                        self.form_line(s, f"e2_makeunique_{kind}", b, f, how, "safe")
                    # internal names (IN-6, rule `no_synonym_target`): refused
                    self.form_line(s, f"e2_internal_{kind}", b, "synonym_id::" + x, "synonym_prefix", "safe",
                                   dec=dict(rule="set", expect=[], e_rule="or_error", or_error=E_INPUT))
                    # many separators: a separator flipped is read (fold equal), doubled is not
                    if sum(x.count(c) for c in "-_.") >= 1:
                        for how, f in (("sep_us", x.replace("-", "_").replace(".", "_")),
                                       ("sep_dot", x.replace("-", ".").replace("_", ".")),
                                       ("sep_dash", x.replace("_", "-").replace(".", "-"))):
                            if f != x:
                                self.form_line(s, f"e2_separators_{kind}", b, f, how, "set")
                        self.form_line(s, f"e2_separators_{kind}", b, x.replace("-", "--", 1) if "-" in x else x + "--",
                                       "sep_double", "safe")
                if kind in ("gv", "gu"):
                    g = C.base_id(x)
                    vv = x.rsplit(".", 1)[1] if "." in x else None
                    for how, f in (("lower", x.lower()), ("title", x[:1] + x[1:].lower())):
                        self.form_line(s, f"e2_idcase_{kind}", b, f, how, "set")
                    self.form_line(s, f"e2_innerspace_{kind}", b, x[:4] + " " + x[4:], "space_after_prefix", "safe")
                    self.form_line(s, "e2_prefix", b, "ensembl:" + x, "prefix_ensembl", "safe")
                    if vv:
                        for how, f in (("v_dash", f"{g}-{vv}"), ("v_us", f"{g}_{vv}"), ("v_space", f"{g} .{vv}"),
                                       ("v_letter", f"{g}.v{vv}"), ("v_zero", f"{g}.0{vv}")):
                            self.form_line(s, f"e2_versionform_{kind}", b, f, how, "safe")
                        # versions of every length (IN-8, rule `bounded_arguments`; rule `version_newer_than_graph` reads a newer version as the gene)
                        for k in (3, 10, 100, 1000, 4300, 4301, 5000, 100000):
                            f = f"{g}." + "9" * k
                            if k <= 3:
                                self.form_line(s, "e2_longversion", b, f, f"version_{k}", "set")
                            elif k == 10:
                                self.form_line(s, "e2_longversion", b, f, f"version_{k}", "set",
                                               dec=dict(rule="set", expect=b["expect"], e_rule="or_error", or_error=E_INPUT))
                            else:
                                self.form_line(s, "e2_longversion", b, f, f"version_{k}", "safe",
                                               dec=dict(rule="set", expect=[], expect_error=E_INPUT))
                if kind in ("tx", "tl"):
                    unv = "." not in x
                    for how, f in (("lower", x.lower()), ("title", x[:1] + x[1:].lower())):
                        if unv:  # IN-5: unversioned transcript and translation ids in another case (rule `one_case_schemes`)
                            self.form_line(s, f"e2_txcase_{kind}_unversioned", b, f, how, "safe", dec="set")
                        else:
                            self.form_line(s, f"e2_txcase_{kind}_versioned", b, f, how, "set")
                if kind == "ncbi":
                    for how, f in (("zero1", "0" + x), ("zero2", "00" + x), ("float", x + ".0"), ("plus", "+" + x),
                                   ("geneid", "GeneID:" + x), ("ncbigene", "NCBIGene:" + x), ("sci", f"{int(x):.6e}" if x.isdigit() else None),
                                   ("sci_upper", f"{int(x):.2E}" if x.isdigit() else None)):
                        if f:
                            self.form_line(s, "e2_numbers_ncbi", b, f, how, "safe")
                if kind in ("refseq", "uniprot"):
                    self.form_line(s, f"e2_acc_{kind}", b, x.lower(), "lower", "safe")
                    self.form_line(s, f"e2_acc_{kind}", b, x.replace("_", ""), "no_underscore", "safe")
                    self.form_line(s, f"e2_prefix", b, ("RefSeq:" if kind == "refseq" else "UniProtKB:") + x, "prefix", "safe")
        # strings that are no identifier at all: nothing (very long ones refused by rule `bounded_arguments`)
        name = bases["name"][0]["id"] if bases["name"] else "TP53"
        junk = {
            "inj_sql1": "' OR 1=1 --", "inj_sql2": "\"; DROP TABLE genes; --", "inj_sql3": f"{name}' OR '1'='1",
            "inj_html": "<script>alert(1)</script>", "inj_html2": "<img src=x onerror=alert(1)>", "inj_path": "../../etc/passwd",
            "inj_path2": "..\\..\\windows\\win.ini", "inj_fmt": "%s%s%s%n", "inj_fmt2": "{0}{1}", "inj_fmt3": "%(x)s",
            "inj_jndi": "${jndi:ldap://x/a}", "inj_shell": "$(rm -rf /)", "inj_tick": "`id`", "inj_url": "http://x.org/?a=1&b=2",
            "inj_json": "{\"identifier\": \"TP53\"}", "inj_null": "\\x00", "inj_crlf": "a\r\nSet-Cookie: x=1",
            "wild_star": "*", "wild_pct": "%", "wild_us": "_", "wild_re": ".*", "wild_re2": f"^{name}$", "wild_q": name[:-1] + "?",
            "wild_glob": name[:-1] + "*", "emoji_only": "\U0001f9ec\U0001f9ec", "rtl_only": "‮‭", "zw_only": "​​",
            "mixed_scripts": "ТР53", "arabic": "جين", "cjk": "基因", "combining": "é́́", "surrogates": "\ud83e",
            "num_sci1": "1e5", "num_sci2": "1.00E+05", "num_neg": "-1", "num_pi": "3.14159", "num_hex": "0x1BF5",
            "num_thousands": "7,157", "num_space": "7 157", "num_inf": "Infinity", "num_roman": "Ⅶ",
            "sep_many": "-".join("ABCDEFGHIJKLMNOP"), "sep_many2": "_".join("abcdefghijklmnopqrstuvwxyz"),
            "sep_only": "-_.-_.-_.", "sep_mixed": "A-B_C.D-E_F.G-H_I.J-K_L.M-N_O.P",
        }
        for how, f in junk.items():
            if how.startswith("num_") and self.names.known(f):
                continue
            self.add(s, "e2_junk_" + how.split("_")[0], f"junk|{how}", f, unit=how, rule="set", expect=[],
                     open="an identifier the truth knows" if self.names.known(f) else None, note=how)
        # very long strings (rule `bounded_arguments`: longer than any real identifier is refused)
        longs = {"a1000": "A" * 1000, "a100000": "A" * 100_000, "a1000000": "A" * 1_000_000,
                 "name_x100000": name + "X" * 100_000, "name_rep1000": name * 1000, "dash100000": "-" * 100_000,
                 "sep200000": "a-" * 100_000, "ens_zero100000": "ENSG" + "0" * 100_000}
        for how, f in longs.items():
            if len(f) >= 10_000:
                dec = dict(rule="set", expect=[], expect_error=E_INPUT)
            else:
                dec = dict(rule="set", expect=[], e_rule="or_error", or_error=E_INPUT)
            self.add(s, "e2_long", f"long|{how}", f, unit=how, rule="set", expect=[], dec=dec, note=how)
        if bases["name"]:
            b = bases["name"][0]
            for how, f in (("pad_ws100000", " " * 100_000 + b["id"]), ("pad_ws1000", " " * 1000 + b["id"] + "\t" * 1000)):
                dec = dict(rule="set", expect=b["expect"], e_rule="or_error", or_error=E_INPUT)
                self.add(s, "e2_long_padded", f"long|{how}", f, unit=how, rule="set", expect=b["expect"], dec=dec, note=how)

    # ---------------- E2: mouse symbols in capitals and lower case (IN-2, rule `single_case_symbols`) ----------------

    def e2_case(self):
        if self.org != "mus_musculus":
            return
        s = "e2_case"
        mgi = H.biomart_pairs(self.org, SNAP, "mgi", 1)
        if mgi is None:
            return
        g2v, v2g = mgi
        for sym in sorted(v2g):
            genes = v2g[sym]
            if len(genes) != 1:
                continue
            (g,) = genes
            if g not in self.cur.versioned:
                continue
            gv = self.cur.versioned[g]
            for how, f in (("upper", sym.upper()), ("lower", sym.lower())):
                if f == sym:
                    continue
                others = {x for x in self.names.exact115.get(f, set()) if C.base_id(x) != g}
                folded_others = {x for x in self.names.owners_folded(f) if x != g}
                if others:
                    # the docs' rule (an exact synonym before a case match): the other gene(s); rule `single_case_symbols`/80: the symbol
                    self.add(s, f"e2_case_{how}_competing", f"{how}|{sym}", f, unit=sym, rule="set", expect=sorted(others),
                             dec=dict(expect=[gv]), note=f"{how} of the MGI symbol {sym}; another gene carries it exactly")
                elif folded_others:
                    self.add(s, f"e2_case_{how}_folded", f"{how}|{sym}", f, unit=sym, rule="set", expect=[gv],
                             open="another gene carries the string in another case or separator form", note=sym)
                else:
                    self.add(s, f"e2_case_{how}_alone", f"{how}|{sym}", f, unit=sym, rule="set", expect=[gv], note=sym)

    # ---------------- E2 over every name: the classes of IN-3 and IN-11 measured whole ----------------

    def e2_all(self):
        """Every display name of one current gene at 115 (BioMart), with a trailing `;` or `.`, in single quotes (IN-11,
        rule `stray_punctuation`), and with a make-unique suffix `-1`, `.1`, `_1` (IN-3, rule `made_unique_copies`). `now`: the gene or nothing;
        `decided`: the gene for the punctuation forms, the gene or nothing for the suffixes. A form that is itself an
        identifier of another gene in Ensembl's tables is open."""
        s = "e2_all"
        for name, genes in sorted(self.names.display115.items()):
            if len(genes) != 1:
                continue
            (gv,) = genes
            if C.base_id(gv) not in self.cur.genes:
                continue
            b = dict(id=name, expect=[gv])
            for how, f in (("semicolon", name + ";"), ("dot", name + "."), ("q_single", f"'{name}'")):
                self.form_line(s, "e2_all_punct", b, f, how, "safe", dec="set")
            for how, f in (("mu_dash1", name + "-1"), ("mu_dot1", name + ".1"), ("mu_us1", name + "_1")):
                self.form_line(s, "e2_all_makeunique", b, f, how, "safe")

    # ---------------- E3 types ----------------

    def e3(self):
        s = "e3"
        names = self.harness("name_unique_to_gene", 3)
        ids = self.harness("gene_versioned_to_self", 2)
        ncbi = self.harness("ncbi_to_gene", 1)
        if len(names) < 3 or len(ids) < 2:
            return
        n0 = names[0]
        num = int(ncbi[0]["id"]) if ncbi and ncbi[0]["id"].isdigit() else 7157
        bad = {"none": None, "int": num, "int0": 0, "float": float(num), "nan": {"$": "float", "v": "nan"},
               "inf": {"$": "float", "v": "inf"}, "true": True, "false": False, "bytes": {"$": "bytes_hex", "v": n0["id"].encode().hex()},
               "list": [n0["id"]], "dict": {"$": "dict", "v": {"id": n0["id"]}}, "tuple": {"$": "tuple", "v": [n0["id"]]},
               "np_int64": {"$": "np", "v": ["int64", num]}, "np_float_nan": {"$": "np", "v": ["float64", "nan"]},
               "pd_na": {"$": "pd_na"}, "pd_nat": {"$": "pd_nat"}, "decimal": {"$": "decimal", "v": str(num)},
               "object": {"$": "object"}}
        for how, v in bad.items():
            self.add(s, "e3_single_nonstring", f"single|{how}", json.dumps(v), id_codec="enc", unit=how, rule="set",
                     expect=[], expect_error=E_TYPE, note=how)
        # a NumPy string is a str: the answer of the plain string
        ref = self.add(s, "e3_single_npstr_ref", f"single|str|{n0['id']}", n0["id"], rule="set", expect=n0["expect"])
        self.add(s, "e3_single_npstr", f"single|np_str|{n0['id']}", json.dumps({"$": "np", "v": ["str_", n0["id"]]}),
                 id_codec="enc", rule="set", expect=n0["expect"], note="np.str_")
        _ = ref
        items = [names[0], names[1], ids[0], names[2], ids[1]]
        batch = [q["id"] for q in items]
        ok = ["list", "tuple", "np_object", "np_str", "pd_series", "pd_series_string", "pd_index", "pd_categorical"]
        for cont in ok:
            for pos, q in enumerate(items):
                # 71690ab takes them all (its docs); d4c2dda's batch call takes lists and tuples only (TypeError)
                if cont in ("list", "tuple"):
                    j = dict(rule="set", expect=q["expect"])
                else:
                    j = dict(rule="set", expect=q["expect"], e_rule="or_error", or_error="d4c2dda:TypeError")
                self.add(s, "e3_batch_container", f"batch|{cont}|{pos}", q["id"], unit=cont, call="convert_many",
                         batch=batch, batch_container=cont, pos=pos, **j)
        for cont in ("generator", "set"):
            # an iterable that is not a sequence: 71690ab's docs name lists, tuples, arrays, Index and Series
            self.add(s, "e3_batch_iterable", f"batch|{cont}", batch[0], unit=cont, call="convert_many", batch=batch,
                     batch_container=cont, pos=None, rule="set", expect=[],
                     open="an iterable the docs do not name; an answer or an argument-type error are both documented",
                     e_rule="or_error", or_error=E_TYPE)
        self.add(s, "e3_batch_str", "batch|str", batch[0], unit="str", call="convert_many", batch=[batch[0]],
                 batch_container="str", pos=None, rule="set", expect=[], expect_error=E_TYPE)
        for how, v in (("none", None), ("nan", {"$": "float", "v": "nan"}), ("int", num), ("pd_na", {"$": "pd_na"}),
                       ("bytes", {"$": "bytes_hex", "v": b"TP53".hex()})):
            for pos in (0, 2):
                b = list(batch[:3])
                b.insert(pos, v)
                self.add(s, "e3_batch_nonstring_item", f"batchitem|{how}|{pos}", json.dumps(v), unit=how, call="convert_many",
                         batch=b, pos=pos, rule="set", expect=[], expect_error=E_TYPE, note=how)
        self.add(s, "e3_batch_empty", "batch|empty", "", unit="empty", call="convert_many", batch=[], pos=None, rule="set",
                 expect=[], open="an empty batch: an empty list of results is documented")
        # a pandas Series of strings holding a missing value (pd.NA, the 'string' dtype): never a gene
        b = [names[0]["id"], {"$": "pd_na"}, ids[0]["id"]]
        self.add(s, "e3_batch_pd_string_na", "batch|pdstring_na", json.dumps(b[1]), unit="pd_string_na", call="convert_many",
                 batch=b, batch_container="pd_series_string", pos=1, rule="set", expect=[], expect_error=E_TYPE)

    # ---------------- E4 another organism ----------------

    def e4(self, others):
        s = "e4"
        cur = self.cur
        r = random.Random(f"{SEED}:{self.org}:e4")
        for o in others:
            oc = T.Current(o)
            ids = sorted(oc.versioned.values())
            for gv in ids:
                self.add(s, "e4_gene_versioned", f"{o}|gv|{gv}", gv, unit=o, rule="set", expect=[], strategy="best",
                         open="an identifier of this organism in Ensembl's tables" if self.names.known(gv) else None)
            for g in r.sample(sorted(oc.genes), min(5000, len(oc.genes))):
                self.add(s, "e4_gene_unversioned", f"{o}|gu|{g}", g, unit=o, rule="set", expect=[], strategy="best",
                         open="an identifier of this organism in Ensembl's tables" if self.names.known(g) else None)
            for t in r.sample(sorted(oc.tx), min(3000, len(oc.tx))):
                v, _ = oc.tx[t]
                self.add(s, "e4_transcript", f"{o}|tx|{t}", f"{t}.{v}" if v else t, unit=o, rule="set", expect=[], strategy="best")
            for t in r.sample(sorted(oc.tl), min(1000, len(oc.tl))):
                v, _ = oc.tl[t]
                self.add(s, "e4_translation", f"{o}|tl|{t}", f"{t}.{v}" if v else t, unit=o, rule="set", expect=[], strategy="best")
            # symbols: the name of this organism's gene where the truth has exactly one (exact or by case alone)
            on = T.Names(o)
            syms = sorted(k for k, v in on.display115.items() if len(v) == 1)
            for sym in r.sample(syms, min(3000, len(syms))):
                mine = self.names.display115.get(sym, set())
                exact = self.names.exact115.get(sym, set())
                owners = self.names.owners_folded(sym)
                if len(mine) == 1 and {C.base_id(x) for x in exact} <= {C.base_id(x) for x in mine} \
                        and owners <= {C.base_id(x) for x in mine}:
                    self.add(s, "e4_symbol_exact", f"{o}|sym|{sym}", sym, unit=o, rule="set", expect=sorted(mine))
                elif not exact and len(owners) == 1:
                    (g,) = owners
                    gv = cur.versioned.get(g)
                    disp = {d for d, v in self.names.display115.items() if gv in v}
                    if gv and any(T.fold(d) == T.fold(sym) for d in disp):
                        self.add(s, "e4_symbol_case", f"{o}|sym|{sym}", sym, unit=o, rule="set", expect=[gv])
                    else:
                        self.add(s, "e4_symbol_other", f"{o}|sym|{sym}", sym, unit=o, rule="set", expect=[],
                                 open="another name source of one gene of this organism carries it")
                elif not owners:
                    self.add(s, "e4_symbol_absent", f"{o}|sym|{sym}", sym, unit=o, rule="set", expect=[],
                             open="not an identifier of this organism in any table the truth reads")
                else:
                    self.add(s, "e4_symbol_shared", f"{o}|sym|{sym}", sym, unit=o, rule="set", expect=[],
                             open="several genes of this organism carry it in some form")
            # NCBI gene ids of the other organism (IN-9): NCBI ids are unique over taxa, so none is an NCBI id here;
            # a bare number is read in every numeric database (rule `bare_numbers`): open both ways, the share answered reported
            ent = H.biomart_pairs(o, SNAP, "entrez", 1)
            if ent:
                nums = sorted(x for x in ent[1] if x.isdigit())
                for n in r.sample(nums, min(1500, len(nums))):
                    self.add(s, "e4_ncbi_foreign", f"{o}|ncbi|{n}", n, unit=o, rule="set", expect=[], strategy="best",
                             open="a number of another organism's NCBI gene: read in every numeric database (rule `bare_numbers`)"
                             if not self.names.known(n) else "an identifier of this organism in Ensembl's tables")

    # ---------------- E5 arguments ----------------

    def e5(self):
        s = "e5"
        first, last = WINDOW[self.org]
        names = self.harness("name_unique_to_gene", 1)
        ids = self.harness("gene_versioned_to_self", 1)
        if not names or not ids:
            return
        n0, i0 = names[0], ids[0]
        ref = {}

        def canon(key, ident, expect, **kw):
            ref[key] = self.add(s, "e5_reference", f"ref|{key}", ident, rule="set", expect=expect, strategy="best", **kw)["qid"]

        canon("name", n0["id"], n0["expect"])
        canon("id", i0["id"], i0["expect"])
        canon("id_db", i0["id"], [], final_database="UniProtKB Gene Name")
        self.sets[s][-1]["open"] = "the gene's UniProt gene name (the reference of its spellings)"
        canon("id_first", i0["id"], [], to_release=first)
        self.sets[s][-1]["open"] = "the reference of the release spellings"

        def rel(key, ident, against, **kw):
            return self.add(s, kw.pop("stratum", "e5_equal"), key, ident, rule="relation", expect=[], strategy="best",
                            relation={"op": "equal", "with": ref[against], "with_view": "best"}, **kw)

        def err(key, ident, ee, stratum="e5_error", **kw):
            return self.add(s, stratum, key, ident, rule="set", expect=[], expect_error=ee, strategy="best", **kw)
        x = i0["id"]
        # releases (to_release, from_release)
        for arg in ("to_release", "from_release"):
            for how, v, ee in (("below", first - 1, E_RELEASE), ("above", last + 1, E_RELEASE), ("zero", 0, E_RELEASE),
                               ("negative", -1, E_RELEASE), ("huge", 10 ** 9, E_RELEASE), ("str", str(last), E_TYPE),
                               ("float", float(last), E_TYPE), ("true", True, E_TYPE), ("false", False, E_TYPE),
                               ("nan", {"$": "float", "v": "nan"}, E_TYPE), ("list", [last], E_TYPE)):
                err(f"{arg}|{how}", x, ee, stratum=f"e5_release_{how}", kwargs={arg: v})
            for how, v in (("np_int64", {"$": "np", "v": ["int64", last]}), ("np_int32", {"$": "np", "v": ["int32", last]})):
                rel(f"{arg}|{how}", x, "id", kwargs={arg: v}, stratum="e5_release_numpy")
        err("to_release|np_float", x, E_TYPE, stratum="e5_release_float", kwargs={"to_release": {"$": "np", "v": ["float64", last]}})
        err("to_release|np_bool", x, E_TYPE, stratum="e5_release_true", kwargs={"to_release": {"$": "np", "v": ["bool_", True]}})
        rel("to_release|first", x, "id_first", kwargs={"to_release": {"$": "np", "v": ["int64", first]}}, stratum="e5_release_numpy")
        # backward and in place: valid (an answer, never an error)
        self.add(s, "e5_backward", "backward|last_first", x, rule="set", expect=[], strategy="best", from_release=last,
                 to_release=first, open="a valid backward conversion (judged by family A)")
        self.add(s, "e5_backward", "backward|mid", x, rule="set", expect=[], strategy="best", from_release=last,
                 to_release=(first + last) // 2, open="a valid backward conversion (judged by family A)")
        # final_database spellings
        for how, v in (("lower", "uniprotkb gene name"), ("upper", "UNIPROTKB GENE NAME"), ("space", "  UniProtKB Gene Name \t")):
            rel(f"final_database|{how}", x, "id_db", kwargs={"final_database": v}, stratum="e5_database_spelling")
        rel("final_database|ensembl_gene", x, "id", kwargs={"final_database": "ensembl_gene"}, stratum="e5_database_spelling")
        rel("final_database|ENSEMBL_GENE", x, "id", kwargs={"final_database": "ENSEMBL_GENE"}, stratum="e5_database_spelling")
        for how, v, ee in (("unknown", "NotADatabase", E_DB), ("empty", "", E_DB), ("blank", "   ", E_DB),
                           ("partial", "UniProt", E_DB), ("int", 5, E_TYPE), ("list", ["UniProtKB Gene Name"], E_TYPE)):
            err(f"final_database|{how}", x, ee, stratum=f"e5_database_{'type' if ee == E_TYPE else 'unknown'}", kwargs={"final_database": v})
        # rule `no_synonym_target`: the synonym prefix is no final database (71690ab offers it)
        for db in ("synonym_id::NCBI gene (formerly Entrezgene)",):
            self.add(s, "e5_database_synonym", f"final_database|synonym|{db}", x, rule="set", expect=[], strategy="best",
                     kwargs={"final_database": db}, open="71690ab offers the synonym databases as targets (rule `no_synonym_target` not built)",
                     dec=dict(open=None, expect_error=E_DB))
        # from_database
        rel("from_database|name_db", n0["id"], "name", kwargs={"from_database": NAME_DB[self.org].lower()}, stratum="e5_database_spelling")
        for how, v, ee in (("unknown", "NotADatabase", E_DB), ("empty", "", E_DB), ("int", 5, E_TYPE)):
            err(f"from_database|{how}", n0["id"], ee, stratum=f"e5_database_{'type' if ee == E_TYPE else 'unknown'}",
                kwargs={"from_database": v})
        # budgets (71690ab on): positive finite numbers; states an integer
        for arg in ("budget_states", "budget_seconds"):
            for how, v, ee in (("zero", 0, E_INVALID), ("negative", -1, E_INVALID), ("nan", {"$": "float", "v": "nan"}, None),
                               ("inf", {"$": "float", "v": "inf"}, None), ("neginf", {"$": "float", "v": "-inf"}, None),
                               ("true", True, E_TYPE), ("str", "100", E_TYPE), ("list", [1], E_TYPE)):
                if ee is None:
                    ee = E_TYPE if arg == "budget_states" else E_INVALID
                err(f"{arg}|{how}", x, ee.split("|")[0], stratum=f"e5_budget_{'type' if ee.startswith('Argument') else 'invalid'}",
                    kwargs={arg: v}, generations=["convert"])
        err("budget_states|float", x, "ArgumentTypeError", stratum="e5_budget_type", kwargs={"budget_states": 1.5}, generations=["convert"])
        rel("budget_states|huge", x, "id", kwargs={"budget_states": 10 ** 12}, generations=["convert"], stratum="e5_budget_valid")
        rel("budget_seconds|huge", x, "id", kwargs={"budget_seconds": 1e9}, generations=["convert"], stratum="e5_budget_valid")
        rel("budget_states|np", x, "id", kwargs={"budget_states": {"$": "np", "v": ["int64", 100000]}}, generations=["convert"],
            stratum="e5_budget_valid")
        self.add(s, "e5_budget_tiny", "budget_seconds|tiny", x, rule="set", expect=i0["expect"], strategy="best",
                 kwargs={"budget_seconds": 1e-9}, generations=["convert"], note="an answer or the budget error (a timeout)")
        self.add(s, "e5_budget_tiny", "budget_states|one", n0["id"], rule="set", expect=n0["expect"], strategy="best",
                 kwargs={"budget_states": 1}, generations=["convert"], note="an answer or the budget error (a timeout)")
        # switches
        for arg in ("explain", "allow_bridges"):
            for how, v in (("int1", 1), ("str", "yes"), ("none", None)):
                err(f"{arg}|{how}", x, "ArgumentTypeError", stratum="e5_switch_type", kwargs={arg: v}, generations=["convert"])
        rel("explain|np_bool", x, "id", kwargs={"explain": {"$": "np", "v": ["bool_", True]}}, generations=["convert"],
            stratum="e5_switch_valid")
        # infer_source: modes and switches
        batch = [n0["id"], x]
        for mode in ("complete", "release", "assembly", "assembly_release"):
            self.add(s, "e5_infer_mode", f"infer|mode|{mode}", mode, call="infer_source", batch=batch, kwargs={"mode": mode},
                     rule="set", expect=[], strategy="best", open="a valid mode (source inference is judged by family D)")
        for how, v in (("unknown", "x"), ("none", None), ("int", 1), ("old_name", "ensembl_release")):
            self.add(s, "e5_infer_mode_invalid", f"infer|mode|{how}", str(v), call="infer_source", batch=batch,
                     kwargs={"mode": v}, rule="set", expect=[], strategy="best", expect_error=E_INVALID,
                     generations=["convert"] if how == "old_name" else None)
        self.add(s, "e5_infer_switch", "infer|winner|str", "x", call="infer_source", batch=batch,
                 kwargs={"report_only_winner": "yes"}, rule="set", expect=[], strategy="best", expect_error="ArgumentTypeError",
                 generations=["convert"])
        for how, b in (("empty", []), ("junk", ["☃not_a_gene☃", "zz_nothing_zz"])):
            self.add(s, "e5_infer_nothing", f"infer|{how}", how, call="infer_source", batch=b, rule="set", expect=[],
                     strategy="best", expect_error=E_UNKNOWN_ID)
        self.add(s, "e5_infer_container", "infer|str", "x", call="infer_source", batch=[n0["id"]], batch_container="str",
                 rule="set", expect=[], strategy="best", expect_error=E_TYPE, generations=["convert"])
        # the other calls of the api (71690ab on)
        gen = ["convert"]
        api = [("names_of_int", "names_of", [5], {}, "ArgumentTypeError"),
               ("names_of_release", "names_of", [i0["expect"][0]], {"release": last + 1}, "UnknownReleaseError"),
               ("names_of_release_str", "names_of", [i0["expect"][0]], {"release": str(last)}, "ArgumentTypeError"),
               ("genes_named_none", "genes_named", [None], {}, "ArgumentTypeError"),
               ("genes_named_fold_str", "genes_named", [n0["id"]], {"fold_case": "yes"}, "ArgumentTypeError"),
               ("gene_to_db_unknown", "gene_to_database", [i0["expect"][0], "NotADatabase"], {}, "UnknownDatabaseError"),
               ("gene_to_db_release", "gene_to_database", [i0["expect"][0], "UniProtKB Gene Name"], {"release": first - 1},
                "UnknownReleaseError"),
               ("final_databases_release", "final_databases", [last + 1], {}, "UnknownReleaseError"),
               ("final_databases_str", "final_databases", [str(last)], {}, "ArgumentTypeError")]
        for how, m, args, kw, ee in api:
            self.add(s, "e5_api_error", f"api|{how}", how, call="api", method=m, args=args, kwargs=kw, rule="set", expect=[],
                     strategy="best", expect_error=ee, generations=gen)
        for how, m, args, kw in (("names_of", "names_of", [i0["expect"][0]], {}), ("genes_named", "genes_named", [n0["id"]], {}),
                                 ("describe", "describe", [], {}), ("releases", "releases", [], {})):
            self.add(s, "e5_api_valid", f"api|{how}", how, call="api", method=m, args=args, kwargs=kw, rule="set", expect=[],
                     strategy="best", open="a valid call (no error expected)", generations=gen)
        # classify: the strategy
        self.add(s, "e5_classify", "classify|strategy", "x", call="classify", batch=batch, kwargs={"strategy": "most"},
                 rule="set", expect=[], strategy="best", expect_error="InvalidArgumentError", generations=gen)
        self.add(s, "e5_classify", "classify|valid", "x", call="classify", batch=batch, kwargs={"strategy": "all"},
                 rule="set", expect=[], strategy="best", open="a valid call (no error expected)", generations=gen)

    # ---------------- E8 static relations ----------------

    def e8(self):
        s = "e8"
        # identity: every gene id at every release it exists at, to that release (the GTF of that release)
        self.identity = True
        # versioned and unversioned forms of a current id give one gene
        for b in self.harness("gene_versioned_to_self", 3000):
            a = self.add(s, "e8_versioned_pair", f"vpair|v|{b['id']}", b["id"], rule="set", expect=b["expect"], strategy="best")
            self.add(s, "e8_versioned_pair", f"vpair|u|{b['id']}", C.base_id(b["id"]), rule="relation", expect=[], strategy="best",
                     relation={"op": "equal", "with": a["qid"], "with_view": "best"})
        # a unique name and its case variants (where no other gene carries the variant in any form)
        for b in self.harness("name_unique_to_gene", 3000):
            x = b["id"]
            gl = {C.base_id(e) for e in b["expect"]}
            a = None
            for how, f in (("lower", x.lower()), ("upper", x.upper())):
                if f == x or not (self.names.owners_folded(f) <= gl):
                    continue
                if a is None:
                    a = self.add(s, "e8_case_pair", f"case|ref|{x}", x, rule="set", expect=b["expect"], strategy="best")
                self.add(s, "e8_case_pair", f"case|{how}|{x}", f, rule="relation", expect=[], strategy="best", unit=x,
                         relation={"op": "equal", "with": a["qid"], "with_view": "best"},
                         open=None)
        # white space trimmed equals untrimmed (every kind)
        for kind in ("name_unique_to_gene", "gene_versioned_to_self", "ncbi_to_gene", "transcript_to_gene", "synonym_to_gene"):
            for b in self.harness(kind, 400):
                a = self.add(s, "e8_whitespace_pair", f"ws|ref|{kind}|{b['id']}", b["id"], rule="set", expect=b["expect"], strategy="best")
                self.add(s, "e8_whitespace_pair", f"ws|pad|{kind}|{b['id']}", f" \t{b['id']} \n", rule="relation", expect=[],
                         strategy="best", unit=C.base_id(b["id"]), relation={"op": "equal", "with": a["qid"], "with_view": "best"})
        # infer_source unchanged by order and by duplicates
        r = random.Random(f"{SEED}:{self.org}:e8infer")
        for rel in sorted(r.sample(H.gtf_releases(self.org), min(6, len(H.gtf_releases(self.org))))):
            g = H.gtf_genes(self.org, rel, "primary") or {}
            for what in ("ids", "names"):
                pool = sorted(g) if what == "ids" else sorted({v["gene_name"] for v in g.values() if v.get("gene_name")})
                for n in (100, 2000):
                    if len(pool) < n:
                        continue
                    lst = r.sample(pool, n)
                    shuffled = lst[:]
                    r.shuffle(shuffled)
                    doubled = lst + lst[: n // 2]
                    key = f"{what}|{rel}|{n}"
                    for winner in (True, False):
                        a = self.add(s, "e8_infer_order", f"infer|ref|{key}|{winner}", key, call="infer_source", batch=lst,
                                     kwargs={"mode": "release", "report_only_winner": winner}, rule="set", expect=[], strategy="best",
                                     open="the reference of the relation (source inference is judged by family D)")
                        self.add(s, "e8_infer_order", f"infer|shuffled|{key}|{winner}", key, call="infer_source", batch=shuffled,
                                 kwargs={"mode": "release", "report_only_winner": winner}, rule="relation", expect=[], strategy="best",
                                 relation={"op": "equal", "with": a["qid"], "with_view": "best"})
                        if winner:
                            self.add(s, "e8_infer_duplicates", f"infer|doubled|{key}", key, call="infer_source", batch=doubled,
                                     kwargs={"mode": "release", "report_only_winner": True}, rule="relation", expect=[],
                                     strategy="best", relation={"op": "equal", "with": a["qid"], "with_view": "best"})

    def identity_lines(self):
        """Every gene id at every release it exists at (GTF), versioned where the GTF has versions, to that release."""
        for rel in H.gtf_releases(self.org):
            if not (WINDOW[self.org][0] <= rel <= WINDOW[self.org][1]):
                continue
            for g, r in sorted((H.gtf_genes(self.org, rel) or {}).items()):
                v = r.get("gene_version") or ""
                ident = f"{g}.{v}" if v else g
                yield dict(qid=f"e|e8_identity|{rel}|{ident}", family=FAMILY, kind="e8_identity",
                           stratum=f"e8_identity_{'versioned' if v else 'unversioned'}_{r['region']}", unit=g, id=ident,
                           from_release=rel, to_release=rel, final_database=None, strategy="best", rule="set",
                           cmp="versioned" if v else "base", expect=[ident], broad=None, locus=None)

    # ---------------- the input hunt's reproducers ----------------

    def hunt(self):
        s = "hunt_inputs"
        o = self.org
        rep = []
        if o == "homo_sapiens":
            rep = [("IN-1", "nan", "single"), ("IN-1", "NA", "single"), ("IN-1", "NaT", "single"), ("IN-1", "unknown", "single"),
                   ("IN-3", "SRF.1", "mu"), ("IN-3", "MB.1", "mu"), ("IN-3", "C4B.1", "mu"), ("IN-3", "AIP_1", "mu"),
                   ("IN-4", "﻿TP53", "zw"), ("IN-4", "TP53​", "zw"), ("IN-4", "HLA–A", "dash"), ("IN-4", '"TP53"', "punct"),
                   ("IN-4", "TP５３", "full"), ("IN-5", "enst00000016171", "tx"), ("IN-6", "synonym_id::TP53", "internal"),
                   ("IN-6", "ENSG00000141510.Void", "void"), ("IN-6", "ENSG00000141510.void", "void"),
                   ("IN-8", "ENSG00000141510." + "9" * 5000, "longv"), ("IN-9", "56386", "num"), ("IN-11", "ADORA3;", "punct")]
        elif o == "mus_musculus":
            rep = [("IN-1", "-", "single"), ("IN-1", "nan", "single"), ("IN-1", "unknown", "single"), ("IN-2", "AR", "case"),
                   ("IN-2", "CCK", "case"), ("IN-2", "AGT", "case"), ("IN-2", "GRIN1", "case"), ("IN-2", "nos2", "case"),
                   ("IN-10", "AY172335_1", "sep"), ("IN-11", "'Agt'", "punct")]
        elif o == "sus_scrofa":
            rep = [("IN-1", "N/A", "single"), ("IN-2", "afp", "case"), ("IN-2", "Has2", "case")]
        truth = {}
        for x in ("TP53", "HLA-A", "SRF", "MB", "C4B", "AIP", "ADORA3", "Agt", "Ar", "Cck", "Grin1", "Nos2", "AFP", "HAS2"):
            g = self.names.display115.get(x)
            if g and len(g) == 1:
                truth[x] = sorted(g)
        tx = self.cur.tx.get("ENST00000016171")
        for fid, x, how in rep:
            key = f"{fid}|{x[:60]}|{len(x)}"
            base_name = {"SRF.1": "SRF", "MB.1": "MB", "C4B.1": "C4B", "AIP_1": "AIP", "﻿TP53": "TP53", "TP53​": "TP53",
                         "HLA–A": "HLA-A", '"TP53"': "TP53", "TP５３": "TP53", "ADORA3;": "ADORA3", "'Agt'": "Agt",
                         "AR": "Ar", "CCK": "Cck", "AGT": "Agt", "GRIN1": "Grin1", "nos2": "Nos2", "afp": "AFP",
                         "Has2": "HAS2", "synonym_id::TP53": "TP53"}.get(x)
            E = truth.get(base_name, []) if base_name else []
            if how == "single":
                j52 = x.strip().casefold() in ("nan", "na", "n/a", "none", "null", "nat", "unknown")
                self.add(s, f"hunt_{fid}", key, x, rule="set", expect=[], unit=fid,
                         open="71690ab reads it as a name (rule `missing_values` not built)" if j52 else
                         "a placeholder string (the decided builder rule (rule `missing_values`) drops strings without a letter or digit)",
                         dec=dict(open=None, expect_error=E_MISSING) if j52 else None)
            elif how == "case" and E and o == "mus_musculus":
                self.add(s, f"hunt_{fid}", key, x, rule="set", expect=E, unit=fid,
                         open="71690ab reads an exact synonym before a case match (rule `single_case_symbols` not built)",
                         dec=dict(open=None))
            elif how in ("mu", "punct", "case") and E:
                dec = dict(rule="set", expect=E) if how == "punct" else None
                self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=E, unit=fid, dec=dec)
            elif how in ("zw", "dash", "full") and E:
                self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=E, unit=fid, dec=dict(rule="set", expect=E))
            elif how == "tx" and tx:
                self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=[tx[1]], unit=fid, dec=dict(rule="set", expect=[tx[1]]))
            elif how == "internal" and E:
                self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=E, unit=fid,
                         dec=dict(rule="set", expect=[], e_rule="or_error", or_error=E_INPUT))
            elif how == "void":
                gv = self.cur.versioned.get("ENSG00000141510")
                self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=[gv] if gv else [], unit=fid,
                         dec=dict(rule="set", expect=[], e_rule="or_error", or_error=E_INPUT))
            elif how == "longv":
                gv = self.cur.versioned.get("ENSG00000141510")
                self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=[gv] if gv else [], unit=fid,
                         dec=dict(rule="set", expect=[], expect_error=E_INPUT))
            elif how == "num":
                self.add(s, f"hunt_{fid}", key, x, rule="set", expect=[], unit=fid,
                         open="a bare number: read in every numeric database (rule `bare_numbers`)")
            elif how == "sep":
                g1 = self.names.exact115.get("AY172335-1", set())
                g2 = self.names.exact115.get("AY172335.1", set())
                if g1 and g2:
                    self.add(s, f"hunt_{fid}", key, x, rule="safe", expect=sorted(g1 | g2), unit=fid,
                             dec=dict(rule="set_exact", expect=sorted(g1 | g2)))
        # IN-7: name lists of a known release (rule `inference_per_release`: the band of releases holds the release)
        for rel in {"homo_sapiens": (55, 75, 80, 105, 110), "mus_musculus": (60, 75, 110)}.get(o, ()):
            g = H.gtf_genes(o, rel, "primary") or {}
            nm = sorted({v["gene_name"] for v in g.values() if v.get("gene_name")})
            if nm:
                self.add(s, "hunt_IN-7", f"IN-7|names|{rel}", f"names of GTF {rel}", call="infer_source", batch=nm,
                         kwargs={"mode": "release"}, rule="set", expect=[str(rel)], strategy="best", unit=f"IN-7|{rel}",
                         e_rule="contains", generations=["convert"])

    # ---------------- publishing ----------------

    def publish(self, which, dry):
        sources = [os.path.join(C.TRUTH, self.org, "SOURCES.tsv"), "harness published sets (accuracy/queries)"]
        out = {}
        for s in which:
            ln = self.sets.get(s, [])
            if not ln:
                continue
            desc = {"e1": "E1 missing values", "e2": "E2 malformed strings", "e2_case": "E2 mouse symbols in another case (IN-2)",
                    "e2_all": "E2 every display name with stray punctuation (IN-11) and make-unique suffixes (IN-3)",
                    "e3": "E3 types", "e4": "E4 another organism's identifiers", "e5": "E5 arguments at and beyond their bounds",
                    "e8": "E8 static relations", "hunt_inputs": "the input hunt's reproducers"}[s]
            m = L.publish(FAMILY, self.org, s, ln, desc + " (judge: 71690ab's docs)", __file__, judge="now", sources=sources, dry=dry)
            out[s] = m["n"]
            rev = []
            by = {q["qid"]: q for q in ln}
            for qid, ch in self.dec.get(s, {}).items():
                q = dict(by[qid])
                for k, v in ch.items():
                    if v is None:
                        q.pop(k, None)
                    else:
                        q[k] = v
                q["revised_by"] = ["decided"]
                rev.append(q)
            if rev:
                m = L.publish(FAMILY, self.org, s + "_dec", rev, desc + " (judge: the decisions not built yet: 52, 56, 79-87)",
                              __file__, revises=[s], judge="decided", sources=sources, dry=dry)
                out[s + "_dec"] = m["n"]
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="e1,e2,e2_case,e3,e4,e5,e8,hunt_inputs,e8_identity")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    which = a.sets.split(",")
    g = Gen(a.organism)
    for s in which:
        if s == "e4":
            g.e4([o for o in C.ORGANISMS if o != a.organism])
        elif s == "hunt_inputs":
            g.hunt()
        elif s == "e8_identity":
            continue
        else:
            getattr(g, s)()
    print(json.dumps(g.publish([s for s in which if s != "e8_identity"], a.dry)))
    if "e8_identity" in which:
        m = L.publish(FAMILY, a.organism, "e8_identity", g.identity_lines(), "E8 identity: every gene id at every release "
                      "it exists at (GTF), to that release", __file__, judge="now",
                      sources=[os.path.join(C.TRUTH, a.organism, "SOURCES.tsv")], dry=a.dry, fmt="parquet")
        print(json.dumps({"e8_identity": m["n"], "counts": m["counts"]}))


if __name__ == "__main__":
    main()
