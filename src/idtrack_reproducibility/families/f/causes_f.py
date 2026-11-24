"""Family F's plugin of F0's classifier (`../causes.py --plugin f/causes_f.py`): F's kinds, and F's own causes.

F's lines are A's and B's lines asked through the harmonizer (one batch per dataset, the best view): the residual
reading's tests (F0) read them as A's and B's lines are read, the id forms as history and the symbols as names (a mix
dataset's line by its id). The explained path is the single call's (`explain.py`); F adds the causes that only a
batch of a dataset has:

- `copy_of_another_feature` (class b, rule `made_unique_copies`): the feature is a made-unique form `X.1`/`X-1` of
  another feature X of its dataset, whose answer it takes in a batch (the single call reads the form on its own and
  gives another answer, `same_answer` false): a wrong copy has X's answer, and X's own line names the cause;
- `batch_differs_from_single` (class a, rule `batch_agrees_with_single`): any other feature whose explained single
  call gives another answer than the harmonizer's batch gave;
- `lrg_record_as_gene` (class a, rule `lrg_copies`): the answer names an LRG record (Ensembl's `LRG_<n>` scheme) as
  the feature's gene; LRG records follow rule `copies_apart` (reported apart, `alternative_locus_genes`), so the
  harmonizer makes an output feature named by the record (human `IKBKG` and `ORAI1` read at 65: `LRG_70.1`,
  `LRG_93.1`; the single call `API.convert` answers so too);
- `best_keeps_part_of_shared_name` (class b, rule `ties`): the harmonizer reads the best view; a name several genes
  carry, judged `set_exact` (the truth wants every carrier), is answered with a non-empty part of them, or of them and
  the wider set Ensembl's history also allows (`broad`), but not with all of them (the top-ranked; the others are
  candidates of the `all` view, which the harmonizer does not read). Placed after F0's
  `support_ranks_shared_name`, which explains the part where the support of the name sources ranks it.
"""
import json

from idtrack_reproducibility.families import causes as K
from idtrack_reproducibility.families.causes import base



def _same_answer(c):
    return c.same_answer


def t_copy_of_another_feature(c):
    line = c.r
    if line.get("stratum") != "f1_made_unique" or _same_answer(c) is not False:
        return None
    return f"a copy of a feature of the dataset: the batch answers the copy as the feature it copies ({line.get('reason', '')[:80]})"


def t_batch_differs_from_single(c):
    if _same_answer(c) is not False:
        return None
    q = json.loads(c.r["query"]) if c.r.get("query") else {}
    return f"the single call of {q.get('id')!r} answers otherwise than the dataset's batch"


def t_best_keeps_part_of_shared_name(c):
    if c.r.get("variant", "").split(",")[0] != "best" or not c.gotb or not c.missing or len(c.expb) < 2:
        return None
    wider = {base(x) for x in c.broad}
    if c.extra - wider:
        return None
    within = f", {len(c.extra)} of them in the wider set Ensembl's history also allows" if c.extra else ""
    return f"best keeps {len(c.gotb)} genes of the {len(c.expb)} carrying the name{within}"


LRG_RECORD = __import__("re").compile(r"^LRG_\d+")  # Ensembl's identifier scheme of LRG records


def t_lrg_record_as_gene(c):
    lrg = sorted(x for x in c.got if LRG_RECORD.match(str(x)))
    return f"the answer names the LRG record(s) {','.join(lrg)[:80]} as the gene" if lrg else None


# The classifier runs as `python -m ...causes`: its tests are those of `__main__`, a module apart from the imported
# `causes`; F's kinds and causes go into both (as family D's plugin does).
import __main__  # noqa: E402

for mod in {id(m): m for m in (K, __main__) if hasattr(m, "register") and hasattr(m, "TESTS")}.values():
    mod.register_kinds(history=("f1_ids", "f1_versioned"), names=("f1_symbols",))
    mod.register("copy_of_another_feature", "b", "rule `made_unique_copies`",
                 "a made-unique copy takes the answer of the feature it copies", t_copy_of_another_feature,
                 before="rule_decided_after_tree")
    mod.register("batch_differs_from_single", "a", "rule `batch_agrees_with_single`",
                 "the harmonizer's batch answers a feature otherwise than its single call", t_batch_differs_from_single,
                 before="rule_decided_after_tree")
    mod.register("lrg_record_as_gene", "a", "rule `lrg_copies`",
                 "an LRG record answered as the feature's gene (LRG records are copies reported apart)", t_lrg_record_as_gene,
                 before="rule_decided_after_tree")
    mod.register("best_keeps_part_of_shared_name", "b", "rule `ties`",
                 "the harmonizer reads `best`: of a name several genes carry it keeps the top-ranked part",
                 t_best_keeps_part_of_shared_name, before="synonym_spelling_merged")
