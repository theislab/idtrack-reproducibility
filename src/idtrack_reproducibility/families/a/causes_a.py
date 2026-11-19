"""Family A's causes (a plugin of F0's causes.py: `--plugin families/a/causes_a.py`).

- `core_version_on_versionless_release` (class d, the truth): A*'s `TruthA` adds the core genes the GTF lacks (human
  48 to 51) with the core table's version, while the GTF of those releases writes no versions; the comparison of such
  a release (`cmp: base`) strips the answer's version but not the expected one, so the expected `<id>.<v>` never
  equals the answer (the same gene); and an answer within the broad walk by stable id is wrong instead of
  undetermined for the same reason.
"""
import os
import sys

from idtrack_reproducibility.families import causes as K


def core_version_on_versionless_release(c):
    if not c.got:
        return None
    versioned = [x for x in c.exp + c.broad if "." in x]
    if not versioned:
        return None
    if c.gotb <= c.expb and ("missing" in (c.r.get("reason") or "") or "not expected" in (c.r.get("reason") or "")):
        return f"{c.id} -> {','.join(sorted(c.got))[:100]}: expected with a version the release's GTF does not write"
    broadb = {K.base(x) for x in c.broad} | c.expb
    if c.broad and c.gotb <= broadb:
        # within the broad walk by stable id (an undetermined answer), wrong only because the broad set carries the
        # core table's version of a core-only gene
        return f"{c.id} -> {','.join(sorted(c.got))[:100]}: within the broad walk; its core-only genes carry a version"
    return None


def lrg_version_not_normalised(c):
    """The judge compares by stable id at a release whose GTF writes no versions, but its normalisation of an id
    (`common.base_id`) knows `LRG_<n>.<v>` only, not LRG gene ids `LRG_<n>_g1.<v>`: the answer `LRG_10_g1.1` is the
    expected `LRG_10_g1`."""
    if c.got and all(x.startswith("LRG_") for x in c.got) and \
            {x.rsplit(".", 1)[0] for x in c.got} == {x.rsplit(".", 1)[0] if x.rsplit(".", 1)[-1].isdigit() else x for x in c.exp}:
        return f"{c.id} -> {','.join(c.got)}: LRG gene id with a version, expected without"
    return None


def transcript_read_elsewhere(c):
    """A transcript or translation id asked from a release R is answered by reading it at another release or as a
    synonym of another database (the explained path's match and release), giving another gene than the GTF's gene of
    the id at R."""
    if not str(c.kind).startswith("a7_"):
        return None
    m = str(c.res.get("match") or "")
    if c.read_at != c.fr or m.startswith("synonym"):
        return f"{c.id} from {c.fr}: read at {c.read_at} (match {m}); got {','.join(sorted(c.gotb)) or 'nothing'}"
    return None


# causes.py runs as `__main__`: register there (its TESTS are the ones applied) and in the importable module
import __main__  # noqa: E402

for mod in {id(m): m for m in (K, __main__) if hasattr(m, "register") and hasattr(m, "TESTS")}.values():
    mod.register("core_version_on_versionless_release", "d", "A*'s TruthA (core-only genes, human 48-51)",
                 "expected id carries the core table's version at a release whose GTF writes none: the same gene",
                 core_version_on_versionless_release, before="external_bridge")
    mod.register("lrg_version_not_normalised", "d", "judge: common.base_id",
                 "an LRG gene id `LRG_<n>_g1.<v>` is not reduced to its stable id at a release compared by stable id",
                 lrg_version_not_normalised, before="external_bridge")
    mod.register("transcript_read_elsewhere", "a", "family A's reading of its wrong answers",
                 "a transcript/translation id asked from R is read at another release or as a synonym: another gene",
                 transcript_read_elsewhere, before="external_bridge")
