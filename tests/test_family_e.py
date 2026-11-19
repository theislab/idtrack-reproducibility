"""Family E (robustness): the encoding of non-string inputs and the folding of names, on invented values."""
import math

from idtrack_reproducibility.families.e import adapter_e as AE
from idtrack_reproducibility.families.e import gen_e as GE
from idtrack_reproducibility.families.e import truth_e as TE


def test_non_string_values_round_trip_through_their_encoding():
    for v in (float("nan"), float("inf"), float("-inf")):
        back = AE.decode(GE.enc(v))
        assert (math.isnan(back) and math.isnan(v)) or back == v
    assert GE.enc("TP53") == "TP53" and AE.decode("TP53") == "TP53"
    assert AE.decode({"$": "bytes_hex", "v": "414243"}) == b"ABC"


def test_names_are_folded_as_the_case_and_separator_steps_fold():
    assert TE.fold("Hla-A") == TE.fold("HLA_A") == TE.fold("hla.a") == "hla-a"
    assert TE.fold("TP53") != TE.fold("TP5-3")


def test_expected_error_classes_name_both_api_generations():
    for e in (GE.E_RELEASE, GE.E_DB, GE.E_INVALID, GE.E_UNKNOWN_ID):
        assert "|" in e and e.split("|")[1].startswith("d4c2dda:")
