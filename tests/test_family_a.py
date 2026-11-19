"""Family A (history travel, exhaustive): its report's axes and the exceptions filter, on invented inputs."""
from idtrack_reproducibility.families.a import gen_a as GA
from idtrack_reproducibility.families.a import report_a as RA


def test_distance_buckets_are_symmetric_and_cover_every_distance():
    assert RA.bucket(0) == "0" and RA.bucket(-1) == "1" and RA.bucket(1) == "1"
    assert RA.bucket(5) == "2-5" and RA.bucket(-40) == "21-40" and RA.bucket(67) == "41+"


class FakeExceptions:
    def __init__(self, rows):
        self._rows, self.cache = rows, {}

    def rows(self, rel):
        return self._rows

    def project(self, place, rel):
        return dict(place, projected=tuple(sorted(self.cache[rel])))


def test_rows_without_coordinates_are_left_out_and_remembered():
    good = dict(start="1", end="10", ref_start="5", ref_end="14")
    bad = dict(start="na", end="na", ref_start="5", ref_end="14")
    ex = GA.ExceptionsA(FakeExceptions({"HAP1": [good, bad], "HAP2": [bad]}))
    assert ex.rows(3) == {"HAP1": [good]}
    assert ex.na[3] == {"HAP1", "HAP2"} and ex.has(3)
    assert ex.project(dict(seqname="HAP1"), 3)["projected"] == ("HAP1",)  # projected with the filtered rows only
