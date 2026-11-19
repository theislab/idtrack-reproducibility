"""Family B's plugin of F0's classifier (`../causes.py --plugin b/causes_b.py`): B's kinds, and B's own causes.

The residual reading's tests (F0) read name lines through BioMart at the release asked about; B's lines are read at
every release, so B adds tests that read its own truth (`truth_b.NamesB`: GTF names, the core xrefs of the release,
BioMart, the reach) where BioMart has no table. Each test returns a detail or None; the first that holds names the
cause (F0's order; B's are placed before F0's generic `name_in_other_source`).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import causes as K

K.register_kinds(names=("b1_name_at_release", "b1_name_to_current", "b1_name_default", "b4_name_to_current",
                        "b_string_default", "b5_synonym", "b5_synonym_to_current", "b6_form", "b7_form"),
                 databases=("b8_name_db",))
