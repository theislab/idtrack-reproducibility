"""Survey of rosetta's per-feature tables of its 33 authors' matrices (read only):
rosetta_storage/scratch/data_discovery_to_unification/r3_build_matrix_locus_dev/out/m02/<dataset>.tsv.gz.
Per dataset and source: species (by the ids it reached), features, the kinds of features (route), source release,
same assembly, counts. Usage: survey_m02.py <out tsv>"""
import sys
from pathlib import Path

import pandas as pd

from idtrack_reproducibility import config as _CFG

M02 = Path(_CFG.get().source("rosetta_m02"))
rows = []
for f in sorted(M02.glob("*.tsv.gz")):
    t = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    ds = f.name[:-len(".tsv.gz")]
    for src, part in t.groupby("source"):
        read = part["route"].str.split("+").str[0]
        ens = part["author_id"].str.match(r"^ENS[A-Z]*G\d{11}")
        rows.append({"dataset": ds, "source": src, "features": len(part), "ensembl_ids": int(ens.sum()), "other": int((~ens).sum()),
                     "source_release": ",".join(sorted(set(part["source_release"]))), "same_assembly": ",".join(sorted(set(part["same_assembly"]))),
                     "counts": int(pd.to_numeric(part["counts"]).sum()), "routes": ";".join(f"{k}:{v}" for k, v in read.value_counts().items()),
                     "classes": ";".join(f"{k}:{v}" for k, v in part["class"].value_counts().items()),
                     "species_hint": ",".join(sorted(set(part["target_id"].str.extract(r"^(ENS[A-Z]*G)")[0].dropna())))})
out = pd.DataFrame(rows)
out.to_csv(sys.argv[1], sep="\t", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 60, "display.max_rows", 100)
print(out[["dataset", "features", "ensembl_ids", "other", "source_release", "same_assembly", "counts", "species_hint"]].to_string())
print("total features", out["features"].sum(), "ensembl", out["ensembl_ids"].sum(), "other", out["other"].sum())
