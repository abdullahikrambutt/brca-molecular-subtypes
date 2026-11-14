"""Create small synthetic files with the same format as the TCGA Xena files.

Used only by `python run_all.py --test` to check that the pipeline runs
end to end. The numbers mean nothing biologically.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import config as C  # noqa: E402

assert C.TEST_MODE, "Run via `python run_all.py --test` so real data folders are untouched."
rng = np.random.default_rng(0)

n_tum, n_norm, n_genes = 500, 30, 3000
pam50_names = ["ESR1", "PGR", "FOXA1", "ERBB2", "GRB7", "KRT5", "KRT14", "KRT17",
               "FOXC1", "MKI67", "BIRC5", "CCNB1", "MYBL2", "BCL2", "MLPH", "NAT1"]
genes = pam50_names + [f"GENE{i}" for i in range(n_genes - len(pam50_names))]
classes = np.array(["LumA", "LumB", "Basal", "Her2", "Normal"])
y = rng.choice(classes, n_tum, p=[0.40, 0.20, 0.20, 0.12, 0.08])

base = np.clip(rng.normal(6, 2.5, n_genes), 0, None)
X = base + rng.normal(0, 0.8, (n_tum, n_genes))
for k, cls in enumerate(classes):
    sig = np.r_[rng.choice(len(pam50_names), 4, replace=False),
                rng.choice(np.arange(len(pam50_names), n_genes), 120, replace=False)]
    X[np.ix_(y == cls, sig)] += rng.normal(0, 2.0, len(sig))
X = np.clip(X, 0, None)
Xn = np.clip(base + rng.normal(0, 0.8, (n_norm, n_genes)), 0, None)

tum_ids = [f"TCGA-XX-{i:04d}-01" for i in range(n_tum)]
norm_ids = [f"TCGA-XX-{i:04d}-11" for i in range(n_norm)]
expr = pd.DataFrame(np.vstack([X, Xn]).T, index=genes, columns=tum_ids + norm_ids)
expr.index.name = "sample"
with gzip.open(C.DATA_RAW / C.FILES["expression"], "wt") as f:
    expr.round(4).to_csv(f, sep="\t")

stages = rng.choice(["Stage I", "Stage IA", "Stage IIA", "Stage IIB", "Stage IIIA",
                     "Stage IIIC", "Stage IV", "Stage X"], n_tum,
                    p=[0.1, 0.1, 0.3, 0.2, 0.12, 0.08, 0.05, 0.05])
age = rng.normal(58, 13, n_tum).round()
pam = y.astype(object)
pam[rng.random(n_tum) < 0.05] = np.nan
clin = pd.DataFrame({"PAM50Call_RNAseq": pam, "age_at_initial_pathologic_diagnosis": age,
                     "pathologic_stage": stages}, index=pd.Index(tum_ids, name="sampleID"))
clin = pd.concat([clin, pd.DataFrame(index=pd.Index(norm_ids, name="sampleID"))])
clin.to_csv(C.DATA_RAW / C.FILES["clinical"], sep="\t")

risk = (pd.Series(y).map({"LumA": 0, "LumB": 0.6, "Basal": 0.9, "Her2": 0.8, "Normal": 0.2}).values
        + pd.Series(stages).str.contains("III|IV").values * 0.9 + (age - 58) / 40)
surv = pd.DataFrame(index=pd.Index(tum_ids, name="sample"))
for ep, scale in [("PFI", 3000), ("OS", 5000)]:
    event_time = rng.exponential(scale * np.exp(-risk))
    censor = rng.uniform(100, 4000, n_tum)
    surv[ep] = (event_time <= censor).astype(int)
    surv[f"{ep}.time"] = np.minimum(event_time, censor).round()
surv.to_csv(C.DATA_RAW / C.FILES["survival"], sep="\t")
print(f"Synthetic data written to {C.DATA_RAW}")
