"""Run the full pipeline.

    python run_all.py          # real TCGA data (downloads ~100 MB first)
    python run_all.py --test   # quick check with small synthetic data
"""
import os
import subprocess
import sys
import time

STEPS = [
    "scripts/01_download_data.py",
    "scripts/02_preprocess.py",
    "scripts/03_unsupervised_subtypes.py",
    "scripts/04_survival_analysis.py",
    "scripts/05_predict_subtype.py",
    "scripts/06_survival_prediction.py",
]

env = dict(os.environ)
if "--test" in sys.argv:
    env["BRCA_TEST_MODE"] = "1"
    STEPS[0] = "tests/make_synthetic_data.py"
    print("TEST MODE: using synthetic data (outputs go to *_test folders)\n")

for step in STEPS:
    print(f"\n{'=' * 70}\n{step}\n{'=' * 70}", flush=True)
    start = time.time()
    result = subprocess.run([sys.executable, step], env=env)
    if result.returncode != 0:
        sys.exit(f"\nStopped: {step} failed. Fix the error above and rerun.")
    print(f"[{step} finished in {time.time() - start:.0f} s]")

print("\nAll steps finished. See figures/ and results/.")
