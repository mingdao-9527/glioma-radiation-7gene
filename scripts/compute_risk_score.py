#!/usr/bin/env python
"""Compute the 7-gene radiation-response risk score for a new glioma cohort.

Usage:
    python compute_risk_score.py --expr expression.csv --out risk_scores.csv [--log2]

expression.csv : genes x samples matrix (first column = gene symbol), log2-transformed
                 expression (log2(TPM+1), log2(FPKM+1) or log2(RSEM+1)); pass --log2 to
                 log-transform raw values. Microarray data should already be log2 RMA.
Procedure (identical to the manuscript):
    1. keep the 7 signature genes (HIST1H2BH is H2BC9 in current HGNC nomenclature);
    2. z-score each gene ACROSS SAMPLES WITHIN THE COHORT (mean 0, SD 1, ddof=0);
    3. risk = sum(coef_g * z_g); high/low risk = split at the within-cohort median.
The score is cohort-relative; do not apply z-scores from another cohort.
"""
import argparse, sys, os
import numpy as np, pandas as pd

ALIASES = {"H2BC9": "HIST1H2BH", "HIST1H2BH": "HIST1H2BH"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True)
    ap.add_argument("--out", default="risk_scores.csv")
    ap.add_argument("--log2", action="store_true", help="apply log2(x+1) to the input")
    ap.add_argument("--coef", default=os.path.join(os.path.dirname(__file__), "..", "data", "signature_coefficients.csv"))
    a = ap.parse_args()
    coef = pd.read_csv(a.coef).set_index("gene")["coef"]
    X = pd.read_csv(a.expr, index_col=0)
    X.index = [ALIASES.get(str(g).upper(), str(g).upper()) for g in X.index]
    X = X[~X.index.duplicated(keep="first")]
    if a.log2:
        X = np.log2(X.clip(lower=0) + 1)
    missing = [g for g in coef.index if g not in X.index]
    if missing:
        sys.exit("Missing signature genes in the expression matrix: %s" % missing)
    sub = X.loc[coef.index].astype(float)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0), axis=0)
    risk = (z.mul(coef, axis=0)).sum(axis=0)
    out = pd.DataFrame({"sample": risk.index, "risk_score": risk.values})
    out["risk_group"] = np.where(out.risk_score > out.risk_score.median(), "High", "Low")
    out.to_csv(a.out, index=False)
    print("Scored %d samples; median cut %.4f; written to %s" % (len(out), out.risk_score.median(), a.out))

if __name__ == "__main__":
    main()
