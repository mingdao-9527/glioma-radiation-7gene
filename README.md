# Radiation-response 7-gene score in diffuse glioma — analysis code

Code accompanying: *A radiation-response-gene-derived 7-gene score refines molecular prognostic
stratification of diffuse glioma across eight cohorts but does not predict radiotherapy benefit:
an integrative multi-omics and randomised-cohort analysis* (Guo J, Qu L, Wang Z, Zhang G, Zhang G; Department of
Radiation Oncology and Center for Neurodegenerative Diseases, the Second People's Hospital of Liaocheng,
Liaocheng, Shandong, China).

## Apply the signature to your own cohort
```
python scripts/compute_risk_score.py --expr my_expression.csv --out risk_scores.csv
```
`my_expression.csv`: genes x samples, log2 expression (use `--log2` for raw counts/TPM).
Coefficients (per within-cohort z-score of log2 expression) are in `data/signature_coefficients.csv`:
SMC4 +0.1905, E2F7 +0.2052, PSMC2 +0.1772, HIST1H2BH/H2BC9 +0.1279, HOXC10 +0.1224, BMP2 -0.2444, MMS19 -0.0841.
The score is cohort-relative (z-scored within cohort); the high/low cut is the within-cohort median.

## Pipeline map (step -> script)
| Step | Script | Notes |
|---|---|---|
| Cohort download / harmonisation (TCGA via UCSC Xena, CGGA) | run as notebook cells (see Supplementary Methods) | sources and URLs in the manuscript Data availability |
| WGCNA (signed, beta = 17) | `scripts/run_wgcna.R` | WGCNA compiled from CRAN source |
| 101-combination ML grid, validation, RT x risk interaction | `scripts/sig_lib.py` | Cox/C-index/time-AUC/KM/DCA helpers and the grid definition |
| Nomogram (rms) | `scripts/nomogram.R` | |
| Drug sensitivity (oncoPredict / GDSC2) | `scripts/run_oncopredict_gdsc2.R`, `scripts/CALCPHENOTYPE.R` | `CALCPHENOTYPE.R` is the unmodified function from CRAN oncoPredict 1.3.1 (GPL) |
| IDH-mutant grade 2-3 subgroup analysis | `scripts/idhmut_analysis.py` | |
| WHO 2021 molecular grading, PFI/DSS, CGGA + TMZ Cox | `scripts/secondary.py` | |
| Additional validation cohorts (EORTC 22033, REMBRANDT RMA, GSE43378, GSE4412, CGGA-301) | `scripts/extracohorts_pipeline_code.py` | includes the R RMA snippet |
| Figures | `scripts/figlib.py` | palette, fixed-canvas panel export, grid assembler |

Consensus clustering, immune (ssGSEA/ESTIMATE), genomics/pan-cancer, single-cell and cis-MR
steps were executed as interactive notebook cells; their parameters, sources and outputs are
documented in the Supplementary Methods and Supplementary Tables S1-S27.

## Environment
See `environment.yml` (Python 3.11: pandas, numpy, scipy, lifelines, scikit-survival, scikit-learn,
statsmodels, gseapy, matplotlib; R 4.4/4.5: WGCNA, rms, survival, ridge, oncoPredict dependencies).

## License
MIT (see LICENSE). `scripts/CALCPHENOTYPE.R` retains its original GPL licence from oncoPredict.
