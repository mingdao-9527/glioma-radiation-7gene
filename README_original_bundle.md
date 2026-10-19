# Code bundle — analysis scripts

Scripts recovered from the project artifact store (no data included). Every script expects the processed cohort files described in the manuscript Supplementary Methods (cohorts_README.md) in its working directory.

| Analysis step | Script |
|---|---|
| Signature track: feature filtering, 101-combination selector × model grid, random-gene controls, benchmark (RSI, Kim 2012, IDH, grade), validation metrics, RT × risk interaction (Python). | `sig_lib.py` |
| Signed WGCNA on TCGA (soft threshold, blockwiseModules, module–trait correlation, hub genes) (R). | `run_wgcna.R` |
| rms nomogram of score + age + grade + IDH in irradiated TCGA, bootstrap validation, external application to CGGA (R). | `nomogram.R` |
| oncoPredict calcPhenotype on GDSC2 for TCGA/CGGA, 5-fold CV quality floor (R; sources CALCPHENOTYPE.R). | `run_oncopredict_gdsc2.R` |
| oncoPredict 1.3.1 calcPhenotype source (third-party, CRAN/GPL) sourced by run_oncopredict_gdsc2.R because the package could not be installed. | `CALCPHENOTYPE.R` |
| Revision: CGGA multivariable Cox with TMZ and TCGA secondary endpoints (PFI/DSS); WHO 2021-adjusted models (Python). | `secondary.py` |
| Revision: IDH-mutant grade 2–3 subgroup — validation, RT × score interaction (per cohort, pooled, IPTW), power, de novo 16-gene grid and random-gene nulls (Python). | `idhmut_analysis.py` |
| Revision: screening and processing of GSE107850 (EORTC 22033), GSE108474 (REMBRANDT; RMA from CEL), GSE43378, GSE4412 and CGGA-301; metrics, KM, interaction tests (Python; RMA step calls R affy). | `extracohorts_pipeline_code.py` |
| Figures track: shared plotting helpers (panel labels, per-panel export, composites). | `figlib.py` |

Not recoverable as standalone scripts (run as notebook cells inside the analysis sessions; their outputs are the Supplementary Tables): cohort download/harmonisation (Cohorts track), curated gene pool and Cox screening (RadiationGenes track), consensus clustering (Subtyping track), ssGSEA/ESTIMATE/immune analyses (Immune track), mutation/CNV/pan-cancer (Genomics track), GBmap single-cell scoring (SingleCell track), cis-MR/coloc (Causal track), CNS WHO 2021 reclassification and GDC treatment retrieval (MolecularGrade track; the Cox/endpoint part is in secondary.py).

Software: Python 3.11–3.13 (pandas, numpy, scipy, statsmodels, scikit-learn 1.8, scikit-survival 0.27, lifelines 0.30, gseapy 1.3.1, scanpy 1.11.5), R 4.4–4.5 (WGCNA 1.73, rms, coloc 5.2.3, affy 1.84.0, ridge, sva, preprocessCore).
