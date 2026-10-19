# Reproduces predicted_ic50_<cohort>.csv and gdsc2_cv5_quality.csv (track "Drugs").
# oncoPredict 1.3.1 CALCPHENOTYPE.R sourced from the CRAN tarball (see drugs_README.md, deviation 1).
suppressPackageStartupMessages({library(ridge); library(car); library(sva); library(preprocessCore); library(arrow); library(glmnet); library(pls)})
source("oncopredict_src/oncoPredict/R/CALCPHENOTYPE.R")
trE <- readRDS("gdsc2/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")  # 17419 x 805
trR <- readRDS("gdsc2/DataFiles/DataFiles/Training Data/GDSC2_Res.rds")                                          # 805 x 198 (ln IC50)
load_test <- function(path){ d <- as.data.frame(read_parquet(path)); ic <- intersect(c("sample","Gene_Name"), names(d))[1]
  rn <- d[[ic]]; d[[ic]] <- NULL; m <- as.matrix(d); rownames(m) <- rn; m }
dir.create("op_output", showWarnings=FALSE)
for (co in c("TCGA","CGGA693","CGGA325")) {
  te <- load_test(sprintf("op_input/%s_expr_primary.parquet", co))   # genes x primary samples, log2(RSEM+1)
  set.seed(2026)
  pr <- calcPhenotype(trainingExprData=trE, trainingPtype=trR, testExprData=te, batchCorrect="eb",
                      powerTransformPhenotype=TRUE, removeLowVaryingGenes=0.2, minNumSamples=10, selection=1,
                      printOutput=FALSE, pcr=FALSE, removeLowVaringGenesFrom="rawData", report_pc=FALSE, cc=FALSE,
                      rsq=FALSE, folder=FALSE, parallel=TRUE, cores=10)
  write.csv(data.frame(sample=rownames(pr), pr, check.names=FALSE), sprintf("op_output/predicted_ic50_%s.csv", co), row.names=FALSE)
}
# 5-fold CV on the GDSC2 training set (quality floor)
set.seed(2026)
cv <- predictionAccuracybyCV(trainingExprData=trE, trainingPtype=trR, testExprData=NULL, batchCorrect="none",
            powerTransformPhenotype=TRUE, removeLowVaryingGenes=0.2, minNumSamples=10, selection=1, printOutput=FALSE,
            pcr=FALSE, removeLowVaringGenesFrom="rawData", cvFold=5, parallel=TRUE, cores=5)
stat <- function(d, f){ y <- cv$realPtype[,d]; p <- cv$cvPtype[,d]; ok <- !is.na(y) & !is.na(p); f(p[ok], y[ok]) }
cvdf <- data.frame(drug=colnames(trR), n_cell_lines=colSums(!is.na(trR)),
                   cv_R2=sapply(colnames(trR), stat, f=function(p,y) 1 - sum((p-y)^2)/sum((y-mean(y))^2)),
                   cv_pearson=sapply(colnames(trR), stat, f=function(p,y) cor(p,y)),
                   cv_spearman=sapply(colnames(trR), stat, f=function(p,y) cor(p,y,method="spearman")))
write.csv(cvdf, "op_output/gdsc2_cv5_quality.csv", row.names=FALSE)
