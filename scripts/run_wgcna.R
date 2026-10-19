
suppressMessages({library(WGCNA); library(survival)})
options(stringsAsFactors=FALSE); allowWGCNAThreads(8)
set.seed(1)
e <- read.csv("wgcna/tcga_top5000_expr.csv", row.names=1, check.names=FALSE)
datExpr <- as.data.frame(t(e))
gsg <- goodSamplesGenes(datExpr, verbose=0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
tree <- hclust(dist(datExpr), method="average")
cutH <- quantile(tree$height,0.99)*1.15
cut <- cutreeStatic(tree, cutHeight=cutH, minSize=10)
keep <- cut==as.integer(names(which.max(table(cut))))
removed <- rownames(datExpr)[!keep]
datExpr <- datExpr[keep, ]
cat("samples kept:", nrow(datExpr), " removed:", length(removed), "\n")
sft <- pickSoftThreshold(datExpr, powerVector=c(1:20,seq(22,30,2)), networkType="signed", RsquaredCut=0.85, verbose=0)
fi <- sft$fitIndices
fi$signedR2 <- -sign(fi$slope)*fi$SFT.R.sq
write.csv(fi, "wgcna/wgcna_soft_threshold.csv", row.names=FALSE)
beta <- sft$powerEstimate
if (is.na(beta)) { cand <- fi$Power[fi$signedR2 > 0.85]; beta <- if(length(cand)) min(cand) else 30 }
cat("beta=", beta, " R2=", fi$SFT.R.sq[fi$Power==beta], "\n")
net <- blockwiseModules(datExpr, power=beta, networkType="signed", TOMType="signed",
        minModuleSize=30, mergeCutHeight=0.25, numericLabels=FALSE, pamRespectsDendro=FALSE,
        maxBlockSize=ncol(datExpr), deepSplit=2, reassignThreshold=0, verbose=0)
MEs <- orderMEs(moduleEigengenes(datExpr, net$colors)$eigengenes)
tr <- read.csv("wgcna/tcga_traits.csv", row.names=1)
tr <- tr[rownames(datExpr), , drop=FALSE]
# OS trait: also use martingale residual of null Cox and log-rank-based ME association
cox_rows <- do.call(rbind, lapply(colnames(MEs), function(m){
  f <- coxph(Surv(tr$OS_time, tr$OS_event) ~ scale(MEs[[m]]))
  s <- summary(f); data.frame(module=m, cox_HR_perSD=s$conf.int[1,1], cox_p=s$coefficients[1,5])}))
cor_mt <- cor(MEs, tr, use="pairwise.complete.obs")
n_mt <- crossprod(!is.na(as.matrix(tr)) * 1, !is.na(as.matrix(MEs)) * 1)  # not needed; use per-column n
p_mt <- cor_mt
for (j in colnames(tr)) { nj <- sum(!is.na(tr[[j]])); p_mt[,j] <- corPvalueStudent(cor_mt[,j], nj) }
long <- do.call(rbind, lapply(colnames(tr), function(j) data.frame(module=rownames(cor_mt), trait=j, cor=cor_mt[,j], p=p_mt[,j], n=sum(!is.na(tr[[j]])))))
long <- merge(long, cox_rows, by="module")
long$module_size <- as.integer(table(net$colors)[sub("^ME","",long$module)])
write.csv(long, "wgcna/wgcna_module_trait.csv", row.names=FALSE)
kME <- signedKME(datExpr, MEs)
ga <- data.frame(gene=colnames(datExpr), module=net$colors)
ga$kME_own <- sapply(seq_len(nrow(ga)), function(i) kME[i, paste0("kME", ga$module[i])])
write.csv(cbind(ga, kME), "wgcna/wgcna_module_assignment.csv", row.names=FALSE)
hub <- do.call(rbind, lapply(split(ga, ga$module), function(d){ d <- d[order(-d$kME_own),]; d$rank <- seq_len(nrow(d)); head(d, 30)}))
write.csv(hub, "wgcna/wgcna_hub_genes.csv", row.names=FALSE)
write.csv(data.frame(sample=rownames(MEs), MEs), "wgcna/wgcna_module_eigengenes.csv", row.names=FALSE)
writeLines(c(paste0("beta=",beta), paste0("R2=",fi$SFT.R.sq[fi$Power==beta]), paste0("samples=",nrow(datExpr)), paste0("removed_outliers=",paste(removed,collapse=";")), paste0("n_modules=",length(unique(net$colors)))), "wgcna/wgcna_run_info.txt")
print(table(net$colors))
