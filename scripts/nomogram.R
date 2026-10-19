suppressMessages({library(rms); library(survival)})
d <- read.csv("sig/risk_scores_all_cohorts.csv", stringsAsFactors=FALSE)
d <- d[!is.na(d$grade) & !is.na(d$IDH_status) & !is.na(d$age),]
d$grade <- factor(d$grade, levels=c("G2","G3","G4")); d$IDH <- factor(ifelse(d$IDH_status=="Wildtype","Wildtype","Mutant"), levels=c("Mutant","Wildtype"))
d$OS.years <- d$OS.time/365.25
tr <- d[d$cohort=="TCGA" & d$RT=="yes" & !is.na(d$RT),]
dd <- datadist(tr); options(datadist="dd")
f <- cph(Surv(OS.years, OS) ~ risk + age + grade + IDH, data=tr, x=TRUE, y=TRUE, surv=TRUE, time.inc=3)
print(f)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1,x), function(x) surv(3,x), function(x) surv(5,x)), funlabel=c("1-year OS","3-year OS","5-year OS"), lp=TRUE)
# export point tables
rows <- list()
for (v in c("risk","age","grade","IDH")) { tab <- nom[[v]]; rows[[v]] <- data.frame(variable=v, value=as.character(tab[[1]]), points=tab$points) }
pt <- do.call(rbind, rows)
tp <- data.frame(variable="total_points", value=as.character(nom$total.points$x), points=nom$total.points$x, lp=nom$lp$x[match(nom$total.points$x, nom$total.points$x)])
tot <- nom$total.points$x; lpv <- nom$lp
# map total points -> lp -> survival
lp_tab <- data.frame(variable="lp_from_total_points", total_points=lpv$x, lp=lpv$x.real)
lp_tab$surv1 <- surv(1, lp_tab$lp); lp_tab$surv3 <- surv(3, lp_tab$lp); lp_tab$surv5 <- surv(5, lp_tab$lp)
write.csv(pt, "sig/nomogram_points.csv", row.names=FALSE)
write.csv(lp_tab, "sig/nomogram_total_points_to_survival.csv", row.names=FALSE)
coefs <- data.frame(term=names(coef(f)), coef=coef(f), se=sqrt(diag(vcov(f))))
write.csv(coefs, "sig/nomogram_model_coefficients.csv", row.names=FALSE)
# internal validation (bootstrap 200) and calibration at 3y
set.seed(1); v <- validate(f, B=200); write.csv(as.data.frame(unclass(v)), "sig/nomogram_internal_validation.csv")
cal_rows <- list()
for (yr in c(1,3,5)) {
  fy <- cph(Surv(OS.years, OS) ~ risk + age + grade + IDH, data=tr, x=TRUE, y=TRUE, surv=TRUE, time.inc=yr)
  cal <- tryCatch(calibrate(fy, u=yr, m=60, B=200, cmethod="KM"), error=function(e) NULL)
  if (!is.null(cal)) cal_rows[[length(cal_rows)+1]] <- data.frame(cohort="TCGA_RT", year=yr, mean_predicted=cal[,"mean.predicted"], KM=cal[,"KM"], KM_corrected=cal[,"KM.corrected"], std_err=cal[,"std.err"])
}
# external: linear predictor applied to CGGA, calibration (KM observed by predicted-tertile/quintile) and C-index
ext_rows <- list(); cidx <- list()
for (co in c("TCGA","CGGA693","CGGA325")) for (popn in c("RT","all")) {
  te <- d[d$cohort==co,]; if (popn=="RT") te <- te[!is.na(te$RT) & te$RT=="yes",]
  te$lp <- predict(f, newdata=te, type="lp")
  cc <- rcorr.cens(-te$lp, Surv(te$OS.years, te$OS)); cidx[[length(cidx)+1]] <- data.frame(cohort=co, population=popn, n=nrow(te), C_index=cc["C Index"], se=cc["S.D."]/2)
  for (yr in c(1,3,5)) {
    te$pred <- surv(yr, te$lp); te$q <- cut(te$pred, quantile(te$pred, seq(0,1,0.2)), include.lowest=TRUE, labels=FALSE)
    for (q in sort(unique(te$q))) { s <- te[te$q==q,]; km <- survfit(Surv(OS.years, OS) ~ 1, data=s); sm <- summary(km, times=yr, extend=TRUE)
      ext_rows[[length(ext_rows)+1]] <- data.frame(cohort=co, population=popn, year=yr, quintile=q, n=nrow(s), mean_predicted=mean(s$pred), KM=sm$surv, KM_lo=sm$lower, KM_hi=sm$upper) }
  }
}
write.csv(do.call(rbind, cal_rows), "sig/nomogram_calibration_internal.csv", row.names=FALSE)
write.csv(do.call(rbind, ext_rows), "sig/nomogram_calibration_external.csv", row.names=FALSE)
write.csv(do.call(rbind, cidx), "sig/nomogram_cindex.csv", row.names=FALSE)
print(do.call(rbind, cidx)); print(v[1:3,])
