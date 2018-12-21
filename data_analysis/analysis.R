library(SuperLearner)
library(hal9001)
library(tmle)
library(dplyr)
library(drtmle)

code_dir <- "./"
data_dir <- "./"
save_dir <- "./"
# source in functions
source(paste0(code_dir,"oat_sl2.R"))
# load in data 
load(paste0(data_dir,"cebu_data_new.RData"))

# set a seed
set.seed(1234)

# call oat_sl to get nuisance estimators, saving files along the way
nuisance <- oat_sl(W = cebu_data$W, A = cebu_data$A, Y = cebu_data$Y,
                   full_adaptive_cv = TRUE, save_files = TRUE, 
                   save_dir = save_dir)

# save nuisance results
save(nuisance, file = paste0(save_dir, "nuisance.RData"))

## pump results into drtmle
# this fit will yield the standard tmle and one-step estimators/ci's
tmle_std <- drtmle::drtmle(
	Y = cebu_data$Y, A = cebu_data$A, W = cebu_data$W,
	# the next line are options native to the drtmle estimator
	# which we ignore. these options are chosen to make the code
	# run fast through the sections that pertain to the drtmle estimator
	SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
	# the treatment values of interest
	a_0 = c(0, 1, 2),
	# the outcome regression estimators formatted to match
	# the ordering of a_0
	Qn = list(nuisance$Q0W, nuisance$Q1W, nuisance$Q2W),
	# the propensity score estimators formatted to match
	# the ordering of a_0
	gn = list(nuisance$G0W, nuisance$G1W, nuisance$G2W),
	tolg = 1e-5
)

# this fit will yield cvtmle and cv one-step estimators/ci's
tmle_std_cv <- drtmle::drtmle(
	Y = cebu_data$Y, A = cebu_data$A, W = cebu_data$W,
	SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
	a_0 = c(0, 1, 2),
	# notice that we pump in CV versions of the nuisance parameters
	Qn = list(nuisance$cv_Q0W, nuisance$cv_Q1W, nuisance$cv_Q2W),
	gn = list(nuisance$cv_G0W, nuisance$cv_G1W, nuisance$cv_G2W),
	# actually not sure if this is necessary, but it couldn't hurt
	cvFolds = nuisance$fold_vec,
	tolg = 1e-5
)

# now we pump in the reduced dimension propensity score to get
# OATMLE and OA-one-step
oat_tmle_std <- drtmle::drtmle(
	Y = cebu_data$Y, A = cebu_data$A, W = cebu_data$W,
	SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
	a_0 = c(0, 1, 2),
	Qn = list(nuisance$Q0W, nuisance$Q1W, nuisance$Q2W),
	# note we now use GQW
	gn = list(nuisance$GQ0W, nuisance$GQ1W, nuisance$GQ2W),
	tolg = 1e-5
)
# now we pump in the "cross-validated" versions of GQW
oat_tmle_std_cv <- drtmle::drtmle(
	Y = cebu_data$Y, A = cebu_data$A, W = cebu_data$W,
	SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
	a_0 = c(0, 1, 2),
	# notice that we pump in CV versions of the nuisance parameters
	Qn = list(nuisance$cv_Q0W, nuisance$cv_Q1W, nuisance$cv_Q2W),
	gn = list(nuisance$cv_GQ0W, nuisance$cv_GQ1W, nuisance$cv_GQ2W),
	cvFolds = nuisance$fold_vec,
	tolg = 1e-5
)

# confidence intervals
ci_tmle <- tmle_std$tmle$est + t(c(0, -1.96, 1.96) %o% sqrt(diag(tmle_std_cv$tmle$cov)))
ci_ctmle <- oat_tmle_std$tmle$est + t(c(0, -1.96, 1.96) %o% sqrt(diag(oat_tmle_std_cv$tmle$cov)))
ci_os <- tmle_std$aiptw$est + t(c(0, -1.96, 1.96) %o% sqrt(diag(tmle_std_cv$aiptw$cov)))
ci_cos <- oat_tmle_std$aiptw$est + t(c(0, -1.96, 1.96) %o% sqrt(diag(oat_tmle_std_cv$aiptw$cov)))

# hypothesis tests
one_wald_test <- function(est_mat, cov_mat){
	# transformation matrix
	A <- matrix(c(1,-1,0,
	              1,0,-1), nrow = 2, byrow = TRUE)
	# covariance matrix for wald test
	Sigma <- A %*% cov_mat %*% t(A)
	diff_mat <- A %*% est_mat
	library(MASS)
	# wald statistic
	T <- t(diff_mat) %*% solve(Sigma) %*% diff_mat 
	# chi-square p-value
	pval <- pchisq(T, df = 2, lower.tail = FALSE)
	return(pval)
}

# p-values
pval_tmle <- one_wald_test(est_mat = matrix(tmle_std$tmle$est),
                           cov_mat = tmle_std_cv$tmle$cov)
pval_ctmle <- one_wald_test(est_mat = matrix(oat_tmle_std$tmle$est),
                           cov_mat = oat_tmle_std_cv$tmle$cov)                          
pval_os <- one_wald_test(est_mat = matrix(tmle_std$aiptw$est),
                           cov_mat = tmle_std_cv$aiptw$cov)
pval_cos <- one_wald_test(est_mat = matrix(oat_tmle_std$aiptw$est),
                           cov_mat = oat_tmle_std_cv$aiptw$cov)

# glm results
fit <- glm(cebu_data$Y ~ ., data = data.frame(A1 = as.numeric(cebu_data$A == 1), 
est_glm <- c(mean(predict(fit, newdata = data.frame(A1 = 0, A2 = 0, cebu_data$W))),
             mean(predict(fit, newdata = data.frame(A1 = 1, A2 = 0, cebu_data$W))),
             mean(predict(fit, newdata = data.frame(A1 = 0, A2 = 2, cebu_data$W))))             
                                              A2 = as.numeric(cebu_data$A == 2), cebu_data$W))

# bootstrap confidence intervals
one_boot <- function(cebu_data, n = length(cebu_data$W[,1])){
	samp_idx <- sample(1:n, replace = TRUE)
	boot_dat <- list()
	boot_dat$W <- cebu_data$W[samp_idx,]
	boot_dat$A <- cebu_data$A[samp_idx]
	boot_dat$Y <- cebu_data$Y[samp_idx]
	# glm results
	boot_fit <- glm(boot_dat$Y ~ ., data = data.frame(A1 = as.numeric(boot_dat$A == 1), 
	                                              A2 = as.numeric(boot_dat$A == 2), boot_dat$W))
	return(c(mean(predict(boot_fit, newdata = data.frame(A1 = 0, A2 = 0, boot_dat$W))),
             mean(predict(boot_fit, newdata = data.frame(A1 = 1, A2 = 0, boot_dat$W))),
             mean(predict(boot_fit, newdata = data.frame(A1 = 0, A2 = 2, boot_dat$W)))))
}
boot_samples <- replicate(500, one_boot(cebu_data = cebu_data))
ci_glm <- cbind(est_glm, t(apply(boot_samples, 1, quantile, p = c(0.025, 0.975))))

# p-value based on wald test
library(sandwich)
fit_cov_mat <- vcovHC(fit, type = "HC0")[c("A1","A2"),c("A1","A2")]
fit_coef <- coef(fit)
fit_est_mat <- matrix(fit_coef[c("A1","A2")], ncol = 1)
wald_stat <- t(fit_est_mat) %*% solve(fit_cov_mat) %*% fit_est_mat
pval_glm <- pchisq(wald_stat, df = 2, lower.tail = FALSE)


# make results table
one_row <- function(ci, pval){
	ci_pretty <- t(apply(ci, 1, formatC, format = "f", digits = 2))
	pval_pretty <- formatC(pval, digits = 2)
	collapse_ci <- function(one_row){
		paste0(one_row[1], " (", one_row[2], ", ", one_row[3],")")
	}
	ci_collapsed <- c(apply(ci_pretty, 1, collapse_ci)[c(2,1,3)], pval_pretty)
	return(ci_collapsed)
}

full_results <- rbind(
  one_row(ci_tmle, pval_tmle),
  one_row(ci_ctmle, pval_ctmle),
  one_row(ci_os, pval_os),
  one_row(ci_cos, pval_cos),
  one_row(ci_glm, pval_glm)
)
full_results <- cbind(c("TMLE","CTMLE","OS", "COS", "LM"), full_results)
colnames(full_results) <- c("Method", "Pre-term", "Full-term", "Post-term", "p-value")

# format for latex
x_results <- xtable(full_results)
print(x_results, include.rownames = FALSE)

#------------------------------------------
# for appendix results on super learner
#------------------------------------------
# outcome regression 
load(paste0(save_dir, "sl_or.RData"))
# format risks/coef for latex
# xtable(print(sl_or), digits = 3)

# cv outcome regression 
load(paste0(save_dir, "cv_sl_or.RData"))
# re-name method so summary.CV.SuperLEarner works
cv_sl_or$call[["method"]] <- "method.CC_LS"
# format risks/coef for latex
tmp <- print(summary(cv_sl_or))
tmp_extra <- paste0("(",formatC(tmp$Min,format = "f",digits = 3) ,", ", formatC(tmp$Max, format = "f", digits = 3),")")
tmp <- cbind(tmp[,-which(colnames(tmp) %in% c("Min","Max"))], tmp_extra)
tmp <- cbind(tmp, c("--", "--", formatC(sl_or$coef, format = "f", digits = 3)))
colnames(tmp) <- c("Algorithm", "Avg. risk", "SE risk)", "Range risk", "SL coef.")
x_tmp <- xtable(tmp, digits = 3)
print(x_tmp, include.rownames = FALSE)

# propensity score fits
load(paste0(save_dir, "sl_ps.RData"))
load(paste0(save_dir, "sl_ps2.RData"))

xtable(print(sl_ps), digits = 3)
xtable(print(sl_ps2), digits = 3)

# cv propensity score models
load(paste0(save_dir, "cv_sl_ps.RData"))
load(paste0(save_dir, "cv_sl_ps2.RData"))
# re-name method so summary.CV.SuperLEarner works
cv_sl_ps$call[["method"]] <- "method.CC_nloglik"
cv_sl_ps2$call[["method"]] <- "method.CC_nloglik"

# format risks/coef for latex
tmp <- print(summary(cv_sl_ps))
tmp_extra <- paste0("(",formatC(tmp$Min,format = "f",digits = 3) ,", ", formatC(tmp$Max, format = "f", digits = 3),")")
tmp <- cbind(tmp[,-which(colnames(tmp) %in% c("se","Min","Max"))], tmp_extra)
colnames(tmp) <- c("Algorithm", "Avg. risk", "Range risk")
x_tmp <- xtable(tmp, digits = 3)
print(x_tmp, include.rownames = FALSE)
# format risks/coef for latex
tmp <- print(summary(cv_sl_ps2))
tmp_extra <- paste0("(",formatC(tmp$Min,format = "f",digits = 3) ,", ", formatC(tmp$Max, format = "f", digits = 3),")")
tmp <- cbind(tmp[,-which(colnames(tmp) %in% c("se","Min","Max"))], tmp_extra)
colnames(tmp) <- c("Algorithm", "Avg. risk", "Range risk")
x_tmp <- xtable(tmp, digits = 3)
print(x_tmp, include.rownames = FALSE)




