library(drtmle)
library(dplyr)
do_once <- function(dat, dat_name, permute_y = FALSE){
  if (permute_y) dat$Y <- sample(dat$Y)
  ## get nuisance estimates using oat_hal
  # output is a data.frame with nuisance estimators evaluated at observed
  # data and a vector of CV-folds
  nuisance <- oat_sl(
    W = dat$W, A = dat$A, Y = dat$Y,
    V = 5, lambda_seq = exp(seq(-1, -8, length = 1e2))
  )
  ## pump results into drtmle
  # this fit will yield the standard tmle and one-step estimators/ci's
  tmle_std <- drtmle::drtmle(
    Y = dat$Y, A = dat$A, W = dat$W,
    # the next line are options native to the drtmle estimator
    # which we ignore. these options are chosen to make the code
    # run fast through the sections that pertain to the drtmle estimator
    SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
    # the treatment values of interest
    a_0 = c(0, 1),
    # the outcome regression estimators formatted to match
    # the ordering of a_0
    Qn = list(nuisance$Q0W, nuisance$Q1W),
    # the propensity score estimators formatted to match
    # the ordering of a_0
    gn = list(1 - nuisance$G1W, nuisance$G1W),
    tolg = 1e-5
  )

  # this fit will yield cvtmle and cv one-step estimators/ci's
  tmle_std_cv <- drtmle::drtmle(
    Y = dat$Y, A = dat$A, W = dat$W,
    SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
    a_0 = c(0, 1),
    # notice that we pump in CV versions of the nuisance parameters
    Qn = list(nuisance$cv_Q0W, nuisance$cv_Q1W),
    gn = list(1 - nuisance$cv_G1W, nuisance$cv_G1W),
    # actually not sure if this is necessary, but it couldn't hurt
    cvFolds = nuisance$fold_vec,
    tolg = 1e-5
  )

  # now we pump in the reduced dimension propensity score to get
  # OATMLE and OA-one-step
  oat_tmle_std <- drtmle::drtmle(
    Y = dat$Y, A = dat$A, W = dat$W,
    SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
    a_0 = c(0, 1),
    Qn = list(nuisance$Q0W, nuisance$Q1W),
    # note we now use GQW
    gn = list(1 - nuisance$GQW, nuisance$GQW),
    tolg = 1e-5
  )
  # now we pump in the "cross-validated" versions of GQW
  oat_tmle_std_cv <- drtmle::drtmle(
    Y = dat$Y, A = dat$A, W = dat$W,
    SL_Qr = "SL.mean", SL_gr = "SL.mean", maxIter = 0,
    a_0 = c(0, 1),
    # notice that we pump in CV versions of the nuisance parameters
    Qn = list(nuisance$cv_Q0W, nuisance$cv_Q1W),
    gn = list(1 - nuisance$cv_GQW, nuisance$cv_GQW),
    cvFolds = nuisance$fold_vec,
    tolg = 1e-5
  )
  ## extract results using ci method
  # standard tmle + ci
  rslt_tmle_std <- drtmle::ci(tmle_std, est = "tmle", contrast = c(-1,1))$tmle[1, ]
  # standard one-step + ci
  rslt_aiptw_std <- drtmle::ci(tmle_std, est = "aiptw", contrast = c(-1,1))$aiptw[1, ]
  # cvtmle + ci (centered at cvtmle)
  rslt_cvtmle_std <- drtmle::ci(tmle_std_cv, est = "tmle", contrast = c(-1,1))$tmle[1, ]
  # cvone-step + ci (centered at one-step)
  rslt_cvaiptw_std <- drtmle::ci(tmle_std_cv, est = "aiptw", contrast = c(-1,1))$aiptw[1, ]
  # OATMLE + ci
  rslt_oat_tmle_std <- drtmle::ci(oat_tmle_std, est = "tmle", contrast = c(-1,1))$tmle[1, ]
  # OA onestep + ci
  rslt_oat_aiptw_std <- drtmle::ci(oat_tmle_std, est = "aiptw", contrast = c(-1,1))$aiptw[1, ]
  # cv OATMLE + ci (centered at cvOATMLE)
  rslt_cvoat_tmle_std <- drtmle::ci(oat_tmle_std_cv, est = "tmle", contrast = c(-1,1))$tmle[1, ]
  # cv OA one-step + ci (centered at cvOA one-step)
  rslt_cvoat_aiptw_std <- drtmle::ci(oat_tmle_std_cv, est = "aiptw", contrast = c(-1,1))$aiptw[1, ]

  ## build confidence intervals that are centered at the standard (not cv) estimators
  # ci for standard tmle with cv-variance estimates
  # here I've elected not to use the targeted nuisance estimators (i.e., getting cv-aiptw variance estimates
  # instead of cvtmle)
  rslt_tmle_std_cvci <- rslt_tmle_std[1] + c(-1,1) * diff(rslt_cvaiptw_std[2:3])/2
  # ci for standard aiptw with cv-variance estimates
  rslt_aiptw_std_cvci <- rslt_aiptw_std[1] + c(-1,1) * diff(rslt_cvaiptw_std[2:3])/2
  # ci for OATMLE with cv-variance estimators (as above, using untargeted nuisance estimators in the variance)
  rslt_oat_tmle_std_cvci <- rslt_oat_tmle_std[1] + c(-1,1) * diff(rslt_cvoat_aiptw_std[2:3])/2
  # ci for OA one-step with cv-variance estimators
  rslt_oat_aiptw_std_cvci <- rslt_oat_aiptw_std[1] + c(-1,1) * diff(rslt_cvoat_aiptw_std[2:3])/2

  ## format output
  # could include cvtmle and cvone-step in here as well

  # diff(rslt_cvtmle_std[2:3]) / 1.96 / 2
  # diff(rslt_cvaiptw_std[2:3]) / 1.96 / 2
  # diff(rslt_cvoat_tmle_std[2:3]) / 1.96 / 2
  # diff(rslt_cvoat_aiptw_std[2:3]) / 1.96 / 2

  df_result <- list()
  df_result <- c(
    df_result,
    list(data.frame(
      method = "tmle",
      Psi = rslt_tmle_std[1],
      sd = diff(rslt_tmle_std[2:3]) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_tmle_std[2], rslt_tmle_std[3])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "tmle_cv_variance",
      Psi = rslt_tmle_std[1],
      sd = diff(rslt_tmle_std_cvci) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_tmle_std_cvci[1], rslt_tmle_std_cvci[2])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "onestep",
      Psi = rslt_aiptw_std[1],
      sd = diff(rslt_aiptw_std[2:3]) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_aiptw_std[2], rslt_aiptw_std[3])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "onestep_cv_variance",
      Psi = rslt_aiptw_std[1],
      sd = diff(rslt_aiptw_std_cvci) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_aiptw_std_cvci[1], rslt_aiptw_std_cvci[2])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "oat",
      Psi = rslt_oat_tmle_std[1],
      sd = diff(rslt_oat_tmle_std[2:3]) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_oat_tmle_std[2], rslt_oat_tmle_std[3])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "oat_cv_variance",
      Psi = rslt_oat_tmle_std[1],
      sd = diff(rslt_oat_tmle_std_cvci) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_oat_tmle_std_cvci[1], rslt_oat_tmle_std_cvci[2])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "onestep_oat",
      Psi = rslt_oat_aiptw_std[1],
      sd = diff(rslt_oat_aiptw_std[2:3]) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_oat_aiptw_std[2], rslt_oat_aiptw_std[3])
    ))
  )
  df_result <- c(
    df_result,
    list(data.frame(
      method = "onestep_oat_cv_variance",
      Psi = rslt_oat_aiptw_std[1],
      sd = diff(rslt_oat_aiptw_std_cvci) / 1.96 / 2,
      is_cover = between(Psi_0, rslt_oat_aiptw_std_cvci[1], rslt_oat_aiptw_std_cvci[2])
    ))
  )
  psi_plugin <- mean(nuisance$Q1W - nuisance$Q0W)
  df_result <- c(
    df_result,
    list(data.frame(
      method = "plugin",
      Psi = psi_plugin,
      sd = 999,
      is_cover = TRUE
    ))
  )

  df_result <- do.call(rbind, df_result)
  rownames(df_result) <- NULL
  df_result$bias <- df_result$Psi - Psi_0
  df_result$mse <- df_result$bias ^ 2
  df_result$n <- nrow(dat$W)
  df_result$n_w <- ncol(dat$W)
  df_result$dat_name <- dat_name
  df_result$it2 <- it2
  # return(list(df_result = df_result, Y = dat$Y))
  return(df_result)
}
