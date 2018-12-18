library(SuperLearner)
#' OAT SL
#'
#' Fit a custom version of super learner for outcome adaptive TMLE estimation
#'
#' @param W A data.frame of predictors
#' @param A A numeric binary (0/1) treatment vector
#' @param Y A numeric outcome vector
#' @param V The number of cross-validation folds (must be at least 3 or \code{glmnet} will complain)
#' @param outcome_family "gaussian" or "binomial" (implies loss function to use)
#' @param full_adaptive_cv If TRUE then the adaptive PS is tuned using cross-validated versions of Q (i.e.,
#' more formal cross-validation). If FALSE then the adaptive PS is tuned based on Q that was fit using
#' the full data.
#' @param SL.library Library to use for super learner
#' @param lambda_seq Sequence of lambdas passed to HAL fit for CTMLE propensity score
#'
#' @return A data.frame. \code{QAW} = the HAL fit at CV-selected lambda for Qbar evaluated at (A_i, W_i),
#' for i = 1,...,n ; \code{Q1W} = the HAL fit at CV-selected lambda for Qbar evaluated at (1, W_i) ;
#' \code{Q0W} = the HAL fit at CV-selected lambda for Qbar evaluated at (0, W_i) ;
#' \code{G1W} = the HAL fit at CV-selected lambda for G evaluated at W_i ;
#' \code{GQW} = the HAL fit of A ~ Q1W + Q0W evaluated at Q(1,W_i), Q(0,W_i) ;
#' \code{cv_} = the cross-validated HAL fit at CV-selected lambda of the above nuisance parameters
#' @importFrom hal9001 enumerate_basis make_design_matrix make_copy_map
oat_sl <- function(
  W, A, Y, V = 5, outcome_family = "gaussian",
  full_adaptive_cv = FALSE,
  SL.library = c(
    # "SL.mean", "SL.glm", "SL.step.forward", "SL.earth", "SL.randomForest", "SL.glmnet"
    "SL.glm", "SL.earth", "SL.randomForest", "SL.glmnet"
    # "SL.mean", "SL.glm", "SL.step.forward", "SL.earth", "SL.glmnet"
    # "SL.mean", "SL.glm", "SL.glmnet"
  ),
  lambda_seq = exp(seq(-1, -13, length = 500))
) {
    n <- length(Y)
    # make a vector of cv-folds
    chunk2 <- function(x, n) split(x, cut(
        sample(seq_along(x), size = length(x), replace = FALSE),
        n,
        labels = FALSE
    ))
    fold_idx <- chunk2(seq_len(n), V)
    fold_vec <- rep(seq_len(V), unlist(lapply(fold_idx, length)))
    # fold_idx <- unlist(fold_idx, use.names = FALSE)

    # fit a super learner for OR
    newX_or <- rbind(
      data.frame(W, A = A),
      data.frame(W, A = 1),
      data.frame(W, A = 0)
    )

    sl_or <- SuperLearner(Y = Y, X = data.frame(W, A = A),
                          newX = newX_or,
                          SL.library = SL.library, family = outcome_family,
                          cvControl = list(V = V, validRows = fold_idx),
                          method = ifelse(outcome_family == "gaussian",
                                          "method.CC_LS",
                                          "method.CC_nloglik"))

    # extract predictions
    idx_aw <- 1:n
    idx_1w <- (n+1):(2*n)
    idx_0w <- (2*n+1):(3*n)
    QAW <- sl_or$SL.predict[idx_aw]
    Q1W <- sl_or$SL.predict[idx_1w]
    Q0W <- sl_or$SL.predict[idx_0w]

    # cross-validated version
    cv_sl_or <- CV.SuperLearner(Y = Y, X = data.frame(W, A = A),
                          SL.library = SL.library, family = outcome_family,
                          innerCvControl = list(list(V = V),list(V = V),list(V = V),
                                           list(V = V),list(V = V)),
                          cvControl = list(V = V, validRows = fold_idx),
                          control = list(saveFitLibrary = TRUE),
                          method = ifelse(outcome_family == "gaussian",
                                          "method.CC_LS",
                                          "method.CC_nloglik"))
    cv_QAW <- as.numeric(cv_sl_or$SL.predict)
    cv_Q1W <- rep(NA, n)
    cv_Q0W <- rep(NA, n)
    for(v in 1:5){
        this_idx <- cv_sl_or$folds[[v]]
        cv_Q1W[this_idx] <-
            predict(cv_sl_or$AllSL[[v]], newdata = data.frame(W, A = 1)[this_idx,])$pred
        cv_Q0W[this_idx] <-
            predict(cv_sl_or$AllSL[[v]], newdata = data.frame(W, A = 0)[this_idx,])$pred
    }

    # fit a super learner for PS
    sl_ps <- SuperLearner(Y = A,
                          X = W,
                          newX = W,
                          SL.library = SL.library, family = binomial(),
                          cvControl = list(V = V, validRows = fold_idx),
                          method = "method.CC_nloglik")

    # extract predictions
    G1W <- as.numeric(sl_ps$SL.predict)

    # cross-validated version
    cv_sl_ps <- CV.SuperLearner(Y = A, X = W,
                          SL.library = SL.library, family = binomial(),
                          innerCvControl = list(list(V = V),list(V = V),list(V = V),
                                           list(V = V),list(V = V)),
                          cvControl = list(V = V, validRows = fold_idx),
                          control = list(saveFitLibrary = TRUE),
                          method = "method.CC_nloglik")
    cv_G1W <- as.numeric(cv_sl_ps$SL.predict)

    # HAL A ~ Q1W + Q0W
    # first check that some variation in outcome regressions
    sanity_check <- !(all(Q1W == Q0W) & (length(unique(Q1W)) == 1))

    # make data matrix Q1W, Q0W
    if(sanity_check){
        if(!full_adaptive_cv){
            X <- as.matrix(cbind(Q1W, Q0W))
            basis_list <- hal9001::enumerate_basis(X, degrees = NULL)
            x_basis <- hal9001:::make_design_matrix(X, basis_list)
            copy_map <- hal9001:::make_copy_map(x_basis)
            unique_columns <- as.numeric(names(copy_map))
            x_basis <- x_basis[, unique_columns]

            hal_lasso <- glmnet::cv.glmnet(x = x_basis, y = A, nfolds = V,
                    family = "binomial", lambda = lambda_seq, foldid = fold_vec,
                    keep = TRUE, standardize = FALSE)
            lambda_goat_idx <- which(hal_lasso$lambda == hal_lasso$lambda.min)

            # get predictions in whole sample
            GQW_cvselect <- predict(hal_lasso, s = "lambda.min", type = "response", newx = x_basis)

            # cross-validated predictions
            GQW_cv_cvselect <- hal_lasso$fit.preval[,lambda_goat_idx]
        }else{
            # cross-validation routine based on CV-Q
            cv_out <- sapply(seq_len(V), one_red_hal,
                             fold_vec = fold_vec,
                             Q1W = cv_Q1W,
                             Q0W = cv_Q0W,
                             A = A,
                             lambda_seq = lambda_seq,
                             simplify = FALSE)
            # outcome regression
            risks <- colMeans(Reduce("rbind",lapply(cv_out, "[[", "risk")))
            # take smallest lambda that has smallest risk
            lambda_idx <- which.min(risks)[1]

            # re-fit on full data
            full_hal <- one_red_hal(fold = NULL, fold_vec = NULL,
                         A = A, Q1W = Q1W, Q0W = Q0W,
                         lambda_seq = lambda_seq)

            # GQW at cv selected lambda
            GQW_cvselect <- full_hal$GQW[, lambda_idx]
            # cross-validated Q1W at cv selected lambda
            GQW_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
                x$GQW[, lambda_idx]
            }))
        }
    }else{
        GQW_cvselect <- GQW_cv_cvselect <- rep(mean(A), n)
    }

    # format outcome
    out <- data.frame(QAW = QAW,
                      Q1W = Q1W,
                      Q0W = Q0W,
                      G1W = G1W,
                      GQW = as.numeric(GQW_cvselect),
                      cv_QAW = cv_QAW,
                      cv_Q1W = cv_Q1W,
                      cv_Q0W = cv_Q0W,
                      cv_G1W = cv_G1W,
                      cv_GQW = as.numeric(GQW_cv_cvselect),
                      fold_vec = fold_vec)
    return(out)
}


one_red_hal <- function(fold, fold_vec, outcome_family = "binomial",
                        A, Q1W, Q0W, lambda_seq = exp(seq(-1,-13,length=10000))){
    n <- length(A)
    n_valid <- sum(fold == fold_vec)
    n_train <- n - n_valid

    if(!is.null(fold)){
        x_fit <- cbind(Q1W, Q0W)[fold_vec != fold, ]
        y_fit <- A[fold_vec != fold]
    # if called fitting to full data
    }else{
        x_fit <- cbind(Q1W, Q0W)
        y_fit <- A
    }

    # for outcome regression
    basis_list <- hal9001::enumerate_basis(x_fit, NULL)
    x_basis <- hal9001:::make_design_matrix(x_fit, basis_list)
    copy_map <- hal9001:::make_copy_map(x_basis)
    unique_columns <- as.numeric(names(copy_map))
    # subset to non-duplicated columns
    x_basis_fit <- x_basis[, unique_columns]

    hal_lasso <- glmnet::glmnet(x = x_basis_fit, y = y_fit,
    family = outcome_family, lambda = lambda_seq,
    standardize = FALSE)

    # predictions on validation sample
    if(!is.null(fold)){
        new_x_fit1 <- cbind(Q1W, Q0W)[fold_vec == fold, ]
    }else{
        new_x_fit1 <- cbind(Q1W, Q0W)
    }

    # make HAL design for getting predictions to select lambda
    new_x_basis1 <- hal9001:::make_design_matrix(new_x_fit1, basis_list)
    new_x_basis1 <- as.matrix(new_x_basis1[, unique_columns])

    # get predictions
    beta_hat <- as.matrix(hal_lasso$beta)
    # intercept
    alpha_hat <- hal_lasso$a0
    # prediction matrices
    pred_matrix1 <- cbind(rep(1, ifelse(is.null(fold), n, n_valid)), new_x_basis1) %*% rbind(alpha_hat, beta_hat)
    if(outcome_family == "binomial"){
        pred_matrix1 <- apply(pred_matrix1, 2, plogis)
    }
    risk <- apply(pred_matrix1, 2, function(x){
            mean(ifelse(A[fold_vec == fold] == 1, -log(x), -log(1 - x)))
    })

    # format output
    out <- list()
    out$GQW <- pred_matrix1
    out$risk <- NULL
    if(!is.null(fold)){
        out$risk <- risk
    }
    return(out)
}
