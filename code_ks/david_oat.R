#' OAT HAL
#'
#' Fit a custom version of HAL for outcome adaptive TMLE estimation
#'
#' @param W A data.frame of predictors
#' @param A A numeric binary (0/1) treatment vector
#' @param Y A numeric outcome vector
#' @param V The number of cross-validation folds (must be at least 3 or \code{glmnet} will complain)
#' @param outcome_family "gaussian" or "binomial" (implies loss function to use)
#' @param lambda_seq A sequence of lambda values to be passed to \code{glmnet}.
#' In general, the \code{glmnet} defaults do not work well for the purposes of HAL.
#' The default is something that has yielded good results in the past, but needs
#' to be validated for each simulation.
#'
#' @return A data.frame. \code{QAW} = the HAL fit at CV-selected lambda for Qbar evaluated at (A_i, W_i),
#' for i = 1,...,n ; \code{Q1W} = the HAL fit at CV-selected lambda for Qbar evaluated at (1, W_i) ;
#' \code{Q0W} = the HAL fit at CV-selected lambda for Qbar evaluated at (0, W_i) ;
#' \code{G1W} = the HAL fit at CV-selected lambda for G evaluated at W_i ;
#' \code{GQW} = the HAL fit of A ~ Q1W + Q0W evaluated at Q(1,W_i), Q(0,W_i) ;
#' \code{cv_} = the cross-validated HAL fit at CV-selected lambda of the above nuisance parameters
#' @importFrom hal9001 enumerate_basis make_design_matrix make_copy_map


oat_hal <- function(W, A, Y, V = 10, outcome_family = "gaussian",
                    lambda_seq = exp(seq(-1, -13, length=10000))){

    n <- length(A)
    # make a vector of cv-folds
    base_rep <- rep(round(n/V), V)
    mod_rep <- n%%V
    base_rep[seq_len(mod_rep)] <- base_rep[seq_len(mod_rep)] + 1
    fold_vec <- rep(seq_len(V), base_rep)
    fold_idx <- sapply(seq_len(V), function(x){ which(fold_vec == x) }, simplify = FALSE)

    # cross-validation routine
    cv_out <- sapply(seq_len(V), one_hal,
                     fold_vec = fold_vec,
                     W = W, A = A, Y = Y,
                     lambda_seq = lambda_seq,
                     outcome_family = outcome_family,
                     simplify = FALSE)
    # cv risks
    # outcome regression
    or_risks <- colMeans(Reduce("rbind",lapply(cv_out, "[[", "riskQ")))
    # take smallest lambda that has smallest risk
    lambda_or_idx <- which.min(or_risks)[1]
    # propensity score
    ps_risks <- colMeans(Reduce("rbind",lapply(cv_out, "[[", "riskG")))
    lambda_ps_idx <- which.min(ps_risks)[1]
    # re-fit on full data
    full_hal <- one_hal(fold = NULL, fold_vec = NULL,
                     W = W, A = A, Y = Y,
                     lambda_seq = lambda_seq,
                     outcome_family = outcome_family)
    # QAW at cv selected lambda
    QAW_cvselect <- full_hal$QAW[, lambda_or_idx]
    # Q1W at cv selected lambda
    Q1W_cvselect <- full_hal$Q1W[, lambda_or_idx]
    # Q0W at cv selected lambda
    Q0W_cvselect <- full_hal$Q0W[, lambda_or_idx]
    # cross-validated Q1W at cv selected lambda
    QAW_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
        x$QAW[, lambda_or_idx]
    }))
    # cross-validated Q1W at cv selected lambda
    Q1W_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
        x$Q1W[, lambda_or_idx]
    }))
    # cross-validated Q0W at cv selected lambda
    Q0W_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
        x$Q0W[, lambda_or_idx]
    }))

    # G1W at cv selected lambda
    G1W_cvselect <- full_hal$G1W[, lambda_ps_idx]
    # cross-validated G1W at cv selected lambda
    G1W_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
        x$G1W[, lambda_ps_idx]
    }))

    # HAL A ~ Q1W + Q0W
    # first check that some variation in outcome regressions
    sanity_check <- !(all(Q1W_cvselect == Q0W_cvselect) & (length(unique(Q1W_cvselect)) == 1))

    # make data matrix Q1W, Q0W
    if(sanity_check){
        X <- as.matrix(cbind(Q1W_cvselect, Q0W_cvselect))
        basis_list <- hal9001::enumerate_basis(X, degrees = NULL)
        x_basis <- hal9001:::make_design_matrix(X, basis_list)
        copy_map <- hal9001:::make_copy_map(x_basis)
        unique_columns <- as.numeric(names(copy_map))
        x_basis <- x_basis[, unique_columns]

        hal_lasso <- glmnet::cv.glmnet(x = x_basis, y = A, nfolds = V,
                family = "binomial", lambda = lambda_seq, foldid = fold_vec,
                keep = TRUE)
        lambda_goat_idx <- which(hal_lasso$lambda == hal_lasso$lambda.min)

        # get predictions in whole sample
        GQW_cvselect <- predict(hal_lasso, s = "lambda.min", type = "response", newx = x_basis)

        # cross-validated predictions
        GQW_cv_cvselect <- hal_lasso$fit.preval[,lambda_goat_idx]

    }else{
        GQW_cvselect <- GQW_cv_cvselect <- rep(mean(A), n)
    }

    # format outcome
    out <- data.frame(QAW = QAW_cvselect,
                      Q1W = Q1W_cvselect,
                      Q0W = Q0W_cvselect,
                      G1W = G1W_cvselect,
                      GQW = as.numeric(GQW_cvselect),
                      cv_QAW = QAW_cv_cvselect,
                      cv_Q1W = Q1W_cv_cvselect,
                      cv_Q0W = Q0W_cv_cvselect,
                      cv_G1W = G1W_cv_cvselect,
                      cv_GQW = as.numeric(GQW_cv_cvselect),
                      lambda_or = lambda_seq[lambda_or_idx],
                      lambda_ps = lambda_seq[lambda_ps_idx],
                      fold_vec = fold_vec)

    return(out)
}


#' A helper function for fitting a single fold of HAL. Can be used both for
#' fitting cross-validated HAL and fitting HAL to full data.
#'
#' @param fold The validation fold (a numeric 1:V)
#' @param fold_vec A vector of which fold each observation is in.
#' @param outcome_family "gaussian" or "binomial"
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param lambda_seq See documentation for oat_hal
#'
#' @return See function itself to understand what it's returning.
one_hal <- function(fold, fold_vec, outcome_family = "gaussian",
                    W, A, Y, lambda_seq = exp(seq(-1, -13, length = 10000))){
    n <- length(Y)
    n_valid <- sum(fold == fold_vec)
    n_train <- n - n_valid

    # if called during cross-validation
    if(!is.null(fold)){
        x_fit <- as.matrix(data.frame(A = A, W)[fold_vec != fold, , drop = FALSE])
        x2_fit <- as.matrix(W[fold_vec != fold, , drop = FALSE])
        y_fit <- Y[fold_vec != fold]
        y2_fit <- A[fold_vec != fold]
    # if called fitting to full data
    }else{
        x_fit <- as.matrix(data.frame(A = A, W))
        x2_fit <- as.matrix(W)
        y_fit <- Y
        y2_fit <- A
    }

    # for outcome regression
    basis_list <- hal9001::enumerate_basis(x_fit, NULL)
    x_basis <- hal9001:::make_design_matrix(x_fit, basis_list)
    copy_map <- hal9001:::make_copy_map(x_basis)
    unique_columns <- as.numeric(names(copy_map))
    # subset to non-duplicated columns
    x_basis_fit <- x_basis[, unique_columns]

    # fit outcome regression
    hal_lasso <- glmnet::glmnet(x = x_basis_fit, y = y_fit,
        family = outcome_family, lambda = lambda_seq,
        standardize = FALSE)
    # get outcome regression predictions
    beta_hat <- as.matrix(hal_lasso$beta)
    # intercept
    alpha_hat <- hal_lasso$a0
    rm(x_basis, x_basis_fit)
    gc()
    # for propensity
    basis_list2 <- hal9001::enumerate_basis(x2_fit, NULL)
    x_basis2 <- hal9001:::make_design_matrix(x2_fit, basis_list2)
    copy_map2 <- hal9001:::make_copy_map(x_basis2)
    unique_columns2 <- as.numeric(names(copy_map2))
    x_basis2_fit <- x_basis2[, unique_columns2]

    # fit propensity score
    hal_lasso2 <- glmnet::glmnet(x = x_basis2_fit, y = y2_fit,
        family = "binomial", lambda = lambda_seq, standardize = FALSE)
    # get propensity score predictions
    beta_hat2 <- as.matrix(hal_lasso2$beta)
    # intercept
    alpha_hat2 <- hal_lasso2$a0
    rm(x_basis2, x_basis_fit2)
    gc()

    # predictions on validation sample
    if(!is.null(fold)){
        # for selecting lambda
        new_x_fit1 <- as.matrix(data.frame(A = A, W)[fold_vec == fold, , drop = FALSE])
        # for use later in TMLE
        new_x_fit2 <- as.matrix(rbind(data.frame(A = 1, W)[fold_vec == fold, , drop = FALSE],
                            data.frame(A = 0, W)[fold_vec == fold, , drop = FALSE]))
        # for propensity score
        new_x_fit3 <- as.matrix(W[fold_vec == fold, , drop = FALSE])
    }else{
        # for selecting lambda
        new_x_fit1 <- as.matrix(data.frame(A = A, W))
        # for use later in TMLE
        new_x_fit2 <- as.matrix(rbind(data.frame(A = 1, W),
                            data.frame(A = 0, W)))
        new_x_fit3 <- as.matrix(W)
    }
    # make HAL design for getting predictions to select lambda
    new_x_basis1 <- hal9001:::make_design_matrix(new_x_fit1, basis_list)
    rm(new_x_fit1)
    gc()
    new_x_basis1 <- as.matrix(new_x_basis1[, unique_columns])
    # make HAL design for getting predictions used by TMLE later
    new_x_basis2 <- hal9001:::make_design_matrix(new_x_fit2, basis_list)
    rm(new_x_fit2)
    gc()
    new_x_basis2 <- as.matrix(new_x_basis2[, unique_columns])
    # for propensity score
    new_x_basis3 <- hal9001:::make_design_matrix(new_x_fit3, basis_list2)
    rm(new_x_fit3)
    gc()
    new_x_basis3 <- as.matrix(new_x_basis3[, unique_columns2])


    # prediction matrices
    # for (A,W)
    or_pred_matrix1 <- cbind(rep(1, ifelse(is.null(fold), n, n_valid)), new_x_basis1) %*% rbind(alpha_hat, beta_hat)
    rm(new_x_basis1)
    gc()
    # for rbind((1,W),(0,W))
    or_pred_matrix2 <- cbind(rep(1, ifelse(is.null(fold), 2* n, n_valid*2)), new_x_basis2) %*% rbind(alpha_hat, beta_hat)
    rm(new_x_basis2)
    gc()
    if(outcome_family == "binomial"){
        or_pred_matrix1 <- apply(or_pred_matrix1, 2, plogis)
        or_pred_matrix2 <- apply(or_pred_matrix2, 2, plogis)
    }

    # predictions
    logit_ps_pred_matrix <- cbind(rep(1, ifelse(is.null(fold), n, n_valid)), new_x_basis3) %*% rbind(alpha_hat2, beta_hat2)
    ps_pred_matrix <- apply(logit_ps_pred_matrix, 2, plogis)
    rm(new_x_basis3)
    gc()
    # compute MSE/loglik
    if(outcome_family == "gaussian"){
        or_risk <- apply(or_pred_matrix1, 2, function(x){
            mean((Y[fold_vec == fold] - x)^2)
        })
    }else{
        or_risk <- apply(or_pred_matrix1, 2, function(x){
            mean(ifelse(Y[fold_vec == fold] == 1, -log(x), -log(1 - x)))
        })
    }
    ps_risk <- apply(ps_pred_matrix, 2, function(x){
        mean(ifelse(A[fold_vec == fold] == 1, -log(x), -log(1 - x)))
    })

    # format output
    out <- list()
    out$QAW <- or_pred_matrix1
    out$riskQ <- NULL
    out$riskG <- NULL
    if(!is.null(fold)){
        out$riskQ <- or_risk
        out$riskG <- ps_risk
    }
    if(is.null(fold)){
        Q1W_idx <- 1:n
        Q0W_idx <- (n+1):(2*n)
    }else{
        Q1W_idx <- 1:n_valid
        Q0W_idx <- (n_valid + 1):(2*n_valid)
    }
    out$Q1W <- or_pred_matrix2[Q1W_idx, , drop = FALSE]
    out$Q0W <- or_pred_matrix2[Q0W_idx, , drop = FALSE]
    out$G1W <- ps_pred_matrix
    return(out)
}

# some sanity testing
# n <- 500
# W <- data.frame(W1 = runif(n, -1, 1))
# A <- rbinom(n, 1, plogis(W$W1^2))
# Y <- rbinom(n, 1, plogis(-W$W1^2))
# G1W <- function(W){
#   plogis(W$W1)
# }
# QAW <- function(A, W){
#   plogis(W$W1)
# }

# debug(oat_hal)
# fit <- oat_hal(W = W, A = A, Y = Y, V = 3, outcome_family = "binomial")

#' Function to generate modified Kang and Schafer data
#' @param n Sample size
#' @return A list with \code{W} a data.frame of covariates, \code{A} a binary treatment
#' and \code{Y} a continuous outcome
mod_ks <- function(n){
    Z <- data.frame(matrix(runif(n * 4, min = -2, max = 2), ncol = 4))
    colnames(Z) <- paste0("Z", 1:4)
    # change z1
    Z$Z1 <- runif(n, min = 0.5, max = 2)
    # add instrument
    Z$Z5 <- runif(n, -2, 2)
    # logit of ps
    # logit_G0 <- (- Z$Z1 + 0.5 * Z$Z2 - 0.25 * Z$Z3 - 0.1 * Z$Z4)/2 + Z$Z5
    logit_G0 <- (- Z$Z1 + 0.5 * Z$Z2 - Z$Z3 - 0.1 * Z$Z4) + Z$Z5 + 0.75*Z$Z5^2
    G0 <- plogis(logit_G0)

    # draw A
    A <- rbinom(n, 1, G0)

    # draw Y
    Q0 <- 210 + 27.4 * Z$Z1 + 13.7 * Z$Z2 + 13.7 * Z$Z3 + 13.7 * Z$Z4
    Y <- Q0 + rnorm(n)

    # make W from Z
    W <- make_W_from_Z(Z)

    # output
    out <- list(W = W, A = A, Y = Y, Q1 = Q1W, Q0 = Q0W, pw = G)
}

#' Helper function to go from Z to W in Kang and Schafer set up
#' @param W data.frame with W covariates
#' @return A data.frame with covariates on Z scale
make_Z_from_W <- function(W){
    Z1 <- 2*log(W$W1)
    Z2 <- (W$W2 - 10) * (1 + exp(Z1))
    Z3 <- ((W$W3^(1/3) - 0.6) * 25) / Z1
    Z4 <- W$W4^(1/2) - 20 - Z2
    return(data.frame(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = W$W5))
}

#' Helper function to go from W to Z in Kang and Schafer set up
#' @param Z data.frame with Z covariates
#' @return A data.frame with covariates on W scale
make_W_from_Z <- function(Z){
    W1 <- exp(Z$Z1/2)
    W2 <- Z$Z2/(1 + exp(Z$Z1)) + 10
    W3 <- ((Z$Z1 * Z$Z3) / 25 + 0.6)^3
    W4 <- (Z$Z2 + Z$Z4 + 20)^2
    return(data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4, W5 = Z$Z5))
}

#' Outcome regression for modified Kang and Schafer
#' @param W data.frame of covariates on W scale
#' @return A vector of Q(1, W)
Q1W <- function(W){
    Z <- make_Z_from_W(W)
    return(210 + 27.4 * Z$Z1 + 13.7 * Z$Z2 + 13.7 * Z$Z3 + 13.7 * Z$Z4)
}

#' Outcome regression for modified Kang and Schafer
#' Note because no treatment effect, is the same as Q1W
#' @param W data.frame of covariates on W scale
#' @return A vector of Q(0, W)
Q0W <- function(W){
    Z <- make_Z_from_W(W)
    return(210 + 27.4 * Z$Z1 + 13.7 * Z$Z2 + 13.7 * Z$Z3 + 13.7 * Z$Z4)
}

#' Propensity score for modified Kang and Schafer
#' @param W data.frame of covariates on W scale
#' @return A vector of P(A = 1 | W)
G <- function(W){
    Z <- make_Z_from_W(W)
    logit_G0 <- (- Z$Z1 + 0.5 * Z$Z2 - 0.25 * Z$Z3 - 0.1 * Z$Z4)/2 + Z$Z5
    G0 <- plogis(logit_G0)
    return(G0)
}
