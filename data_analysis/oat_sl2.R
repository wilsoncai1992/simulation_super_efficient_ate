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
    "SL.glm", "SL.earth", "SL.ranger", "SL.glmnet", "SL.gbm", "SL.mean", "SL.step.forward"
  ),
  lambda_seq = exp(seq(-1, -13, length = 200)),
  save_files = TRUE,
  save_dir = "~/Documents/"
) {
    n <- length(Y)
    # make a vector of cv-folds
    base_rep <- rep(round(n/V), V)
    mod_rep <- n%%V
    base_rep[seq_len(mod_rep)] <- base_rep[seq_len(mod_rep)] + 1
    fold_vec <- rep(seq_len(V), base_rep)
    fold_idx <- sapply(seq_len(V), function(x){ which(fold_vec == x) }, simplify = FALSE)

    # fit a super learner for OR
    newX_or <- rbind(
      data.frame(W, A2 = as.numeric(A == 2), A1 = as.numeric(A == 1)),
      data.frame(W, A2 = 1, A1 = 0),
      data.frame(W, A2 = 0, A1 = 1),
      data.frame(W, A2 = 0, A1 = 0)
    )

    cat("Fitting outcome regression SL \n")
    sl_or <- SuperLearner(Y = Y, X = data.frame(W, A2 = as.numeric(A == 2), A1 = as.numeric(A == 1)),
                          newX = newX_or,
                          SL.library = SL.library, family = outcome_family,
                          cvControl = list(V = V, validRows = fold_idx),
                          method = ifelse(outcome_family == "gaussian",
                                          "method.CC_LS2",
                                          "method.CC_nloglik2"))
    if(save_files){
      save(sl_or, file = paste0(save_dir, "sl_or.RData"))
    }
    # extract predictions
    idx_aw <- 1:n
    idx_2w <- (n+1):(2*n)
    idx_1w <- (2*n+1):(3*n)
    idx_0w <- (3*n+1):(4*n)
    QAW <- sl_or$SL.predict[idx_aw]
    Q2W <- sl_or$SL.predict[idx_2w]
    Q1W <- sl_or$SL.predict[idx_1w]
    Q0W <- sl_or$SL.predict[idx_0w]

    cat("Fitting cv outcome regression SL \n")
    # cross-validated version
    cv_sl_or <- CV.SuperLearner(Y = Y, X = data.frame(W, A2 = as.numeric(A==2), A1 = as.numeric(A ==1)),
                          SL.library = SL.library, family = outcome_family,
                          innerCvControl = list(list(V = V),list(V = V),list(V = V),
                                           list(V = V),list(V = V)),
                          cvControl = list(V = V, validRows = fold_idx),
                          control = list(saveFitLibrary = TRUE),
                          method = ifelse(outcome_family == "gaussian",
                                          "method.CC_LS2",
                                          "method.CC_nloglik2"))
    if(save_files){
      save(cv_sl_or, file = paste0(save_dir, "cv_sl_or.RData"))
    }
    cv_QAW <- as.numeric(cv_sl_or$SL.predict)
    cv_Q1W <- rep(NA, n)
    cv_Q2W <- rep(NA, n)
    cv_Q0W <- rep(NA, n)
    for(v in 1:5){
        this_idx <- cv_sl_or$folds[[v]]
        cv_Q1W[this_idx] <-
            predict(cv_sl_or$AllSL[[v]], newdata = data.frame(W, A2 = 0, A1 = 1)[this_idx,])$pred
        cv_Q2W[this_idx] <-
            predict(cv_sl_or$AllSL[[v]], newdata = data.frame(W, A2 = 1, A1 = 0)[this_idx,])$pred
        cv_Q0W[this_idx] <-
            predict(cv_sl_or$AllSL[[v]], newdata = data.frame(W, A2 = 0, A1 = 0)[this_idx,])$pred
    }
    cat("Fitting first propensity score SL \n")
    # fit a super learner for PS A = 2
    sl_ps <- SuperLearner(Y = as.numeric(A == 2),
                          X = W,
                          newX = W,
                          SL.library = SL.library, family = binomial(),
                          cvControl = list(V = V, validRows = fold_idx),
                          method = "method.CC_nloglik2")
    if(save_files){
      save(sl_ps, file = paste0(save_dir, "sl_ps.RData"))
    }
    # extract predictions
    G2W <- as.numeric(sl_ps$SL.predict)

    # fit a super learner for PS A = 1 | A != 2
    # first remove A == 2 folks from validation folds
    A2_idx <- which(A == 2)
    not_A2_idx <- which(A != 2)
    n_not_A2 <- length(not_A2_idx)
    id2 <- rep(NA, n)
    id2[A != 2] <- seq_len(n_not_A2)
    
    new_fold_idx <- lapply(fold_idx, function(x){
      tmp <- id2[x]
      tmp[!is.na(tmp)]
    })
    
    cat("Fitting second propensity score SL \n")
    sl_ps2 <- SuperLearner(Y = as.numeric(A == 1)[ A != 2 ],
                          X = W[A !=2,],
                          newX = W,
                          SL.library = SL.library, family = binomial(),
                          cvControl = list(V = V, validRows = new_fold_idx),
                          method = "method.CC_nloglik2")
    if(save_files){
      save(sl_ps2, file = paste0(save_dir, "sl_ps2.RData"))
    }
    # extract predictions
    G1W <- as.numeric(sl_ps2$SL.predict) * (1 - G2W)

    G0W <- 1 - G1W - G2W

    cat("Fitting first cv propensity score SL \n")
    # cross-validated version
    cv_sl_ps <- CV.SuperLearner(Y = as.numeric(A == 2), X = W,
                          SL.library = SL.library, family = binomial(),
                          innerCvControl = list(list(V = V),list(V = V),list(V = V),
                                           list(V = V),list(V = V)),
                          cvControl = list(V = V, validRows = fold_idx),
                          control = list(saveFitLibrary = TRUE),
                          method = "method.CC_nloglik2")
    if(save_files){
      save(cv_sl_ps, file = paste0(save_dir, "cv_sl_ps.RData"))
    }
    cv_G2W <- as.numeric(cv_sl_ps$SL.predict)

    cat("Fitting second cv propensity score SL \n")
    cv_sl_ps2 <- CV.SuperLearner(Y = as.numeric(A==1)[ A != 2 ], X = W[ A != 2 ,],
                          SL.library = SL.library, family = binomial(),
                          innerCvControl = list(list(V = V),list(V = V),list(V = V),
                                           list(V = V),list(V = V)),
                          cvControl = list(V = V, validRows = new_fold_idx),
                          control = list(saveFitLibrary = TRUE),
                          method = "method.CC_nloglik2")
    if(save_files){
      save(cv_sl_ps2, file = paste0(save_dir, "cv_sl_ps2.RData"))
    }
    cv_G1W <- rep(NA, n)
    cv_G0W <- rep(NA, n)
    for(v in 1:5){
        this_idx <- cv_sl_ps$folds[[v]]
        cv_G1W[this_idx] <-
            predict(cv_sl_ps2$AllSL[[v]], newdata = data.frame(W)[this_idx,])$pred * (1 - cv_G2W[this_idx])
        cv_G0W[this_idx] <-
            1 - cv_G1W[this_idx] - cv_G2W[this_idx]
    }

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
          cat("Fitting HALs \n")
            # cross-validation routine based on CV-Q
          QaW_list <- list(Q0W = matrix(Q0W,ncol = 1), Q1W = matrix(Q1W,ncol = 1), Q2W = matrix(Q2W,ncol = 1))
          GQW_list <- vector(mode = 'list', length = 3)
          cv_GQW_list <- vector(mode = 'list', length = 3)
          for(a in 0:2){
            cv_out <- sapply(seq_len(V), one_red_hal,
                             fold_vec = fold_vec,
                             QaW = QaW_list[[a+1]],
                             A = as.numeric(A == a),
                             lambda_seq = lambda_seq,
                             simplify = FALSE)
            # outcome regression
            risks <- colMeans(Reduce("rbind",lapply(cv_out, "[[", "risk")))
            # take smallest lambda that has smallest risk
            lambda_idx <- which.min(risks)[1]

            # re-fit on full data
            full_hal <- one_red_hal(fold = NULL, fold_vec = NULL,
                         A = as.numeric(A == a), QaW = QaW_list[[a+1]],
                         lambda_seq = lambda_seq)

            # GQW at cv selected lambda
            GQW_list[[a+1]] <- full_hal$GQW[, lambda_idx]
            # cross-validated Q1W at cv selected lambda
            cv_GQW_list[[a+1]] <- Reduce("c", lapply(cv_out, function(x){
                x$GQW[, lambda_idx]
            }))
          }
        }
    }else{
        GQW_cvselect <- GQW_cv_cvselect <- rep(mean(A), n)
    }

    # format outcome
    out <- data.frame(QAW = QAW,
                      Q2W = Q2W,
                      Q1W = Q1W,
                      Q0W = Q0W,
                      G0W = G0W,
                      G1W = G1W,
                      G2W = G2W,
                      GQ0W = as.numeric(GQW_list[[1]]),
                      GQ1W = as.numeric(GQW_list[[2]]),
                      GQ2W = as.numeric(GQW_list[[3]]),
                      cv_QAW = cv_QAW,
                      cv_Q2W = cv_Q2W,
                      cv_Q1W = cv_Q1W,
                      cv_Q0W = cv_Q0W,
                      cv_G2W = cv_G2W,
                      cv_G1W = cv_G1W,
                      cv_G0W = cv_G0W,
                      cv_GQ0W = as.numeric(cv_GQW_list[[1]]),
                      cv_GQ1W = as.numeric(cv_GQW_list[[2]]),
                      cv_GQ2W = as.numeric(cv_GQW_list[[3]]),
                      fold_vec = fold_vec)
    return(out)
}


one_red_hal <- function(fold, fold_vec, outcome_family = "binomial",
                        A, QaW, lambda_seq = exp(seq(-1,-13,length=10000))){
    n <- length(A)
    n_valid <- sum(fold == fold_vec)
    n_train <- n - n_valid

    if(!is.null(fold)){
        x_fit <- matrix(QaW[fold_vec != fold, ], ncol = 1)
        y_fit <- A[fold_vec != fold]
    # if called fitting to full data
    }else{
        x_fit <- matrix(QaW, ncol = 1)
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
        new_x_fit1 <- cbind(QaW)[fold_vec == fold, , drop = FALSE]
    }else{
        new_x_fit1 <- cbind(QaW)
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

# less sensitive version of cc method in SL
method.CC_LS2 <- function () 
{
    computeCoef = function(Z, Y, libraryNames, verbose, obsWeights, 
        errorsInLibrary = NULL, ...) {
        cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x - 
            Y)^2))
        names(cvRisk) <- libraryNames
        compute <- function(x, y, wt = rep(1, length(y))) {
            wX <- sqrt(wt) * x
            wY <- sqrt(wt) * y
            D <- crossprod(wX)
            d <- crossprod(wX, wY)
            A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
            bvec <- c(1, rep(0, ncol(wX)))
            fit <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = A, 
                bvec = bvec, meq = 1)
            invisible(fit)
        }
        modZ <- Z
        naCols <- which(apply(Z, 2, function(z) {
            all(z == 0)
        }))
        anyNACols <- length(naCols) > 0
        if (anyNACols) {
            warning(paste0(paste0(libraryNames[naCols], collapse = ", "), 
                " have NAs.", "Removing from super learner."))
        }
        tol <- 3
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "), 
                " are duplicates of previous learners.", " Removing from super learner."))
        }
        if (anyDupCols | anyNACols) {
            rmCols <- unique(c(naCols, dupCols))
            modZ <- Z[, -rmCols]
        }
        fit <- compute(x = modZ, y = Y, wt = obsWeights)
        coef <- fit$solution
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] = 0
        }
        if (anyDupCols | anyNACols) {
            ind <- c(seq_along(coef), rmCols - 0.5)
            coef <- c(coef, rep(0, length(rmCols)))
            coef <- coef[order(ind)]
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        if (!sum(coef) > 0) 
            warning("All algorithms have zero weight", call. = FALSE)
        list(cvRisk = cvRisk, coef = coef, optimizer = fit)
    }
    computePred = function(predY, coef, ...) {
        predY %*% matrix(coef)
    }
    out <- list(require = "quadprog", computeCoef = computeCoef, 
        computePred = computePred)
    invisible(out)
}

method.CC_nloglik2 <- function () 
{
    computePred = function(predY, coef, control, ...) {
        if (sum(coef != 0) == 0) {
            stop("All metalearner coefficients are zero, cannot compute prediction.")
        }
        plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*% 
            matrix(coef[coef != 0]))
    }
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, 
        verbose, ...) {
        tol <- 3
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        modZ <- Z
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "), 
                " are duplicates of previous learners.", " Removing from super learner."))
            modZ <- modZ[, -dupCols]
        }
        modlogitZ <- trimLogit(modZ, control$trimLogit)
        logitZ <- trimLogit(Z, control$trimLogit)
        cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights * 
            ifelse(Y, plogis(x, log.p = TRUE), plogis(x, log.p = TRUE, 
                lower.tail = FALSE))))
        names(cvRisk) <- libraryNames
        obj_and_grad <- function(y, x, w = NULL) {
            y <- y
            x <- x
            function(beta) {
                xB <- x %*% cbind(beta)
                loglik <- y * plogis(xB, log.p = TRUE) + (1 - 
                  y) * plogis(xB, log.p = TRUE, lower.tail = FALSE)
                if (!is.null(w)) 
                  loglik <- loglik * w
                obj <- -2 * sum(loglik)
                p <- plogis(xB)
                grad <- if (is.null(w)) 
                  2 * crossprod(x, cbind(p - y))
                else 2 * crossprod(x, w * cbind(p - y))
                list(objective = obj, gradient = grad)
            }
        }
        lower_bounds = rep(0, ncol(modZ))
        upper_bounds = rep(1, ncol(modZ))
        if (anyNA(cvRisk)) {
            upper_bounds[is.na(cvRisk)] = 0
        }
        r <- nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)), 
            eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds, 
            ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) - 
                1), eval_jac_g_eq = function(beta) rep(1, length(beta)), 
            opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
        if (r$status < 1 || r$status > 4) {
            warning(r$message)
        }
        coef <- r$solution
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] <- 0
        }
        if (anyDupCols) {
            ind <- c(seq_along(coef), dupCols - 0.5)
            coef <- c(coef, rep(0, length(dupCols)))
            coef <- coef[order(ind)]
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
        return(out)
    }
    list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}
