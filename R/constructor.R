#' Construct ARGO object
#'
#' Wrapper for ARGO. The real work horse is glmnet package and/or linear model.
#'
#' This function takes the time series and exogenous variables (optional) as
#' input, and produces out-of-sample prediction for each time point.
#'
#' @param data response variable as xts, last element can be NA. If the response
#' is later revised, it should be an xts that resembles upper triangular square
#' matrix, with each column being the data available as of date of column name
#' @param exogen exogenous predictors, default is NULL
#' @param N_lag vector of the AR model lags used,
#' if NULL then no AR lags will be used
#'
#' @param use_all_previous boolean variable indicating whether to use "all available data"
#' (when \code{TRUE}) or "a sliding window" (when \code{FALSE}) for training
#' @param N_training number of training points, if \code{use_all_previous} is true,
#' this is the least number of training points required
#' @param alpha penalty between lasso and ridge, alpha=1 represents lasso,
#' alpha=0 represents ridge, alpha=NA represents no penalty
#' @param mc.cores number of cores to compute argo in parallel
#' @param schedule list to specify prediction schedule. Default to have \code{y_gap} as 1, and \code{forecast} as 0, 
#' i.e., nowcasting with past week ILI available from CDC.
#'
#' @return A list of following named objects
#' \itemize{
#' \item \code{pred} An xts object with the same index as input,
#' which contains historical nowcast estimation
#'
#' \item \code{coef} A matrix contains historical coefficient values of the predictors.
#'
#' \item \code{parm} Parameter values passed to argo function.
#'
#' \item \code{penalfac} the value of lambda ratio selected by cross-validation,
#' NULL if \code{lamid} is NULL or has only one level.
#'
#' \item \code{penalregion} the lambda ratios that has a cross validation error
#' within one standard error of minimum cross validation error
#' }
#'
#' @references
#' Yang, S., Santillana, M., & Kou, S. C. (2015). Accurate estimation of influenza epidemics using Google search data via ARGO. Proceedings of the National Academy of Sciences. <doi:10.1073/pnas.1515373112>.
#' @examples
#' GFT_xts <- xts::xts(exp(matrix(rnorm(180), ncol=1)), order.by = Sys.Date() - (180:1))
#' randomx <- xts::xts(exp(matrix(rnorm(180*100), ncol=100)), order.by = Sys.Date() - (180:1))
#' \donttest{
#' argo_result1 <- argo(GFT_xts)
#' argo_result2 <- argo(GFT_xts, exogen = randomx)
#' }
#' @export
argo <- function(data, exogen=xts::xts(NULL), N_lag=1:52, N_training=104,
                 alpha=1, use_all_previous=FALSE, mc.cores=1, schedule = list()){

  if(is.null(schedule$y_gap)){
    schedule$y_gap <- 1 # default information gap is 1
  }
  if(is.null(schedule$forecast)){
    schedule$forecast <- 0 # default is now-cast
  }

  parm <- list(N_lag = N_lag, N_training = N_training,
               alpha = alpha, use_all_previous = use_all_previous,
               schedule = schedule)
  if(ncol(data)==1){
    data_mat <- matrix(rep(data, nrow(data)), nrow=nrow(data))
    colnames(data_mat) <- as.character(index(data))
    for(i in schedule$y_gap:min(ncol(data_mat), ncol(data_mat)-1+schedule$y_gap)){
      data_mat[(i-schedule$y_gap+1):nrow(data_mat),i] <- NA
    }
    data_mat <- xts::xts(data_mat, zoo::index(data))
    data <- data_mat
  }

  lasso.pred <- c()
  lasso.coef <- list()

  if(length(exogen)>0) # exogenous variables must have the same timestamp as y
    if(!all(zoo::index(data)==zoo::index(exogen)))
      stop("error in data and exogen: their time steps must match")

  starttime <- N_training+max(c(N_lag,0)) + 2*schedule$y_gap + schedule$forecast
  endtime <- nrow(data)


  each_iteration <- function(i) {
    if(use_all_previous){
      training_idx <- (schedule$y_gap + schedule$forecast + max(c(N_lag, 0))):(i - schedule$y_gap)
    }else{
      training_idx <- (i - N_training + 1):i - schedule$y_gap
    }

    lagged_y <- sapply(N_lag, function(l)
      as.numeric(diag(data.matrix(data[training_idx - l + 1 - schedule$y_gap - schedule$forecast,training_idx- schedule$forecast]))))

    if(length(lagged_y) == 0){
      lagged_y <- NULL
    }else{
      colnames(lagged_y) <- paste0("lag_", N_lag+schedule$y_gap-1)
    }

    if(length(exogen) > 0 && schedule$y_gap > 0){
      xmat <- lapply(1:schedule$y_gap, function(l)
        as.matrix(exogen[training_idx - schedule$forecast - l + 1, ]))
      xmat <- do.call(cbind, xmat)
      design_matrix <- cbind(lagged_y, xmat)
    }else{
      design_matrix <- cbind(lagged_y)
    }

    y.response <- data[training_idx, i]

    if(is.finite(alpha)){
      lasso.fit <-
        glmnet::cv.glmnet(x=design_matrix,y=y.response,nfolds=10,
                          grouped=FALSE,alpha=alpha)

      lam.s <- lasso.fit$lambda.1se

    }else{
      lasso.fit <- lm(y.response ~ ., data=data.frame(design_matrix))
    }

    if(is.finite(alpha)){
      lasso.coef[[i]] <- as.matrix(coef(lasso.fit, lambda = lam.s))
    }else{
      lasso.coef[[i]] <- as.matrix(coef(lasso.fit))
    }
    lagged_y_next <- matrix(sapply(N_lag, function(l)
      as.numeric(data[i - schedule$y_gap + 1 - l, i])), nrow=1)

    if(length(lagged_y_next) == 0)
      lagged_y_next <- NULL
    if(length(exogen) > 0 && schedule$y_gap > 0){
      xmat.new <- lapply(1:(schedule$y_gap), function(l)
        data.matrix(exogen[i - l + 1, ]))
      xmat.new <- do.call(cbind, xmat.new)
      newx <- cbind(lagged_y_next, xmat.new)
    }else{
      newx <- lagged_y_next
    }
    if(is.finite(alpha)){
      lasso.pred[i] <- predict(lasso.fit, newx = newx, s = lam.s)
    }else{
      colnames(newx) <- names(lasso.fit$coefficients)[-1]
      newx <- as.data.frame(newx)
      lasso.pred[i] <- predict(lasso.fit, newdata = newx)
    }
    result_i <- list()
    result_i$pred <- lasso.pred[i]
    result_i$coef <- lasso.coef[[i]]
    result_i
  }

  result_all <- parallel::mclapply(starttime:endtime, each_iteration,
                                   mc.cores = mc.cores, mc.set.seed = FALSE)

  lasso.pred[starttime:endtime] <- sapply(result_all, function(x) x$pred)
  lasso.coef <- lapply(result_all, function(x) x$coef)

  data$predict <- lasso.pred
  lasso.coef <- do.call("cbind", lasso.coef)
  colnames(lasso.coef) <- as.character(zoo::index(data))[starttime:endtime]
  argo <- list(pred = data$predict, coef = lasso.coef, parm = parm)
  class(argo) <- "argo"
  argo
}


#' ARGO second step
#'
#' Wrapper for ARGO second step. Best linear predictor / Bayesian posterior
#'
#'
#' @param truth prediction target
#' @param argo1.p argo first step prediction
#' @param argo.nat.p argo national level prediction
#'
#' @importFrom Matrix bdiag
#' @import stats
#' @references
#' Shaoyang Ning, Shihao Yang, S. C. Kou. Accurate Regional Influenza Epidemics Tracking Using Internet Search Data. Scientific Reports
#' @examples
#' truth <- xts::xts(exp(matrix(rnorm(180*10), ncol=10)), order.by = Sys.Date() - (180:1))
#' argo1.p <- xts::xts(exp(matrix(rnorm(180*10), ncol=10)), order.by = Sys.Date() - (180:1))
#' argo.nat.p <- xts::xts(exp(matrix(rnorm(180*10), ncol=10)), order.by = Sys.Date() - (180:1))
#' argo2result <- argo2(truth, argo1.p, argo.nat.p)
#' @export
argo2 <- function(truth, argo1.p, argo.nat.p){
  naive.p <- lag(truth, 1)

  common_idx <- zoo::index(na.omit(merge(truth, naive.p, argo1.p, argo.nat.p)))
  if(ncol(argo.nat.p)==1){
    argo.nat.p <- do.call(cbind, lapply(1:10, function(i) argo.nat.p))
  }
  colnames(argo.nat.p) <- colnames(truth)

  Z <- truth - naive.p
  W <- merge(lag(Z, 1), argo1.p - naive.p, argo.nat.p - naive.p)

  arg2zw_result <- argo2zw(as.matrix(na.omit(merge(Z, W))))
  Z.hat <- arg2zw_result$Z.hat
  Z.hat <- xts(Z.hat, as.Date(rownames(Z.hat)))

  argo2.p <- Z.hat + naive.p

  heat.vec <- na.omit(merge(Z, lag(Z), argo1.p - truth, argo.nat.p - truth))
  colnames(heat.vec) <- paste0(rep(c("detla.ili.", "err.argo.", "err.nat.", "lag.detla.ili."), each=10),
                               rep(1:10,3))
  result_wrapper <- list(onestep=argo1.p, twostep=argo2.p, naive=naive.p, truth=truth,
                         heat.vec=heat.vec)
  c(result_wrapper, arg2zw_result)
}

argo2zw <- function(zw.mat){
  if(is.null(rownames(zw.mat)))
    stop("row name must exist as index")
  projection.mat <- list()
  mean.mat <- list()
  zw_used <- list()

  sigma_ww.structured <- sigma_ww.empirical <-
    sigma_zw.structured <- sigma_zw.empirical <-
    heat.vec.structured <-
    sigma_zwzw.structured <- sigma_zwzw.empirical <- list()

  epsilon <- list()
  Z <- zw.mat[,1:10]
  W <- zw.mat[,-(1:10)]
  W1 <- W[,1:10]
  W2 <- W[,11:20]
  W3 <- W[,21:30]

  Z.hat <- Z
  Z.hat[] <- NA
  for(it in 105:nrow(zw.mat)){
    training_idx <- (it-104):(it-1)
    t.now <- rownames(zw.mat)[it]
    if(is.null(t.now))
      t.now <- it

    epsilon[[as.character(t.now)]] <- list()

    sigma_zz <- var(Z[training_idx,])
    sigma_zz_chol <- chol(sigma_zz)
    epsilon[[as.character(t.now)]]$z <- solve(t(sigma_zz_chol), t(data.matrix(Z[training_idx,])))

    zw_used[[as.character(t.now)]] <- cbind(Z[training_idx,], W[training_idx,])

    m1 <- cor(Z[training_idx,], W[training_idx,1:10])
    m2 <- cor(Z[training_idx,])
    rho <- sum(m1*m2)/sum(m2^2)

    d.gt <- diag(diag(var(W2[training_idx,] - Z[training_idx,])))
    epsilon[[as.character(t.now)]]$argoreg <- t(data.matrix((W2 - Z)[training_idx,]))
    epsilon[[as.character(t.now)]]$argoreg <- epsilon[[as.character(t.now)]]$argoreg / sqrt(diag(d.gt))
    sigma.nat <- var((W3 - Z)[training_idx,])
    epsilon[[as.character(t.now)]]$argonat <- t(data.matrix((W3 - Z)[training_idx,]))
    sigma.nat_chol <- chol(sigma.nat)
    epsilon[[as.character(t.now)]]$argonat <- solve(t(sigma.nat_chol), epsilon[[as.character(t.now)]]$argonat)

    sigma_ww <- rbind(
      cbind(sigma_zz, rho*sigma_zz, rho*sigma_zz),
      cbind(rho*sigma_zz, sigma_zz + d.gt, sigma_zz),
      cbind(rho*sigma_zz, sigma_zz, sigma_zz + sigma.nat)
    )

    sigma_zw <- cbind(rho*sigma_zz, sigma_zz, sigma_zz)

    mu_w <- colMeans(W[training_idx,])
    mu_z <- colMeans(Z[training_idx,])

    d_ww <- diag(diag(var(W[training_idx,])))

    # not shrinked
    pred.blp <- mu_z + sigma_zw %*% solve(sigma_ww, W[t.now,] - mu_w)

    Kzz <- solve((1-rho^2)*sigma_zz)
    Kgt <- diag(1/diag(var((W2 - Z)[training_idx,])))
    Knat <- solve(var((W3 - Z)[training_idx,]))

    pred.bayes <- mu_z +
      solve(Kzz+Knat+Kgt,
            Knat%*%(W[t.now,21:30] - mu_w[21:30]) + Kgt%*%(W[t.now,11:20] - mu_w[11:20]) +
              Kzz%*%(rho*(W[t.now,1:10] - mu_w[1:10])))

    if(all(is.finite(pred.blp))){
      stopifnot(all(abs(pred.blp-pred.bayes) < 1e-8))
    }

    # shrinked
    z.hat <- mu_z + sigma_zw %*% solve(sigma_ww + d_ww, (W[t.now,] - mu_w))
    Z.hat[t.now, ] <- t(z.hat)

    projection.mat[[as.character(t.now)]] <- sigma_zw %*% solve(sigma_ww + d_ww)
    mean.mat[[as.character(t.now)]] <- c(mu_z, mu_w)
    sigma_ww.structured[[as.character(t.now)]] <- sigma_ww
    sigma_ww.empirical[[as.character(t.now)]] <- var(W[training_idx,])
    sigma_zw.structured[[as.character(t.now)]] <- sigma_zw
    sigma_zw.empirical[[as.character(t.now)]] <- cov(Z[training_idx,], W[training_idx,])

    sigma_zwzw.structured[[as.character(t.now)]] <- rbind(
      cbind(sigma_zz, sigma_zw),
      cbind(t(sigma_zw), sigma_ww)
    )
    sigma_zwzw.empirical[[as.character(t.now)]] <- var(cbind(Z, W)[training_idx,])

    heat.vec.struc <- rbind(
      cbind(sigma_zz, rho*sigma_zz),
      cbind(rho*sigma_zz, sigma_zz)
    )
    heat.vec.struc <- Matrix::bdiag(heat.vec.struc, d.gt, sigma.nat)
    heat.vec.structured[[as.character(t.now)]] <- as.matrix(heat.vec.struc)
  }
  projection.mat <- sapply(projection.mat, identity, simplify = "array")
  mean.mat <- sapply(mean.mat, identity, simplify = "array")

  sigma_ww.structured <- sapply(sigma_ww.structured, identity, simplify = "array")
  sigma_ww.empirical <- sapply(sigma_ww.empirical, identity, simplify = "array")
  sigma_zw.structured <- sapply(sigma_zw.structured, identity, simplify = "array")
  sigma_zw.empirical <- sapply(sigma_zw.empirical, identity, simplify = "array")
  sigma_zwzw.structured <- sapply(sigma_zwzw.structured, identity, simplify = "array")
  sigma_zwzw.empirical <- sapply(sigma_zwzw.empirical, identity, simplify = "array")
  zw_used <- sapply(zw_used, identity, simplify = "array")

  heat.vec.structured <- sapply(heat.vec.structured, identity, simplify = "array")

  return(list(
    Z.hat=Z.hat,
    heat.vec.structured=heat.vec.structured,
    projection.mat=projection.mat, mean.mat=mean.mat,
    sigma_ww.structured=sigma_ww.structured, sigma_ww.empirical=sigma_ww.empirical,
    sigma_zw.structured=sigma_zw.structured, sigma_zw.empirical=sigma_zw.empirical,
    sigma_zwzw.structured=sigma_zwzw.structured, sigma_zwzw.empirical=sigma_zwzw.empirical,
    zw_used=zw_used, epsilon=epsilon
  ))
}
