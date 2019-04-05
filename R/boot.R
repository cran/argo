#' bootstrap relative efficiency confidence interval
#'
#' This function is used to reproduce the ARGO bootstrap confidence interval
#'
#' @param pred_data A matrix that contains the truth vector and the predictions.
#' It can be data.frame or xts object
#'
#' @param model_good The model to evaluate, must be in the column names of pred_data
#' @param model_bench The model to compare to, must be in the column names of pred_data
#' @param l stationary bootstrap mean block length
#' @param N number of bootstrap samples
#' @param truth the column name of the truth
#' @param sim simulation method, pass to boot::tsboot
#' @param conf confidence level
#' @param type Must be one of "mse" (mean square error),
#' "mape" (mean absolute percentage error), or
#' "mae" (mean absolute error)
#'
#' @return A vector of point estimate and corresponding bootstrap confidence interval
#'
#' @examples
#' GFT_xts = xts::xts(exp(matrix(rnorm(1000), ncol=5)), order.by = Sys.Date() - (200:1))
#' names(GFT_xts) <- paste0("col", 1:ncol(GFT_xts))
#' names(GFT_xts)[1] <- "CDC.data"
#' bootstrap_relative_efficiency(
#'   pred_data = GFT_xts,
#'   model_good = "col2",
#'   model_bench = "col3",
#'   truth="CDC.data",
#'   N = 100
#' )
#' @export
bootstrap_relative_efficiency <-
  function(pred_data, model_good, model_bench, l=50, N = 1e4, truth="CDC.data",
           sim = "geom", conf = 0.95, type=c("mse", "mape", "mae", "mspe", "rmse", "rmspe")){

  if(type[1] %in% c("mse", "mspe")){
    err.fun <- function(tsb) {
      log(mean(tsb[,2]^2)) - log(mean(tsb[,1]^2))
    }
  }else if (type[1] %in% c("mae", "mape")){
    err.fun <- function(tsb) {
      log(mean(abs(tsb[,2]))) - log(mean(abs(tsb[,1])))
    }
  }else if(type[1] %in% c("rmse", "rmspe")){
    err.fun <- function(tsb) {
      (log(mean(tsb[,2]^2)) - log(mean(tsb[,1]^2)))/2
    }
  }




  pred_error <- xts::merge.xts(pred_data[,model_good]-pred_data[,truth],
                               pred_data[,model_bench]-pred_data[,truth])
  if(type[1] %in% c("mape", "mspe", "rmspe")){
    pred_error <- xts::merge.xts(
      (pred_data[,model_good]-pred_data[,truth])/pred_data[,truth],
      (pred_data[,model_bench]-pred_data[,truth])/pred_data[,truth])
  }
  pred_error <- na.omit(pred_error)
  err.fun(pred_error)
  boot.fit <- boot::tsboot(pred_error, err.fun, R = 1e4, l = l, sim = sim)


  CI <- boot::boot.ci(boot.fit, conf = conf, type = c("norm","basic","perc"))

  result <- exp(c(err.fun(pred_error),CI$basic[,4:5], CI$normal[,2:3],
                  CI$percent[,4:5]))
  names(result) <- c("point", "basic_lb", "basic_ub", "normal_lb", "normal_ub",
                     "percent_lb", "percent_ub")
  result
  }

#' wrapper for bootstrap relative efficiency confidence interval
#'
#' This function is used to wrap the \code{bootstrap_relative_efficiency},
#' taking vectorized arguments.
#'
#' @param pred_data A matrix that contains the truth vector and the predictions.
#' It can be data.frame or xts object
#'
#' @param period.all vector of the periods to evaluate relative efficiency
#' @param model_good The model to evaluate, must be in the column names of pred_data
#' @param bench.all vector of the models to compare to, must be in the column names of pred_data
#' @param l stationary bootstrap mean block length
#' @param N number of bootstrap samples
#' @param truth the column name of the truth
#' @param sim simulation method, pass to boot::tsboot
#' @param conf confidence level
#' @param type Must be one of "mse" (mean square error),
#' "mape" (mean absolute percentage error), or
#' "mae" (mean absolute error)
#'
#' @return A vector of point estimate and corresponding bootstrap confidence interval
#'
#' @examples
#' GFT_xts = xts::xts(exp(matrix(rnorm(500), ncol=5)), order.by = Sys.Date() - (100:1))
#' names(GFT_xts) <- paste0("col", 1:ncol(GFT_xts))
#' names(GFT_xts)[1] <- "CDC.data"
#' \donttest{
#' boot_re(
#'   pred_data = GFT_xts,
#'   period.all = c(paste0(zoo::index(GFT_xts)[1], "/", zoo::index(GFT_xts)[50]),
#'                  paste0(zoo::index(GFT_xts)[51], "/", zoo::index(GFT_xts)[100])),
#'   model_good = "col2",
#'   bench.all = c("col3", "col4"),
#'   type = "mse",
#'   truth="CDC.data",
#'   l = 5,
#'   N = 20
#' )
#' }
#' @export
boot_re <- function(pred_data, period.all, model_good, bench.all, type,
                    truth="CDC.data", l=50, N = 1e4, sim = "geom", conf = 0.95){
  rela_effi.all <- list()
  for(period in period.all){
    for(type.each in type){
      rela_effi <- lapply(bench.all, function(model_bench)
        bootstrap_relative_efficiency(
          pred_data[period], model_good, model_bench, truth=truth,
          l=l, sim = sim, type = type.each, conf=conf, N=N))
      rela_effi <- sapply(rela_effi, identity)
      rela_effi <- t(rela_effi)
      colnames(rela_effi) <- c("point estimate",
                               paste(rep(c("basic","normal","percent"), each=2),
                                     c("95% CI lower bound", "95% CI upper bound")))
      rownames(rela_effi) <- bench.all
      rela_effi.all[[paste0(type.each,";",period)]] <- rela_effi[,1:3]
    }
  }
  rela_effi.all
}
