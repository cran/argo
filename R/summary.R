#' performance summary of ARGO applied on CDC's ILI data
#'
#'
#' @param GFT_xts dataframe with all predicted values
#' @param model_names name of predicting models
#' @param legend_names legend for predicting models
#' @param periods vector of periods to zoom into
#' @param whole_period the whole period duration
#'
#' @examples
#' GFT_xts = xts::xts(exp(matrix(rnorm(1000), ncol=10)), order.by = Sys.Date() - (100:1))
#' names(GFT_xts) <- paste0("col", 1:10)
#' names(GFT_xts)[1] <- "CDC.data"
#' summary_argo(
#'   GFT_xts = GFT_xts,
#'   model_names = colnames(GFT_xts)[-1],
#'   legend_names = paste0(colnames(GFT_xts)[-1], "legend"),
#'   periods = c(paste0(zoo::index(GFT_xts)[1], "/", zoo::index(GFT_xts)[49]),
#'               paste0(zoo::index(GFT_xts)[50], "/", zoo::index(GFT_xts)[100])),
#'   whole_period="2009-03/"
#' )
#'
#' @import stats
#'
#' @return A list of summary tables for the input periods, including RMSE, MAE, MAPE, corr
#' @references
#' Yang, S., Santillana, M., & Kou, S. C. (2015). Accurate estimation of influenza epidemics using Google search data via ARGO. Proceedings of the National Academy of Sciences, \href{https://dx.doi.org/10.1073/pnas.1515373112}{doi: 10.1073/pnas.1515373112}.
#' Shaoyang Ning, Shihao Yang, S. C. Kou. Accurate Regional Influenza Epidemics Tracking Using Internet Search Data. Scientific Reports
#'
#' @export
summary_argo <- function(GFT_xts, model_names, legend_names, periods,
                         whole_period="2009-03/2015-10"){
  periods <- c(whole_period, periods)
  corr_by_period <- function(x,y, periods){
    sapply(periods, function(p) {
      re <- try(cor(x[p],y[p], use="everything"), silent = T)#"complete.obs"
      if(class(re)=="try-error")
        return(NaN)
      else
        return(re)
    })
  }

  RMSE_by_period <- function(x,y, periods){
    sapply(periods, function(p) {
      re <- try(sqrt(mean((x[p]-y[p])^2, na.rm=FALSE)),silent = T)
      if(class(re)=="try-error")
        return(NaN)
      else
        return(re)
    })
  }

  RMSPE_by_period <- function(x_truth,y, periods){
    x <- x_truth
    sapply(periods, function(p) {
      re <- try(sqrt(mean((x[p]-y[p])^2 / x_truth[p]^2, na.rm=FALSE)),silent = T)
      if(class(re)=="try-error")
        return(NaN)
      else
        return(re)
    })
  }

  ABSE_by_period <- function(x,y, periods){
    sapply(periods, function(p) {
      re <- try(mean(abs(x[p]-y[p]), na.rm=FALSE),silent = T)
      if(class(re)=="try-error")
        return(NaN)
      else
        return(re)
    })
  }

  corr_martingale_diff <- function(x,y, periods){
    sapply(periods, function(p) {
      re <- try(cor(diff(x[p])[-1],diff(y[p])[-1], use="everything"), silent = T)#"complete.obs"
      if(class(re)=="try-error")
        return(NaN)
      else
        return(re)
    })
  }


  MAPE_by_period <- function(x_truth,y, periods){
    sapply(periods, function(p) {
      re <- try(mean(abs(x_truth[p]-y[p])/x_truth[p], na.rm=FALSE),silent = T)
      if(class(re)=="try-error")
        return(NaN)
      else
        return(re)
    })
  }

  MaxAPE_by_period <- function(x_truth,y, periods){
    sapply(periods, function(p) {
      tryCatch({
        max(abs(x_truth[p]-y[p])/x_truth[p], na.rm=FALSE)
      }, warning = function(w) {
        return(NaN)
      }, error = function(e) {
        return(NaN)
      })
    })
  }

  bias_by_period <- function(x_truth,y, periods){
    sapply(periods, function(p) {
      tryCatch({
        mean(y[p] - x_truth[p], na.rm=FALSE)
      }, warning = function(w) {
        return(NaN)
      }, error = function(e) {
        return(NaN)
      })
    })
  }

  relativebias_by_period <- function(x_truth,y, periods){
    sapply(periods, function(p) {
      tryCatch({
        mean((y[p] - x_truth[p])/x_truth[p], na.rm=FALSE)
      }, warning = function(w) {
        return(NaN)
      }, error = function(e) {
        return(NaN)
      })
    })
  }

  percent_up_by_period <- function(x_truth,y, periods){
    sapply(periods, function(p) {
      tryCatch({
        mean(y[p] > x_truth[p], na.rm=FALSE)
      }, warning = function(w) {
        return(NaN)
      }, error = function(e) {
        return(NaN)
      })
    })
  }

  mse_up_by_period <- function(x_truth,y, periods){
    pred_error <- y-x_truth
    sapply(periods, function(p) {
      tryCatch({
        pred_error_this <- as.numeric(pred_error[p])
        pred_error_this <- pred_error_this[pred_error_this>0]
        mean(pred_error_this^2, na.rm=FALSE)
      }, warning = function(w) {
        return(NaN)
      }, error = function(e) {
        return(NaN)
      })
    })
  }

  mse_down_by_period <- function(x_truth,y, periods){
    pred_error <- y-x_truth
    sapply(periods, function(p) {
      tryCatch({
        pred_error_this <- as.numeric(pred_error[p])
        pred_error_this <- pred_error_this[pred_error_this<0]
        mean(pred_error_this^2, na.rm=FALSE)
      }, warning = function(w) {
        return(NaN)
      }, error = function(e) {
        return(NaN)
      })
    })
  }


  periods_short <- periods



  corr_period_tab <-
    sapply(model_names, function(name)
      corr_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  corr_period_tab <- matrix(corr_period_tab, ncol=length(model_names))
  colnames(corr_period_tab) <- legend_names
  rownames(corr_period_tab) <- periods_short

  RMSE_period_tab <-
    sapply(model_names, function(name)
      RMSE_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  RMSE_period_tab <- matrix(RMSE_period_tab, ncol=length(model_names))
  colnames(RMSE_period_tab) <- legend_names
  rownames(RMSE_period_tab) <- periods_short

  ABSE_period_tab <-
    sapply(model_names, function(name)
      ABSE_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  ABSE_period_tab <- matrix(ABSE_period_tab, ncol=length(model_names))
  colnames(ABSE_period_tab) <- legend_names
  rownames(ABSE_period_tab) <- periods_short

  corr_martingale_diff_period_tab <-
    sapply(model_names, function(name)
      corr_martingale_diff(GFT_xts$CDC.data, GFT_xts[,name], periods))
  corr_martingale_diff_period_tab <- matrix(corr_martingale_diff_period_tab, ncol=length(model_names))
  colnames(corr_martingale_diff_period_tab) <- legend_names
  rownames(corr_martingale_diff_period_tab) <- periods_short

  MAPE_period_tab <-
    sapply(model_names, function(name)
      MAPE_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  MAPE_period_tab <- matrix(MAPE_period_tab, ncol=length(model_names))
  colnames(MAPE_period_tab) <- legend_names
  rownames(MAPE_period_tab) <- periods_short

  MaxAPE_period_tab <-
    sapply(model_names, function(name)
      MaxAPE_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  MaxAPE_period_tab <- matrix(MaxAPE_period_tab, ncol=length(model_names))
  colnames(MaxAPE_period_tab) <- legend_names
  rownames(MaxAPE_period_tab) <- periods_short

  RMSPE_period_tab <-
    sapply(model_names, function(name)
      RMSPE_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  RMSPE_period_tab <- matrix(RMSPE_period_tab, ncol=length(model_names))
  colnames(RMSPE_period_tab) <- legend_names
  rownames(RMSPE_period_tab) <- periods_short

  bias_period_tab <-
    sapply(model_names, function(name)
      bias_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  bias_period_tab <- matrix(bias_period_tab, ncol=length(model_names))
  colnames(bias_period_tab) <- legend_names
  rownames(bias_period_tab) <- periods_short

  percent_up_period_tab <-
    sapply(model_names, function(name)
      percent_up_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  percent_up_period_tab <- matrix(percent_up_period_tab, ncol=length(model_names))
  colnames(percent_up_period_tab) <- legend_names
  rownames(percent_up_period_tab) <- periods_short

  mse_up_period_tab <-
    sapply(model_names, function(name)
      mse_up_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  mse_up_period_tab <- matrix(mse_up_period_tab, ncol=length(model_names))
  colnames(mse_up_period_tab) <- legend_names
  rownames(mse_up_period_tab) <- periods_short

  mse_down_period_tab <-
    sapply(model_names, function(name)
      mse_down_by_period(GFT_xts$CDC.data, GFT_xts[,name], periods))
  mse_down_period_tab <- matrix(mse_down_period_tab, ncol=length(model_names))
  colnames(mse_down_period_tab) <- legend_names
  rownames(mse_down_period_tab) <- periods_short

  RMSE_period_tab_normed <-
    format(round(RMSE_period_tab/RMSE_period_tab[,ncol(RMSE_period_tab)], 3),
           nsmall=3)
  RMSE_period_tab_print <-
    matrix(paste0(RMSE_period_tab_normed, " (",format(round(RMSE_period_tab, 3), nsmall=3),")"),
           nrow=nrow(RMSE_period_tab),
           ncol=ncol(RMSE_period_tab),
           dimnames=dimnames(RMSE_period_tab))
  RMSE_period_tab_print[,-ncol(RMSE_period_tab_print)] <-
    RMSE_period_tab_normed[,-ncol(RMSE_period_tab)]

  RMSPE_period_tab_normed <-
    format(round(RMSPE_period_tab/RMSPE_period_tab[,ncol(RMSPE_period_tab)], 3),
           nsmall=3)
  RMSPE_period_tab_print <-
    matrix(paste0(RMSPE_period_tab_normed, " (",format(round(RMSPE_period_tab, 3), nsmall=3),")"),
           nrow=nrow(RMSPE_period_tab),
           ncol=ncol(RMSPE_period_tab),
           dimnames=dimnames(RMSPE_period_tab))
  RMSPE_period_tab_print[,-ncol(RMSPE_period_tab_print)] <-
    RMSPE_period_tab_normed[,-ncol(RMSPE_period_tab)]


  ABSE_period_tab_normed <-
    format(round(ABSE_period_tab/ABSE_period_tab[,ncol(ABSE_period_tab)], 3),
           nsmall=3)
  ABSE_period_tab_print <-
    matrix(paste0(ABSE_period_tab_normed, " (",format(round(ABSE_period_tab, 3), nsmall=3),")"),
           nrow=nrow(ABSE_period_tab),
           ncol=ncol(ABSE_period_tab),
           dimnames=dimnames(ABSE_period_tab))
  ABSE_period_tab_print[,-ncol(ABSE_period_tab_print)] <-
    ABSE_period_tab_normed[,-ncol(ABSE_period_tab_normed)]


  MAPE_period_tab_normed <-
    format(round(MAPE_period_tab/MAPE_period_tab[,ncol(MAPE_period_tab)], 3),
           nsmall=3)
  MAPE_period_tab_print <-
    matrix(paste0(MAPE_period_tab_normed, " (",format(round(MAPE_period_tab, 3), nsmall=3),")"),
           nrow=nrow(MAPE_period_tab),
           ncol=ncol(MAPE_period_tab),
           dimnames=dimnames(MAPE_period_tab))
  MAPE_period_tab_print[,-ncol(MAPE_period_tab_print)] <-
    MAPE_period_tab_normed[,-ncol(MAPE_period_tab_normed)]


  add_bold_to_print_tab <- function(tab_print, tab_orig, smallerbetter=TRUE){
    if(!smallerbetter)
      tab_orig <- -tab_orig
    for(i in 1:ncol(tab_orig)){
      id <- order(tab_orig[,i])[1]
      tab_print[id,i] <- paste0("\\textbf{",tab_print[id,i],"}")
    }
    return(tab_print)
  }

  corr_period_tab_print <-
    matrix(format(round(corr_period_tab, 3), nsmall=3),
           nrow=nrow(corr_period_tab),
           ncol=ncol(corr_period_tab),
           dimnames=dimnames(corr_period_tab))

  corr_martingale_diff_period_tab_print <-
    matrix(format(round(corr_martingale_diff_period_tab, 3), nsmall=3),
           nrow=nrow(corr_martingale_diff_period_tab),
           ncol=ncol(corr_martingale_diff_period_tab),
           dimnames=dimnames(corr_martingale_diff_period_tab))

  corr_period_tab_print <-
    add_bold_to_print_tab(t(corr_period_tab_print), t(corr_period_tab), FALSE)

  corr_martingale_diff_period_tab_print <-
    add_bold_to_print_tab(t(corr_martingale_diff_period_tab_print), t(corr_martingale_diff_period_tab), FALSE)

  RMSE_period_tab_print <-
    add_bold_to_print_tab(t(RMSE_period_tab_print), t(RMSE_period_tab))

  RMSPE_period_tab_print <-
    add_bold_to_print_tab(t(RMSPE_period_tab_print), t(RMSPE_period_tab))

  ABSE_period_tab_print <-
    add_bold_to_print_tab(t(ABSE_period_tab_print), t(ABSE_period_tab))

  MAPE_period_tab_print <-
    add_bold_to_print_tab(t(MAPE_period_tab_print), t(MAPE_period_tab))


  return(list(corr=t(corr_period_tab),
              rmse=t(RMSE_period_tab),
              mse=t(RMSE_period_tab^2),
              abse=t(ABSE_period_tab),
              mape=t(MAPE_period_tab),
              rmspe=t(RMSPE_period_tab),
              maxape=t(MaxAPE_period_tab),
              bias=t(bias_period_tab),
              percent_up=t(percent_up_period_tab),
              mse_up=t(mse_up_period_tab),
              mse_down=t(mse_down_period_tab),
              corr_diff=t(corr_martingale_diff_period_tab),
              corr_print=(corr_period_tab_print),
              rmse_print=(RMSE_period_tab_print),
              rmspe_print=RMSPE_period_tab_print,
              abse_print=(ABSE_period_tab_print),
              mape_print=(MAPE_period_tab_print),
              corr_diff_print=(corr_martingale_diff_period_tab_print)
  ))

}

