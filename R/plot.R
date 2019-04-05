plot_xts <- function (x, y = NULL, type = "l", auto.grid = TRUE, major.ticks = "auto",
                      minor.ticks = TRUE, major.format = TRUE, bar.col = "grey",
                      candle.col = "white", ann = TRUE, axes = TRUE, add = FALSE, ...)
{
  series.title <- deparse(substitute(x))
  ep <- axTicksByTime(x, major.ticks, format.labels = major.format)
  otype <- type
  if (NCOL(x) > 1)
    warning("only the univariate series will be plotted")
  if (is.null(y))
    xycoords <- xy.coords(.index(x), x[, 1])
  if(add){
    lines(xycoords$x, xycoords$y, type = type, ann = FALSE,
          ...)
    return(invisible(NULL))
  }else{
    plot(xycoords$x, xycoords$y, type = type, axes = FALSE, ann = FALSE,
         ...)
  }
  if (auto.grid) {
    abline(v = xycoords$x[ep], col = "grey", lty = 4)
    grid(NA, NULL)
  }

  dots <- list(...)
  if (axes) {
    if (minor.ticks)
      axis(1, at = xycoords$x, labels = FALSE, col = "#BBBBBB",
           ...)
    axis(1, at = xycoords$x[ep], labels = names(ep), las = 1,
         lwd = 1, mgp = c(3, 2, 0), ...)
    axis(2, ...)
  }
  box()
  if (!"main" %in% names(dots))
    title(main = series.title)
  do.call("title", list(...))
}

#' Time series plot of ARGO applied on CDC's ILI data
#'
#' This function is used to reproduce the ARGO plot.
#'
#' @param GFT_xts dataframe with all predicted values
#' @param GC_GT_cut_date cutting date for switching datasets
#' @param model_names name of predicting models
#' @param legend_names legend for predicting models
#' @param zoom_periods vector of periods to zoom into
#'
#' @importFrom zoo index
#' @import graphics
#'
#' @return a graph on the default plot window
#'
#' @examples
#' GFT_xts = xts::xts(exp(matrix(rnorm(1000), ncol=5)), order.by = Sys.Date() - (200:1))
#' names(GFT_xts) <- paste0("col", 1:ncol(GFT_xts))
#' names(GFT_xts)[1] <- "CDC.data"
#' zoom_periods = c()
#' for (i in 0:5){
#'   zoom_periods = c(
#'     zoom_periods,
#'     paste0(zoo::index(GFT_xts)[i*30+1], "/", zoo::index(GFT_xts)[i*30+30])
#'   )
#' }
#' plot_argo(
#'   GFT_xts = GFT_xts,
#'   GC_GT_cut_date = zoo::index(GFT_xts)[50],
#'   model_names = colnames(GFT_xts)[-1],
#'   legend_names = paste0(colnames(GFT_xts)[-1], "legend"),
#'   zoom_periods = zoom_periods
#' )
#'
#' @export
plot_argo <- function(GFT_xts, GC_GT_cut_date, model_names, legend_names, zoom_periods){
  par_orig <- par()[c("mar","mfrow")]
  cutting_dates <- which(index(GFT_xts) %in% c(as.Date("2010-05-22"), GC_GT_cut_date-1))

  color_code <- c("red","limegreen", "blue", "goldenrod3","grey48")
  term_change_col <- "khaki"
  cdc_lwd <- 2
  line_width <- c(1.5,1,1,1,1)
  line_type <- c(1,1,1,1,2)

  layout(rbind(rep(1,3),rep(2,3),3:5), heights=c(2.5/3.5*6,6/3.5,5))
  par(mar=c(4.1, 2.1, 1.1, 1))
  plot_xts(GFT_xts$CDC.data, ylim=c(-0.1,8.5), type='n', auto.grid = FALSE, main=NA)

  abline(v=xts::.index(GFT_xts)[cutting_dates[1]], col=term_change_col, lwd=2, lty=4)
  rect(xts::.index(GFT_xts)[cutting_dates[2]], -10, xts::.index(GFT_xts)[nrow(GFT_xts)], 100,
       border=NA, col="lightgoldenrodyellow")
  box()
  plot_xts(GFT_xts$CDC.data, lwd=cdc_lwd, add=TRUE)
  for(i in length(model_names):1){
    plot_xts(GFT_xts[,model_names[i]], col=color_code[i], lwd=line_width[i], lty=line_type[i], add=TRUE)
  }
  legend(xts::.index(GFT_xts)[75], 8.6,c("CDC's ILI activity level (weighted)",legend_names),
         col=c(1, color_code),
         lwd=c(cdc_lwd, line_width),
         lty=c(1,line_type), cex=1.5)

  arrows(xts::.index(GFT_xts)[cutting_dates[1]-1],0,xts::.index(GFT_xts)[1],0,
         length=0.1, col="grey23", code=3)
  text(xts::.index(GFT_xts)[cutting_dates[1]%/%2+1], 0,
       "Google Correlate", pos=3, offset=0.3, cex=1.25)
  text(xts::.index(GFT_xts)[cutting_dates[1]%/%2+1], 0,
       "Search terms identified on 2009-03-28", pos=1, offset=0.3, cex=1.25)
  arrows(xts::.index(GFT_xts)[cutting_dates[1]+1],0,xts::.index(GFT_xts)[cutting_dates[2]-1],0,
         length=0.1, col="grey23", code=3)
  text(xts::.index(GFT_xts)[mean(cutting_dates)], 0,
       "Google Correlate", pos=3, offset=0.3, cex=1.25)
  text(xts::.index(GFT_xts)[mean(cutting_dates)], 0,
       "Search terms identified on 2010-05-22", pos=1, offset=0.3, cex=1.25)
  arrows(xts::.index(GFT_xts)[cutting_dates[2]+1],0,xts::.index(GFT_xts)[nrow(GFT_xts)-1],0,
         length=0.1, col="grey23", code=3)
  text(xts::.index(GFT_xts)[(cutting_dates[2]+nrow(GFT_xts))%/%2+1], 0,
       "Google Trend", pos=3, offset=0.3, cex=1.25)
  text(xts::.index(GFT_xts)[(cutting_dates[2]+nrow(GFT_xts))%/%2+1], 8,
       "Search terms\nidentified on\n2010-05-22", pos=1, offset=0, cex=1.25)

  old_mar <- par(mar=c(0.5,2.1,0.5,1))


  plot_xts(GFT_xts$CDC.data-GFT_xts[,model_names[i]], ylim=c(-2,2),main=NA,type='n',
           xlab=NA, xaxt='n',auto.grid = FALSE)


  abline(v=xts::.index(GFT_xts)[cutting_dates[1]], col=term_change_col, lwd=2, lty=4)

  rect(xts::.index(GFT_xts)[cutting_dates[2]], -10, xts::.index(GFT_xts)[nrow(GFT_xts)], 100,
       border=NA, col="lightgoldenrodyellow")

  box()
  mtext("prediction error")
  abline(h=0,lty=3,lwd=cdc_lwd)
  for(i in length(model_names):1){
    plot_xts(-GFT_xts$CDC.data+GFT_xts[,model_names[i]], col=color_code[i],
             lwd=line_width[i], lty=line_type[i], add=TRUE)
  }




  zoom_periods_dates <- strsplit(zoom_periods, "/")
  zoom_periods_dates <- lapply(zoom_periods_dates, as.Date)

  legend_period <- lapply(1:length(zoom_periods_dates), function(i) {
    if(i == 1)
      l <- "H1N1 Flu outbreak"
    else if(i == length(zoom_periods_dates))
      l <- paste0(2008+i, "-",9+i,"\nFlu Season")
    else
      l <- paste0(2008+i, "-",9+i,"\nFlu Season")

    d <- paste(format(zoom_periods_dates[[i]], "%m/%d/%y"), collapse="\n     --\n")
    return(c(l, d))
  })

  par(mar=c(3,2.1,2.1,1))
  zoom_id <- c(1,4,6)
  for(k in 1:length(zoom_id)){
    period_counter <- zoom_id[k]
    plot_period <- zoom_periods[period_counter]
    plot_xts(GFT_xts$CDC.data[plot_period], ylim=c(1.2,8.5),
         main=NA, type='n', auto.grid = FALSE)
    plot_xts(GFT_xts$CDC.data[plot_period], lwd=cdc_lwd, add=TRUE)
    for(i in length(model_names):1){
      plot_xts(GFT_xts[,model_names[i]][plot_period], col=color_code[i], lwd=line_width[i], lty=line_type[i], add=TRUE)
    }
    mtext(paste0("(",letters[k],")"))
    if(period_counter == 1){
      legend("topleft", rev(legend_period[[period_counter]]), bty='n', cex=1.5)
    }else{
      legend("topright", legend_period[[period_counter]][1], bty='n', cex=1.5)
      legend("topleft", legend_period[[period_counter]][2], bty='n', cex=1.5)
    }
  }
  par(par_orig)
  invisible(NULL)
}

#' Heatmap plot of ARGO coefficients applied on CDC's ILI data
#'
#'
#' @param argo_coef The coefficient matrix
#' @param lim the limit to truncate for large coefficients for better presentation
#' @param na.grey whether to plot grey for NA values
#' @param scale margin scale
#'
#' @importFrom zoo index
#' @import graphics
#' @import grDevices
#'
#' @return a graph on the default plot window
#' @examples
#' cor_coef <- matrix(runif(100, -1, 1), ncol=10)
#' colnames(cor_coef) <- as.character(Sys.Date() - 10:1)
#' rownames(cor_coef) <- paste0("row", 1:10)
#' pdf(file.path(tempdir(), "heatmap_argo.pdf"), height=11,width=12)
#' heatmap_argo(cor_coef)
#' dev.off()
#'
#' @export
heatmap_argo <- function(argo_coef, lim=0.1, na.grey=TRUE, scale=1){
  argo_coef_plot <- argo_coef
  argo_coef_plot[which(argo_coef_plot > lim)] <- lim
  argo_coef_plot[which(argo_coef_plot < -lim)] <- -lim
  n_col <- 201
  col.rng.id <- round(1+200*(lim+range(argo_coef_plot, na.rm = TRUE))/(2*lim))


  if(na.grey){
    layout(matrix(c(1,3,2,3),2), heights=c(1.2,24),widths=c(1,5))
  }else{
    layout(1:2, heights=c(1.8,24))
  }

  if(na.grey){
    par(mar = c(0.8,1.3,1.5,1.4)*scale)
    image(x=1, z=matrix(0,1,1), col="grey88",
          xlab=NA, ylab=NA, yaxt='n', xaxt='n')
    axis(3, at=1,line=-1,cex.axis=0.75,
         labels="N/A",
         lwd.ticks=0, lwd=0)
  }


  col.rng <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(n = 201)
  if(na.grey){
    par(mar = c(0.8,2.1,1.5,3.1)*scale)
  }else{
    par(mar = c(0.8,3.1,1.5,3.1)*scale)
  }
  image(x=1:length(col.rng), z=matrix(seq(-lim, lim, length=length(col.rng)), ncol=1), col=col.rng,
        xlab=NA, ylab=NA, yaxt='n', xaxt='n')
  axis(1, at=seq(1, n_col, length=5),line=-1,cex.axis=0.75,
       labels=c(paste("<",-lim), -lim/2, 0, lim/2, paste(">", lim)),
       lwd.ticks=0, lwd=0)
  axis(3, at=seq(1, n_col, length=3),line=-1,cex.axis=0.75,
       labels=paste(c("Negative", "Zero", "Positive"), "coefficient"),
       lwd.ticks=0, lwd=0)

  bg <- argo_coef_plot
  bg[,] <- 0

  par(mar = c(2.1,5.1,1.1,2.1))
  colname2date <- function(str) {
    tryCatch(
      as.Date(str),
      error = function(e) as.Date(paste("1",str), c("%d %b %Y")))
  }
  par(mar = c(2.1,5.1,1.1,2.1)*scale)
  if(na.grey){
    image(x=colname2date(colnames(argo_coef_plot)), y=1:nrow(argo_coef_plot),
          z=t(bg), col="grey88", xlab=NA, ylab=NA, yaxt='n')
    image(x=colname2date(colnames(argo_coef_plot)), y=1:nrow(argo_coef_plot),add=TRUE,
          z=t(argo_coef_plot)[,nrow(argo_coef):1], col=col.rng[col.rng.id[1]:col.rng.id[2]],
          xlab=NA, ylab=NA, yaxt='n')
  }else{
    image(x=colname2date(colnames(argo_coef_plot)), y=1:nrow(argo_coef_plot),
          z=t(argo_coef_plot)[,nrow(argo_coef):1], col=col.rng[col.rng.id[1]:col.rng.id[2]],
          xlab=NA, ylab=NA, yaxt='n')
  }
  axis(2, at=nrow(argo_coef_plot):1, labels=rownames(argo_coef_plot),
       tick=FALSE, cex.axis=0.4, las=2, line=-0.75)
  axis(1, at=seq(as.Date("2004-01-01"), as.Date("2018-01-01"), by="month"),labels=NA,
       lwd.ticks=0.5, lwd=0, tck=-0.01)
}

#' Heatmap plot of correlation matrix
#'
#'
#' @param cor_heat The coefficient matrix to draw heatmap
#' @param lim the limit to truncate for large coefficients for better presentation
#'
#' @importFrom zoo index
#' @import graphics
#'
#' @return a graph on the default plot window
#'
#' @examples
#' cor_coef <- matrix(runif(100, -1, 1), ncol=10)
#' colnames(cor_coef) <- paste0("col", 1:10)
#' rownames(cor_coef) <- paste0("row", 1:10)
#' heatmap_cor(cor_coef)
#'
#' @export
heatmap_cor <- function(cor_heat, lim=1){
  cor_heat_plot <- cor_heat
  cor_heat_plot[which(cor_heat_plot > lim)] <- lim
  cor_heat_plot[which(cor_heat_plot < -lim)] <- -lim
  n_col <- 201
  col.rng.id <- round(1+200*(lim+range(cor_heat_plot, na.rm = TRUE))/(2*lim))

  layout(matrix(1:2), heights=c(1,8))

  n_col <- 201
  col.rng <- colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(n = 201)
  par(mar = c(0.8,3.1,1.5,3.1))
  image(x=1:length(col.rng), z=matrix(seq(-lim, lim, length=length(col.rng)), ncol=1), col=col.rng,
        xlab=NA, ylab=NA, yaxt='n', xaxt='n')
  axis(1, at=seq(1, n_col, length=5),line=-1,cex.axis=0.75,
       labels=c(paste(-lim), -lim/2, 0, lim/2, paste(lim)),
       lwd.ticks=0, lwd=0)
  axis(3, at=seq(1, n_col, length=3),line=-1,cex.axis=0.75,
       labels=paste(c("Negative", "Zero", "Positive"), "value"),
       lwd.ticks=0, lwd=0)

  par(mar = c(4.1,4.1,1.1,2.1))
  image(x=1:ncol(cor_heat_plot), y=1:nrow(cor_heat_plot),
        z=t(cor_heat_plot)[,nrow(cor_heat):1], col=col.rng[col.rng.id[1]:col.rng.id[2]],
        xlab=NA, ylab=NA, yaxt='n', xaxt='n')
  axis(2, at=nrow(cor_heat_plot):1, labels=rownames(cor_heat_plot),
       tick=FALSE, cex.axis=0.4, las=2, line=-0.75)
  axis(1, at=1:ncol(cor_heat_plot), labels=colnames(cor_heat_plot),
       tick=FALSE, cex.axis=0.4, las=2, line=-0.75)
}

