################################################################################
# Author: Shihao Yang                                                          #
# Date: Mar 16, 2019                                                           #
# Main script for argo package, also served as an example of package usage     #
#                                                                              #
# Required package: argo, xts, glmnet, boot, xtable                            #
################################################################################


#' main function for argo
#'
#' main function that reproduce the results in ARGO paper
#' @param save.folder output folder to save graphics. If NULL then do not output graphics.
#'
#' @examples
#' \donttest{
#' argo_main()
#' }
#'
#' @export
argo_main <- function(save.folder=NULL){
  #### loading data ####
  all_data <- load_data()
  tail(all_data$GC09)
  tail(all_data$GC10)
  tail(all_data$GT)
  tail(all_data$CDC)
  tail(all_data$GFT)

  #### start of meat ####
  santillana_etal <- list()
  our_argo <- list()

  # The script takes around 10 minutes to run, generating the entire results
  # Real-case nowcast is available when dtname == "GT"
  for(dtname in c("GC10","GC09","GT")){
    exogen <- all_data[[dtname]]
    yx_merged <- merge(all_data$CDC, exogen, join = "right")
    y <- yx_merged$CDC.data
    santillana_etal[[dtname]] <-
      argo(y, exogen = exogen,
           alpha = 1, N_lag = NULL, use_all_previous = TRUE,  N_training = 104)

    our_argo[[dtname]] <-
      argo(logit(y / 100), exogen = log((exogen + 0.5) / 100),
           alpha = 1, N_lag=1:52, use_all_previous = FALSE, N_training = 104)
  }

  ar3 <- argo(all_data$CDC, alpha = NA, use_all_previous = FALSE, N_lag=1:3,
              N_training = 104)

  ar3_gft <- argo(all_data$CDC[zoo::index(all_data$GFT)], exogen = all_data$GFT,
                  alpha = NA, use_all_previous = FALSE, N_lag=1:3,
                  N_training = 104)


  #### presentation and plot ####
  h1n1_start <- which(colnames(our_argo$GC09$coef)=="2009-03-28")
  h1n1_end <- which(colnames(our_argo$GC09$coef)=="2010-05-15")
  ARGO_coef_GC09_h1n1 <- our_argo$GC09$coef[, h1n1_start:h1n1_end]

  post2009_start <- which(colnames(our_argo$GC10$coef)=="2009-03-28")
  post2009_end <- ncol(our_argo$GC10$coef)
  ARGO_coef_GC10_post2009 <- our_argo$GC10$coef[,post2009_start:post2009_end]

  GT_coef_start <- post2009_end + 1
  GT_coef_end <- ncol(our_argo$GT$coef)
  GT_coef_id <- GT_coef_start:GT_coef_end

  ARGO_coef_blend <- merge(ARGO_coef_GC10_post2009, our_argo$GT$coef[,GT_coef_id],
                           by = "row.names", all=TRUE)
  rownames(ARGO_coef_blend) <- ARGO_coef_blend[,1]
  ARGO_coef_blend <- data.matrix(ARGO_coef_blend[,-1], rownames.force = TRUE)
  ARGO_coef_blend <- ARGO_coef_blend[rownames(ARGO_coef_GC10_post2009),]

  ARGO_coef_blend <- ARGO_coef_blend[,setdiff(colnames(ARGO_coef_blend),colnames(ARGO_coef_GC09_h1n1))]
  ARGO_coef_blend <- merge(ARGO_coef_GC09_h1n1, ARGO_coef_blend, by = "row.names", all=TRUE)
  rownames(ARGO_coef_blend) <- ARGO_coef_blend[,1]
  ARGO_coef_blend <- data.matrix(ARGO_coef_blend[,-1], rownames.force = TRUE)
  phrases_both <- intersect(rownames(our_argo$GC10$coef), rownames(our_argo$GC09$coef))
  phrases_GC10 <- setdiff(rownames(our_argo$GC10$coef), rownames(our_argo$GC09$coef))
  phrases_GC09 <- setdiff(rownames(our_argo$GC09$coef), rownames(our_argo$GC10$coef))
  ARGO_coef_blend <- ARGO_coef_blend[c(phrases_both, phrases_GC10, phrases_GC09),]


  pred_xts_blend <- merge(all_data$CDC, logit_inv(our_argo$GT$pred)*100, all_data$GFT,
                          santillana_etal$GT$pred, ar3_gft$pred, ar3$pred, all=FALSE)
  names(pred_xts_blend) <- c("CDC.data", "ARGO", "GFT", "Santillana", "GFT_AR3", "AR3")

  pred_xts_blend$naive <- c(NA, as.numeric(pred_xts_blend$CDC.data[-nrow(pred_xts_blend)]))

  pred_xts_blend$ARGO[zoo::index(our_argo$GC10$pred)] <- logit_inv(our_argo$GC10$pred)*100
  pred_xts_blend$ARGO["/2010-05-23"] <- logit_inv(our_argo$GC09$pred["/2010-05-23"])*100

  pred_xts_blend$Santillana[zoo::index(santillana_etal$GC10$pred)] <- santillana_etal$GC10$pred
  pred_xts_blend$Santillana["/2010-05-23"] <- santillana_etal$GC09$pred["/2010-05-23"]

  pred_xts_blend <- pred_xts_blend["2009-03-29/"]

  model_names <- c("ARGO", "GFT", "Santillana", "GFT_AR3", "AR3", "naive")
  legend_names <- c("ARGO", "GFT (Oct 2014)", "Santillana et al. (2014)", "GFT+AR(3)","AR(3)", "Naive")

  zoom_periods <- c("2009-03-29/2009-12-27",
                    "2010-10-03/2011-05-22",
                    "2011-10-02/2012-05-20",
                    "2012-09-30/2013-05-19",
                    "2013-09-29/2014-05-18",
                    "2014-09-28/2015-05-17")

  GC_GT_cut_date <- tail(zoo::index(our_argo$GC10$pred),1)+1

  if(!is.null(save.folder)){
    pdf(file.path(save.folder, "final_plot.pdf"), height=11,width=12)
    plot_argo(pred_xts_blend, GC_GT_cut_date, model_names[1:5], legend_names[1:5], zoom_periods)
    dev.off()
    pdf(file.path(save.folder, "final_plot_job_talk.pdf"), height=11,width=12)
    plot_argo(pred_xts_blend, GC_GT_cut_date, model_names[1:5], c("ARGO", "GFT (Oct 2014)", "Google Data Only", "GFT+AR(3)","AR(3)"), zoom_periods)
    dev.off()
    pdf(file.path(save.folder, "final_plot_argo_gft.pdf"), height=11,width=12)
    plot_argo(pred_xts_blend, GC_GT_cut_date, c("ARGO", "GFT"), c("ARGO", "GFT (Oct 2014)"), zoom_periods)
    dev.off()
    pdf(file.path(save.folder, "heatmap.pdf"),width=4, height=10)
    heatmap_argo(ARGO_coef_blend, 0.1)
    dev.off()
  }

  all_tab <- summary_argo(pred_xts_blend, model_names, legend_names, zoom_periods, "2009/2015")

  corr_header <- all_tab$corr_print[1,,drop=F]
  rownames(corr_header) <- "\\textbf{Correlation}\\hfill\\vadjust{}"
  corr_header[,] <- NA

  rmse_header <- all_tab$rmse_print[1,,drop=F]
  rownames(rmse_header) <- "\\textbf{RMSE}\\hfill\\vadjust{}"
  rmse_header[,] <- NA

  abse_header <- all_tab$abse_print[1,,drop=F]
  rownames(abse_header) <- "\\textbf{MAE}\\hfill\\vadjust{}"
  abse_header[,] <- NA

  mape_header <- all_tab$mape_print[1,,drop=F]
  rownames(mape_header) <- "\\textbf{MAPE}\\hfill\\vadjust{}"
  mape_header[,] <- NA

  corr_diff_header <- all_tab$corr_diff_print[1,,drop=F]
  rownames(corr_diff_header) <- "\\textbf{Corr. of increment}\\hfill\\vadjust{}"
  corr_diff_header[,] <- NA


  big_tab_print <- rbind(
    rmse_header, all_tab$rmse_print,
    abse_header, all_tab$abse_print,
    mape_header, all_tab$mape_print,
    corr_header, all_tab$corr_print,
    corr_diff_header, all_tab$corr_diff_print)



  blank_suffix <- sapply(1:(nrow(big_tab_print)/(nrow(all_tab$rmse_print)+1)), function(i)
    paste(rep(" ",i), collapse = ""))
  rownames(big_tab_print) <-
    paste0(rownames(big_tab_print), rep(blank_suffix, each=(nrow(all_tab$rmse_print)+1)))

  print(xtable::xtable(big_tab_print), sanitize.text.function=identity,
        sanitize.rownames.function=identity)

  rela_effi <- sapply(model_names[-1], function(model_bench)
    bootstrap_relative_efficiency(
      na.omit(pred_xts_blend), model_names[1], model_bench, l=52, sim = "geom"))
  rela_effi <- t(rela_effi)
  colnames(rela_effi) <- c("point estimate",
                           paste(rep(c("basic","normal","percent"), each=2),
                                 c("95% CI lower bound", "95% CI upper bound")))
  rownames(rela_effi) <- legend_names[-1]
  print(xtable::xtable(rela_effi))
}


#' main function for argo2
#'
#' main function that reproduce the results in ARGO2 paper
#'
#' @param gt.folder folder with Google Trends files, which should be thousands of csv file such as "US-MA_fever cough.csv" or "US-NY_cold or flu.csv"
#' @param ili.folder folder with ILINet data files: "ILINet_nat.csv" and "ILINet_regional.csv"
#' @param population.file file path to population csv file
#' @param gft.file file path to Google Flu Trends csv file
#' @param save.folder output folder to save graphics. If NULL then do not output graphics.
#'
#' @references Shaoyang Ning, Shihao Yang, S. C. Kou. Accurate Regional Influenza Epidemics Tracking Using Internet Search Data. Scientific Reports
#' @import xts
#' @import glmnet
#' @import parallel
#' @import boot
#'
#' @examples
#'
#' \dontrun{
#' download.file("https://zenodo.org/records/17681160/files/GT2016-10-24.zip",
#' file.path(tempdir(), "gt2016-10-24.zip"))
#' unzip(file.path(tempdir(), "gt2016-10-24.zip"), exdir = tempdir())
#' gt.folder <- file.path(tempdir(), "2016-10-19")
#' argo2_main(
#'   gt.folder=gt.folder,
#'   ili.folder=system.file("regiondata", "ili20161121", package = "argo"),
#'   population.file=system.file("regiondata", "Population.csv", package = "argo"),
#'   gft.file=system.file("regiondata", "GFT.txt", package = "argo")
#' )
#' }
#'
#' @export
argo2_main <- function(gt.folder, ili.folder, population.file, gft.file, save.folder=NULL){
  #### parse each individual data file ####
  # library(xts)
  # library(glmnet)
  # library(argo)
  # library(parallel)
  # library(boot)
  if(!is.null(save.folder)) dir.create(save.folder, showWarnings = FALSE, recursive = TRUE)


  reg_data <- load_reg_data(gt.folder=gt.folder,
                            ili.folder=ili.folder,
                            population.file=population.file,
                            gft.file=gft.file)

  ili_national <- reg_data$ili_national
  ili_regional <- reg_data$ili_regional
  GT_national <- reg_data$GT_national
  GT_regional <- reg_data$GT_regional

  #### first-step argo prediction ####
  transY <- function(y){
    logit((y+1e-6) / 100)
  }

  inv_transY <- function(y){
    100*logit_inv(y)-1e-6
  }


  get_argo1 <- function(terms, period){
    common_idx <- period
    set.seed(1000)

    # national
    argo_national <- argo(data = transY(ili_national[common_idx]),
                          exogen = log(GT_national[common_idx, terms]+1),
                          mc.cores = parallel::detectCores())

    # regional
    argo_result <- list()
    all.pred <- list()
    for(region.id in 1:10){
      j <- paste0("Region.", region.id)
      set.seed(1000)
      argo_result[[j]] <- argo(transY(ili_regional[common_idx, j]),
                               log(GT_regional[[region.id]][common_idx, terms]+1),
                               mc.cores = parallel::detectCores(),
                               N_lag = NULL)

      pred_xts_blend <- inv_transY(argo_result[[j]]$pred)
      pred_xts_blend <- merge(ili_regional[,j], pred_xts_blend, all=FALSE)
      pred_xts_blend$naive <- c(NA, as.numeric(pred_xts_blend[1:(nrow(pred_xts_blend)-1), j]))
      names(pred_xts_blend)[1] <- "CDC.data"
      all.pred[[j]] <- pred_xts_blend
      print(j)
    }
    list(argo_national=argo_national,
         argo_result=argo_result,
         all.pred=all.pred)
  }

  common_idx <- index(merge(ili_national, GT_national, all=FALSE))
  common_idx09 <- common_idx[common_idx < as.Date("2010-05-22")]
  terms <- colnames(GT_national)
  terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))

  # first-step argo for two periods
  argo1.09 <- get_argo1(terms09, common_idx09)
  argo1.post10 <- get_argo1(terms, common_idx)

  #### blend 09 data to all data ####
  argo1.coef <- list()
  all.pred <- argo1.post10$all.pred

  for(region.id in 1:10){
    j <- paste0("Region.", region.id)
    all.pred[[j]][index(argo1.09$all.pred[[j]]),"predict"] <- argo1.09$all.pred[[j]][,"predict"]
    argo1.coef[[j]] <- argo1.post10$argo_result[[j]]$coef
    argo1.coef[[j]][,colnames(argo1.09$argo_result[[j]]$coef)] <- NA
    argo1.coef[[j]][rownames(argo1.09$argo_result[[j]]$coef),colnames(argo1.09$argo_result[[j]]$coef)] <-
      argo1.09$argo_result[[j]]$coef
  }
  argo.nat.p <- inv_transY(argo1.post10$argo_national$pred)
  argo.nat.p[index(argo1.09$argo_national$pred)] <- inv_transY(argo1.09$argo_national$pred)
  argo.nat.coef <- argo1.post10$argo_national$coef
  argo.nat.coef[,colnames(argo1.09$argo_national$coef)] <- NA
  argo.nat.coef[rownames(colnames(argo1.09$argo_national$coef)),colnames(argo1.09$argo_national$coef)] <-
    argo1.09$argo_national$coef

  #### argo second step ####
  argo.pred <- lapply(all.pred, function(x) x[,"predict"])
  argo.pred <- do.call(merge, argo.pred)
  colnames(argo.pred) <- paste0(colnames(argo.pred)[1], ".", 1:10)

  ili_regional2003post <- ili_regional["2003/"]
  argo2_result <- argo2(ili_regional2003post, argo.pred, argo.nat.p)

  var.est <- sapply(1:dim(argo2_result$projection.mat)[3], function(it){
    stopifnot(argo2_result$projection.mat[,,it] ==
                argo2_result$sigma_zwzw.structured[1:10,-(1:10),it]%*%
                solve(argo2_result$sigma_zwzw.structured[-(1:10),-(1:10),it]+diag(diag(argo2_result$sigma_zwzw.empirical[-(1:10),-(1:10),it]))))

    argo2_result$sigma_zwzw.structured[1:10,1:10,it] - 0.5 *(
      argo2_result$sigma_zwzw.structured[1:10,-(1:10),it] %*%
        solve(argo2_result$sigma_zwzw.structured[-(1:10),-(1:10),it]+diag(diag(argo2_result$sigma_zwzw.empirical[-(1:10),-(1:10),it]))) %*%
        argo2_result$sigma_zwzw.structured[-(1:10),1:10,it]
    )
  }, simplify = "array")
  dimnames(var.est)[[3]] <- as.character(index(na.omit(argo2_result$twostep)))
  err.realized <- na.omit(argo2_result$truth - argo2_result$twostep)

  if(!is.null(save.folder)) pdf(file.path(save.folder, "interval_coverage.pdf"))
  interval.coverage <- list()
  for(region.id in 1:10){
    if(!is.null(save.folder)){
      hist(data.matrix(err.realized)[,region.id]/sqrt(var.est[region.id,region.id,]), breaks = 40,
           main=paste("region", region.id),
           probability = TRUE)
      plot.function(dnorm, from=-10, to=10, n=1001, col=2, add=TRUE)
      mtext("(truth - estimate) / sqrt(forecasted variance)")
    }

    interval.coverage[[region.id]] <-
      c(mean(abs(data.matrix(err.realized)[,region.id]/sqrt(var.est[region.id,region.id,])) <= 1.96),
        var(data.matrix(err.realized)[,region.id]/sqrt(var.est[region.id,region.id,])))
  }
  if(!is.null(save.folder)) dev.off()
  interval.coverage <- do.call(rbind, interval.coverage)
  rownames(interval.coverage) <- colnames(err.realized)
  print(xtable::xtable(data.frame(interval.coverage), digits = 3))

  #### benchmark models ####
  gft.pred <- reg_data$GFT

  naive.pred <- lapply(all.pred, function(x) x[,"naive"])
  naive.pred <- do.call(merge, naive.pred)
  colnames(naive.pred) <- paste0(colnames(naive.pred)[1], ".", 1:10)


  common_idx <- index(na.omit(merge(ili_regional2003post, naive.pred)))
  var.pred <- mclapply(1:10, function(region.id){
    var.each <- argo(data = transY(ili_regional2003post)[common_idx, region.id],
                     exogen = transY(naive.pred)[common_idx],
                     N_lag = NULL, N_training = 104, alpha = NA, mc.cores = 1)
    inv_transY(var.each$pred)
  }, mc.cores = parallel::detectCores())
  var.pred <- do.call(merge, var.pred)
  colnames(var.pred) <- paste0("var1.Region.",1:10)

  common_idx <- index(na.omit(merge(ili_regional2003post, gft.pred, naive.pred)))
  var1gft.pred <- mclapply(1:10, function(region.id){
    var1gft.each <- argo(data = transY(ili_regional2003post)[common_idx, region.id],
                         exogen = merge(transY(naive.pred), transY(gft.pred+1))[common_idx],
                         N_lag = NULL, N_training = 104, alpha = NA, mc.cores = 1)
    inv_transY(var1gft.each$pred)
  }, mc.cores = parallel::detectCores())
  var1gft.pred <- do.call(merge, var1gft.pred)
  colnames(var1gft.pred) <- paste0("var1gft.Region.",1:10)


  #### aggreage predicitons ####
  reg.method.p <- sapply(1:10, function(region.id){
    pred.region <-
      merge(ili_regional2003post[,region.id],
            argo2_result$onestep[,region.id],
            argo2_result$twostep[,region.id],
            var.pred[,region.id],
            gft.pred[,region.id],
            var1gft.pred[,region.id],
            lag(ili_regional2003post,1)[,region.id])
    colnames(pred.region) <- c("CDC.data", paste0("argo", 1:2),
                               "var1","GFT","var1gft","naive")
    pred.region <- pred.region[is.finite(pred.region$argo2),]
    data.matrix(pred.region["2009-03-29/2016-10-01"])
  }, simplify = "array")


  #### comparisons ####

  zoom_periods <- c(#"2011-01-22/2014-09-27",
    "2009-03-29/2009-12-27",
    "2010-10-03/2011-05-22",
    "2011-10-02/2012-05-20",
    "2012-09-30/2013-05-19",
    "2013-09-29/2014-05-18",
    "2014-09-28/2015-05-24",
    "2015-10-10/2016-05-22",
    "2009-03-29/2015-08-15")
  eval.period <- "2009-03-29/2016-10-01" #"2009-01-10/2014-09-27"


  tab.allregion <- mclapply(1:10, function(region.id){
    tab <- summary_argo(xts(reg.method.p[,,region.id], as.Date(dimnames(reg.method.p)[[1]])),
                        dimnames(reg.method.p)[[2]][-1], dimnames(reg.method.p)[[2]][-1],
                        zoom_periods, eval.period)
  }, mc.cores = parallel::detectCores())
  tab.allregion <- sapply(c("mse","abse","rmspe","mape","corr", "corr_diff",
                            "bias", "percent_up", "mse_up", "mse_down"), function(type){
                              sapply(1:10, function(region.id){
                                tab.allregion[[region.id]][[type]]
                              }, simplify = "array")
                            }, simplify = "array")
  dimnames(tab.allregion)[[3]] <- paste0("Region.",1:10)
  dim(tab.allregion)

  tab.allregion[,"2009-03-29/2016-10-01",,"bias"]
  tab.allregion[,"2009-03-29/2016-10-01",,"percent_up"]
  tab.allregion[,"2009-03-29/2016-10-01",,"mse_up"]
  tab.allregion[,"2009-03-29/2016-10-01",,"mse_down"]

  testthat::test_that("mse, mse_up, mse_down consistent", {
    x1p1 <- tab.allregion[,,,"mse_up"] * tab.allregion[,,,"percent_up"]
    x1p1[is.na(x1p1)] <- 0
    x1p2 <- tab.allregion[,,,"mse_down"] * (1-tab.allregion[,,,"percent_up"])
    x1p2[is.na(x1p2)] <- 0
    x1 <- x1p1 + x1p2
    x1[x1==0] <- NA
    testthat::expect_equal(x1, tab.allregion[,,,"mse"])
  })
  which(tab.allregion[,,,"percent_up"]==1)

  which(abs(tab.allregion[,,,"mse_up"] * tab.allregion[,,,"percent_up"] +
              tab.allregion[,,,"mse_down"] * (1-tab.allregion[,,,"percent_up"])
            - tab.allregion[,,,"mse"]) > 1e-9)

  print("mse")
  rbind(t(tab.allregion[,"2009-03-29/2016-10-01",,"mse"]),
        colMeans(t(tab.allregion[,"2009-03-29/2016-10-01",,"mse"])))

  rbind(t(tab.allregion[,"2009-03-29/2016-10-01",,"mape"]),
        colMeans(t(tab.allregion[,"2009-03-29/2016-10-01",,"mape"])))

  rbind(t(tab.allregion[,"2009-03-29/2015-08-15",,"mse"]),
        colMeans(t(tab.allregion[,"2009-03-29/2015-08-15",,"mse"])))

  #### bootstrap confidence interval ####
  # overall
  set.seed(1000)
  boot.id <- tsboot(1:dim(reg.method.p)[[1]], identity, R = 1e3, l = 52, sim="geom")

  boot.err <- mclapply(1:boot.id$R, function(b.it) {
    reg.method.p.b <- reg.method.p[boot.id$t[b.it,],,]

    tab.allregion.b <- sapply(1:10, function(region.id){
      tab.b <- summary_argo(xts(reg.method.p.b[,,region.id], as.Date(dimnames(reg.method.p)[[1]])),
                            dimnames(reg.method.p)[[2]][-1], dimnames(reg.method.p)[[2]][-1],
                            NULL, "2009-01-17/2016-10-22")
      tab.b <- do.call(cbind, tab.b[c("rmse","abse","mape")])
      tab.b[,1] <- tab.b[,1]^2
      colnames(tab.b) <- c("mse","abse","mape")
      tab.b[c("argo2","var1","naive"),]
    }, simplify = "array")
    apply(tab.allregion.b,c(1,2),mean)
  }, mc.cores = parallel::detectCores())
  boot.err <- sapply(boot.err, identity, simplify = "array")

  boot.t0 <- sapply(1:10, function(region.id){
    tab.b <- summary_argo(xts(reg.method.p[,,region.id], as.Date(dimnames(reg.method.p)[[1]])),
                          dimnames(reg.method.p)[[2]][-1], dimnames(reg.method.p)[[2]][-1],
                          NULL, "2009-01-17/2016-10-22")
    tab.b <- do.call(cbind, tab.b[c("rmse","abse","mape")])
    tab.b[,1] <- tab.b[,1]^2
    colnames(tab.b) <- c("mse","abse","mape")
    tab.b[c("argo2","var1","naive"),]
  }, simplify = "array")
  boot.t0 <- apply(boot.t0, 1:2, mean)

  effi.t0 <- log(boot.t0[-1,1]) - log(boot.t0[1,1])
  effi.b <- t(t(log(boot.err[-1,1,])) - log(boot.err[1,1,]))
  ci <- apply(2*effi.b - effi.t0, 1, function(x) quantile(x, c(0.025, 0.975)))
  tab.ci <- format(round(exp(rbind(effi.t0, ci)),2), nsmall=2)
  tab.ci <- cbind(tab.ci[1,], paste0("[",tab.ci[2,],", ",tab.ci[3,],"]"))
  colnames(tab.ci) <- c("point estimate", "95% CI")
  # up to 2015-08-15 for GFT-based metrics
  set.seed(1000)
  reg.method.p.gft <- reg.method.p[1:333,,]
  boot.gft.id <- tsboot(1:dim(reg.method.p.gft)[[1]], identity, R = 1e3, l = 52, sim="geom")

  boot.gft.err <- mclapply(1:boot.id$R, function(b.it) {
    reg.method.p.b <- reg.method.p.gft[boot.gft.id$t[b.it,],,]

    tab.allregion.b <- sapply(1:10, function(region.id){
      tab.b <- summary_argo(xts(reg.method.p.b[,,region.id], as.Date(dimnames(reg.method.p.gft)[[1]])),
                            dimnames(reg.method.p.gft)[[2]][-1], dimnames(reg.method.p.gft)[[2]][-1],
                            NULL, "2009-01-17/2016-10-22")
      tab.b <- do.call(cbind, tab.b[c("rmse","abse","mape")])
      tab.b[,1] <- tab.b[,1]^2
      colnames(tab.b) <- c("mse","abse","mape")
      tab.b[c("argo2","GFT","var1gft","var1","naive"),]
    }, simplify = "array")
    apply(tab.allregion.b,c(1,2),mean)
  }, mc.cores = parallel::detectCores())
  boot.gft.err <- sapply(boot.gft.err, identity, simplify = "array")

  boot.gft.t0 <- sapply(1:10, function(region.id){
    tab.b <- summary_argo(xts(reg.method.p.gft[,,region.id], as.Date(dimnames(reg.method.p.gft)[[1]])),
                          dimnames(reg.method.p.gft)[[2]][-1], dimnames(reg.method.p.gft)[[2]][-1],
                          NULL, "2009-01-17/2016-10-22")
    tab.b <- do.call(cbind, tab.b[c("rmse","abse","mape")])
    tab.b[,1] <- tab.b[,1]^2
    colnames(tab.b) <- c("mse","abse","mape")
    tab.b[c("argo2","GFT","var1gft","var1","naive"),]
  }, simplify = "array")
  boot.gft.t0 <- apply(boot.gft.t0, 1:2, mean)

  effi.gft.t0 <- log(boot.gft.t0[-1,1]) - log(boot.gft.t0[1,1])
  effi.gft.b <- t(t(log(boot.gft.err[-1,1,])) - log(boot.gft.err[1,1,]))
  ci.gft <- apply(2*effi.gft.b - effi.gft.t0, 1, function(x) quantile(x, c(0.025, 0.975)))
  tab.ci.gft <- format(round(exp(rbind(effi.gft.t0, ci.gft)),2), nsmall=2)
  tab.ci.gft <- cbind(tab.ci.gft[1,], paste0("[",tab.ci.gft[2,],", ",tab.ci.gft[3,],"]"))
  colnames(tab.ci.gft) <- c("point estimate", "95% CI")
  tab.ci.all <- rbind(tab.ci.gft[1:2,], tab.ci)
  print(xtable::xtable(tab.ci.all))
  #### plotting ####
  if(!is.null(save.folder)){
    pdf(paste0(save.folder,"heatmap_allregs.pdf"),width=4, height=9)
    nterm <- ncol(GT_national)+1
    for(region.id in 1:10){
      j <- paste0("Region.", region.id)
      coef.mat <- argo1.coef[[region.id]]
      num.lag <- length(argo1.post10$argo_result[[j]]$parm$N_lag)
      coef.mat[(num.lag+2):(num.lag+nterm),] <- coef.mat[((num.lag+2):(num.lag+nterm))[order(rowSums(is.na(coef.mat[(num.lag+2):(num.lag+nterm),])))],]
      argo::heatmap_argo(coef.mat, 0.1, na.grey = TRUE, scale = 1)
      title(j)
    }
    dev.off()



    pdf(file.path(save.folder, "heatmap_2nd_step_v2.pdf"), width = 4,height = 7)
    heatmap_argo(argo2_result$mean.mat[1:10,], lim=0.5, na.grey = FALSE, scale=1)
    title("mean of truth - naive.pred")
    heatmap_argo(argo2_result$mean.mat[11:20,], lim=0.5, na.grey = FALSE, scale=1)
    title("mean of argo.pred - naive.pred")
    heatmap_argo(argo2_result$mean.mat[21:30,], lim=0.5, na.grey = FALSE, scale=1)
    title("mean of argo.nat.p - naive.pred")
    heatmap_argo(argo2_result$mean.mat[31:40,], lim=0.5, na.grey = FALSE, scale=1)
    title("mean of naive.pred - naive2.pred")

    for(j in 1:10){
      x <- argo2_result$projection.mat[j,,]
      rownames(x) <- paste0(rep(c("GT","ARGO.nat","momentum"),each=10),".Reg.",1:10)
      # print(range(x))
      heatmap_argo(x, lim=0.1, na.grey = FALSE, scale=1)
      title(paste0("Region.",j))
    }
    dev.off()


    heat.vec <- argo2_result$heat.vec
    heat <- lapply(104:nrow(heat.vec), function(t) cov(heat.vec[(t-103):t,]))
    avg.heat <- Reduce("+", heat)/length(heat)
    colnames(avg.heat) <- rownames(avg.heat) <- paste0(rep(c("detla.ili.","past.detla.ili.", "err.argo.", "err.nat."), each=10),
                                                       rep(1:10,3))

    heat.vec.structured <- apply(argo2_result$heat.vec.structured,1:2,mean)
    colnames(heat.vec.structured) <- rownames(heat.vec.structured) <- colnames(avg.heat)
    pdf(file.path(save.folder, "heatmap_cor.pdf"), width = 6,height = 6)
    avg.heat.cor <- diag(1/sqrt(diag(avg.heat))) %*% avg.heat %*% diag(1/sqrt(diag(avg.heat)))
    heatmap_cor(avg.heat)
    heatmap_cor(avg.heat.cor)
    title("Heat map for model justification (empirical)")
    abline(v=1:3*10+0.5, h=1:3*10+0.5)
    heat.vec.structured.cor <- diag(1/sqrt(diag(heat.vec.structured))) %*% heat.vec.structured %*% diag(1/sqrt(diag(heat.vec.structured)))
    heatmap_cor(heat.vec.structured)
    heatmap_cor(heat.vec.structured.cor)
    title("Heat map for model justification (structured)")
    abline(v=1:3*10+0.5, h=1:3*10+0.5)
    dev.off()


    avg.sigma_ww.structured <- apply(argo2_result$sigma_ww.structured,1:2,mean)
    avg.sigma_ww.empirical <- apply(argo2_result$sigma_ww.empirical,1:2,mean)

    avg.sigma_zw.structured <- apply(argo2_result$sigma_zw.structured,1:2,mean)
    avg.sigma_zw.empirical <- apply(argo2_result$sigma_zw.empirical,1:2,mean)

    avg.sigma_zwzw.structured <- apply(argo2_result$sigma_zwzw.structured,1:2,mean)
    avg.sigma_zwzw.empirical <- apply(argo2_result$sigma_zwzw.empirical,1:2,mean)
    rownames(avg.sigma_zwzw.structured) <- colnames(avg.sigma_zwzw.structured) <-
      rownames(avg.sigma_zwzw.empirical) <- colnames(avg.sigma_zwzw.empirical) <- rep(" ",40)

    cor_zwzw.structured <-
      apply(argo2_result$sigma_zwzw.structured, 3, function(x) diag(1/sqrt(diag(x)))%*%x%*%diag(1/sqrt(diag(x))))
    cor_zwzw.empirical <-
      apply(argo2_result$sigma_zwzw.empirical, 3, function(x) diag(1/sqrt(diag(x)))%*%x%*%diag(1/sqrt(diag(x))))
    cor_zwzw.structured <- sapply(data.frame(cor_zwzw.structured), function(x) matrix(x, nrow=sqrt(length(x))),
                                  simplify = "array")
    cor_zwzw.empirical <- sapply(data.frame(cor_zwzw.empirical), function(x) matrix(x, nrow=sqrt(length(x))),
                                 simplify = "array")


    avg.cor_zwzw.structured <- apply(cor_zwzw.structured,1:2,mean)
    avg.cor_zwzw.empirical <- apply(cor_zwzw.empirical,1:2,mean)
    rownames(avg.cor_zwzw.structured) <- colnames(avg.cor_zwzw.structured) <-
      rownames(avg.cor_zwzw.empirical) <- colnames(avg.cor_zwzw.empirical) <- rep(" ",40)


    pdf(file.path(save.folder, "sigma.pdf"), width = 6,height = 6)
    heatmap_cor(avg.sigma_ww.structured)
    title("Average Sigma_WW(structured)")
    heatmap_cor(avg.sigma_ww.empirical)
    title("Average Sigma_WW(empirical)")

    heatmap_cor(avg.sigma_zw.structured)
    title("Average Sigma_ZW(structured)")
    heatmap_cor(avg.sigma_zw.empirical)
    title("Average Sigma_ZW(empirical)")

    heatmap_cor(rbind(avg.sigma_zw.structured, avg.sigma_ww.structured))
    title("Average Cov( (Z W), W) (structured)")
    abline(h=c(10,20)+0.5, v=c(10,20)+0.5, col="grey50")
    abline(h=30.5, lwd=3)
    heatmap_cor(rbind(avg.sigma_zw.empirical, avg.sigma_ww.empirical))
    title("Average Cov( (Z W), W) (empirical)")
    abline(h=c(10,20)+0.5, v=c(10,20)+0.5, col="grey50")
    abline(h=30.5, lwd=3)

    heatmap_cor(avg.sigma_zwzw.structured)
    title("structured covariance")
    abline(h=c(10,20)+0.5, v=c(20,30)+0.5, col=1)
    abline(h=30.5, v=10.5, lwd=1)
    heatmap_cor(avg.cor_zwzw.structured)
    title("structured correlation")
    abline(h=c(10,20)+0.5, v=c(20,30)+0.5, col=1)
    abline(h=30.5, v=10.5, lwd=1)
    heatmap_cor(avg.sigma_zwzw.empirical)
    title("empirical covariance")
    abline(h=c(10,20)+0.5, v=c(20,30)+0.5, col=1)
    abline(h=30.5, v=10.5, lwd=1)
    heatmap_cor(avg.cor_zwzw.empirical)
    title("empirical correlation")
    abline(h=c(10,20)+0.5, v=c(20,30)+0.5, col=1)
    abline(h=30.5, v=10.5, lwd=1)
    dev.off()
  }
}

#' main function for argox
#'
#' Main function that reproduce the results in ARGOX paper. The datasets are available at Harvard Dataverse <doi:10.7910/DVN/2IVDGK>.
#'
#' @param gt.folder folder with Google Trends files, which should be thousands of csv file such as "US-MA_fever cough.csv" or "US-NY_cold or flu.csv"
#' @param ili.folder folder with ILINet data files: "ILINet_nat.csv" and "ILINet_regional.csv"
#' @param population.file file path to population csv file
#' @param gft.file file path to Google Flu Trends csv file
#' @param mix the weighted avarage mixing of raw state-level Google Trends data. Set to be 0 for stand-alone model. Set to be 1/3 for spatial-pooling model.
#' @param save.folder output folder to save graphics. If NULL then do not output graphics.
#' @param NCORES number of parallel cpu cores to be used.
#'
#' @references Yang, S., Ning, S. & Kou, S.C. Use Internet search data to accurately track state level influenza epidemics. Sci Rep 11, 4023 (2021)
#' @import xts
#' @import glmnet
#' @import parallel
#' @import boot
#'
#'
#' @export
argox_main <- function(gt.folder, ili.folder, population.file, gft.file, mix, save.folder=NULL, NCORES = 8){
  reg_data <- load_reg_data(gt.folder=gt.folder,
                            ili.folder=ili.folder,
                            population.file=population.file,
                            gft.file=gft.file,
                            gt.parser = gt.parser.pub.api)
  names(reg_data$GT_state)[names(reg_data$GT_state)=="501"] <- "US.NY.NYC"
  names(reg_data$ili_state)[names(reg_data$ili_state)=="US.New York City"] <- "US.NY.NYC"
  summary(t(sapply(reg_data$GT_state, dim)))
  summary(t(sapply(reg_data$GT_state_filled, dim)))
  
  ili_national <- reg_data$ili_national
  ili_regional <- reg_data$ili_regional
  ili_state <- reg_data$ili_state
  GT_national <- reg_data$GT_national
  GT_regional <- reg_data$GT_regional
  GT_state <- reg_data$GT_state
  
  state_names <- setdiff(names(reg_data$GT_state), c("US", "US.FL"))
  state_info <- read.csv(population.file)
  GT_regstate_mixed <- list()
  for (each_state in state_names){
    region_id_for_state = state_info[state_info$Abbre==strsplit(each_state, "\\.")[[1]][2], "Region"]
    GT_regstate_mixed[[each_state]] = (reg_data$GT_state[[each_state]] * (1-mix) + GT_regional[[region_id_for_state]] * mix)
  }
  
  
  #### first-step argo prediction ####
  transY <- function(y){
    logit((y+0.1) / 100)
  }
  
  inv_transY <- function(y){
    100*logit_inv(y)-0.1
  }
  
  
  get_argo1 <- function(terms, period_nat, period_reg){
    common_idx <- period_nat
    set.seed(1000)
    
    # national 
    GT_national <- merge(GT_national, ili_national)
    GT_national$ili_national <- NULL
    common_idx_nat_append <- c(common_idx[1] - (52:1)*7, common_idx)
    argo_national <- argo(data = transY(ili_national[common_idx_nat_append]),
                          exogen = log(GT_national[common_idx_nat_append, terms]+1),
                          mc.cores = NCORES)
    
    # regional
    argo_regional <- list()
    pred_regional <- list()
    for(region.id in 1:10){
      j <- paste0("Region.", region.id)
      set.seed(1000)
      argo_regional[[j]] <- argo(transY(ili_regional[common_idx, j]),
                                 log(GT_regional[[region.id]][common_idx, terms]+1),
                                 mc.cores = NCORES,
                                 N_lag = NULL)
      
      pred_xts_blend <- inv_transY(argo_regional[[j]]$pred)
      pred_xts_blend <- merge(ili_regional[,j], pred_xts_blend, all=FALSE)
      pred_xts_blend$naive <- c(NA, as.numeric(pred_xts_blend[1:(nrow(pred_xts_blend)-1), j]))
      names(pred_xts_blend)[1] <- "CDC.data"
      pred_regional[[j]] <- pred_xts_blend
      print(j)
    }
    
    # state-level
    common_idx <- period_reg
    argo_state <- list()
    pred_state <- list()
    for(region.id in 1:length(state_names)){
      j <- state_names[region.id]
      set.seed(1000)
      argo_state[[j]] <- argo(transY(ili_state[common_idx, j]),
                              log(GT_regstate_mixed[[j]][common_idx, terms]+1),
                              mc.cores = NCORES,
                              N_lag = NULL)
      
      pred_xts_blend <- inv_transY(argo_state[[j]]$pred)
      pred_xts_blend <- merge(ili_state[,j], pred_xts_blend, all=FALSE)
      pred_xts_blend$naive <- c(NA, as.numeric(pred_xts_blend[1:(nrow(pred_xts_blend)-1), j]))
      names(pred_xts_blend)[1] <- "CDC.data"
      pred_state[[j]] <- pred_xts_blend
      print(j)
    }
    list(argo_national=argo_national,
         argo_regional=argo_regional,
         pred_regional=pred_regional,
         argo_state=argo_state,
         pred_state=pred_state)
  }
  
  common_idx_nat <- index(merge(ili_national, GT_national, all=FALSE))
  common_idx_reg <- index(merge(ili_state, GT_national, all=FALSE))
  terms <- colnames(GT_national)
  
  # first-step argo for two periods
  argo1.post10 <- get_argo1(terms, common_idx_nat, common_idx_reg)
  
  argo1.coef <- list()
  pred_state <- argo1.post10$pred_state
  
  for(j in state_names){
    argo1.coef[[j]] <- argo1.post10$argo_state[[j]]$coef
  }
  argo.nat.p <- inv_transY(argo1.post10$argo_national$pred)
  
  argo.nat.coef <- argo1.post10$argo_national$coef
  
  #### argo second step ####
  argo.reg.p <- lapply(argo1.post10$pred_regional, function(x) x[,"predict"])
  argo.reg.p <- do.call(merge, argo.reg.p)
  colnames(argo.reg.p) <- names(argo1.post10$pred_regional)
  
  argo.state.p <- lapply(pred_state, function(x) x[,"predict"])
  argo.state.p <- do.call(merge, argo.state.p)
  colnames(argo.state.p) <- names(pred_state)
  
  argo2 <- function(truth, argo.state.p, argo.reg.p, argo.nat.p, state_names, state_info, use_yt2=TRUE){
    naive.p <- truth
    index(naive.p) <- index(truth) + 7
    common_idx <- index(na.omit(merge(truth, naive.p, argo.state.p, argo.nat.p)))
    
    Y <- truth - naive.p
    Yt2 <- Y
    index(Yt2) <- index(Yt2) + 7
    
    argo.nat.p <- argo.nat.p[common_idx]
    argo.reg.p <- argo.reg.p[common_idx]
    naive.p <- naive.p[common_idx]
    truth <- truth[common_idx]
    argo.state.p <- argo.state.p[common_idx]
    
    X <- argo.state.p - naive.p
    
    X.nat <- as.numeric(argo.nat.p) - naive.p
    X.nat <- X.nat[common_idx]
    
    argo.reg.p.dup <- lapply(state_names, function(each_state){
      region_id_for_state = state_info[state_info$Abbre==strsplit(each_state, "\\.")[[1]][2], "Region"]
      argo.reg.p[,region_id_for_state]
    })
    argo.reg.p.dup <- do.call(merge, argo.reg.p.dup)
    names(argo.reg.p.dup) <- state_names
    X.reg <- argo.reg.p.dup - naive.p
    
    projection.mat <- list()
    mean.mat <- list()
    
    Y.pred <- X
    Y.pred[] <- NA
    
    zw_used <- list()
    sigma_ww.structured <- sigma_ww.empirical <-
      sigma_zw.structured <- sigma_zw.empirical <-
      heat.vec.structured <-
      sigma_zwzw.structured <- sigma_zwzw.empirical <- list()
    
    for(it in 105:length(common_idx)){
      training_idx <- common_idx[(it-104):(it-1)]
      t.now <- common_idx[it]
      y <- Y[training_idx,]
      x <- X[training_idx,]
      x.nat <- X.nat[training_idx,]
      x.reg <- X.reg[training_idx,]
      yt2 <- Yt2[training_idx,]
      
      sigma_yy <- var(y)
      
      m1 <- cor(y, yt2)
      m2 <- cor(y)
      rho.l2 <- sum(m1*m2)/sum(m2^2)
      
      autocov.y.yt2 <- rho.l2*sigma_yy
      
      vcov.x_xreg_xnat <- 
        cbind(rbind(sigma_yy+var((argo.state.p - truth)[training_idx,]),sigma_yy, sigma_yy),
              rbind(sigma_yy,sigma_yy+var((argo.reg.p.dup - truth)[training_idx,]), sigma_yy),
              rbind(sigma_yy,sigma_yy,sigma_yy+var((as.numeric(argo.nat.p) - truth)[training_idx,])))
      sigma_zw <- cbind(sigma_yy,sigma_yy,sigma_yy)
      
      if(use_yt2){
        vcov.x_xreg_xnat <- cbind(vcov.x_xreg_xnat, rbind(autocov.y.yt2,autocov.y.yt2,autocov.y.yt2))
        vcov.x_xreg_xnat <- rbind(vcov.x_xreg_xnat, cbind(t(autocov.y.yt2),t(autocov.y.yt2),t(autocov.y.yt2),sigma_yy))
        sigma_zw <- cbind(sigma_zw, autocov.y.yt2)
      }
      
      # not shrinked
      if(use_yt2){
        y.pred.blp <- colMeans(y) +
          sigma_zw %*% 
          solve(vcov.x_xreg_xnat, c(t(X[t.now,])-colMeans(x), t(X.reg[t.now,])-colMeans(x.reg), t(X.nat[t.now,])-colMeans(x.nat), t(Yt2[t.now,]-colMeans(yt2))))
      }else{
        y.pred.blp <- colMeans(y) +
          sigma_zw %*% 
          solve(vcov.x_xreg_xnat, c(t(X[t.now,])-colMeans(x), t(X.reg[t.now,])-colMeans(x.reg), t(X.nat[t.now,])-colMeans(x.nat)))
      }
      
      Kzz <- solve((1-rho.l2^2)*sigma_yy)
      Kgt <- solve(var((argo.state.p - truth)[training_idx,]))
      Kreg <- solve(var((argo.reg.p.dup - truth)[training_idx,]))
      Knat <- solve(var((as.numeric(argo.nat.p) - truth)[training_idx,]))
      
      if(use_yt2){
        y.pred.bayes <- colMeans(y) + 
          solve(Kzz+Knat+Kreg+Kgt,
                Knat%*%(t(X.nat[t.now,])-colMeans(x.nat)) + Kgt%*%(t(X[t.now,])-colMeans(x)) + 
                  Kreg%*%(t(X.reg[t.now,])-colMeans(x.reg)) + 
                  Kzz%*%(rho.l2*(t(Yt2[t.now,]-colMeans(yt2)))))
      }else{
        y.pred.bayes <- colMeans(y) + 
          solve(solve(sigma_yy) + Knat+Kreg+Kgt,  # prior of y is mean 0 vcov sigma_yy
                Knat%*%(t(X.nat[t.now,])-colMeans(x.nat)) + Kgt%*%(t(X[t.now,])-colMeans(x)) + 
                  Kreg%*%(t(X.reg[t.now,])-colMeans(x.reg)))
      }
      
      if(all(is.finite(y.pred.blp))){
        stopifnot(all(abs(y.pred.blp-y.pred.bayes) < 1e-8))
      }
      
      # shrinked  
      if(use_yt2){
        y.pred <- colMeans(y) +
          sigma_zw %*% 
          solve(vcov.x_xreg_xnat + diag(diag(var(cbind(X,X.reg,X.nat,Yt2)[training_idx,]))), 
                c(t(X[t.now,])-colMeans(x), t(X.reg[t.now,])-colMeans(x.reg), t(X.nat[t.now,])-colMeans(x.nat), t(Yt2[t.now,]-colMeans(yt2))))
      }else{
        y.pred <- colMeans(y) +
          sigma_zw %*% 
          solve(vcov.x_xreg_xnat + diag(diag(var(cbind(X,X.reg,X.nat)[training_idx,]))), 
                c(t(X[t.now,])-colMeans(x), t(X.reg[t.now,])-colMeans(x.reg), t(X.nat[t.now,])-colMeans(x.nat)))
      }
      
      Y.pred[t.now, ] <- t(y.pred)
      projection.mat[[as.character(t.now)]] <- sigma_zw %*% solve(vcov.x_xreg_xnat + diag(diag(var(cbind(X,X.reg,X.nat,Yt2)[training_idx,]))))
      mean.mat[[as.character(t.now)]] <- c(colMeans(y), colMeans(x), colMeans(x.reg), colMeans(x.nat), colMeans(yt2))
      
      sigma_ww <- vcov.x_xreg_xnat
      sigma_zz <- sigma_yy
      sigma_ww.structured[[as.character(t.now)]] <- sigma_ww
      sigma_ww.empirical[[as.character(t.now)]] <- var(cbind(X,X.reg,X.nat,Yt2)[training_idx,])
      sigma_zw.structured[[as.character(t.now)]] <- sigma_zw
      sigma_zw.empirical[[as.character(t.now)]] <- cov(Y[training_idx,], cbind(X,X.reg,X.nat,Yt2)[training_idx,])
      
      sigma_zwzw.structured[[as.character(t.now)]] <- rbind(
        cbind(sigma_zz, sigma_zw),
        cbind(t(sigma_zw), sigma_ww)
      )
      sigma_zwzw.empirical[[as.character(t.now)]] <- var(cbind(Y,X,X.reg,X.nat,Yt2)[training_idx,])
      zw_used[[as.character(t.now)]] <- cbind(Y,X,X.reg,X.nat,Yt2)[training_idx,]
    }
    
    projection.mat <- sapply(projection.mat, identity, simplify = "array")
    mean.mat <- sapply(mean.mat, identity, simplify = "array")
    
    argo2.p <- Y.pred + naive.p
    
    err.twostep <- argo2.p - truth
    
    
    heat.vec <- na.omit(merge(truth-naive.p, argo.state.p - truth, argo.reg.p.dup - truth, as.numeric(argo.nat.p) - truth, Yt2))
    colnames(heat.vec) <- paste0(rep(c("CDC.increment.", "err.argo.", "err.reg.", "err.nat.", "err.y2."), each=length(state_names)),
                                 state_names)
    
    sigma_ww.structured <- sapply(sigma_ww.structured, identity, simplify = "array")
    sigma_ww.empirical <- sapply(sigma_ww.empirical, identity, simplify = "array")
    sigma_zw.structured <- sapply(sigma_zw.structured, identity, simplify = "array")
    sigma_zw.empirical <- sapply(sigma_zw.empirical, identity, simplify = "array")
    sigma_zwzw.structured <- sapply(sigma_zwzw.structured, identity, simplify = "array")
    sigma_zwzw.empirical <- sapply(sigma_zwzw.empirical, identity, simplify = "array")
    zw_used <- sapply(zw_used, identity, simplify = "array")
    
    list(onestep=argo.state.p, twostep=argo2.p, naive=naive.p, truth=truth,
         Y.pred=Y.pred, err.twostep=err.twostep,
         heat.vec=heat.vec, projection.mat=projection.mat, mean.mat=mean.mat,
         sigma_ww.structured=sigma_ww.structured, sigma_ww.empirical=sigma_ww.empirical,
         sigma_zw.structured=sigma_zw.structured, sigma_zw.empirical=sigma_zw.empirical,
         sigma_zwzw.structured=sigma_zwzw.structured, sigma_zwzw.empirical=sigma_zwzw.empirical,
         zw_used=zw_used, zw_overall = cbind(Y,X,X.reg,X.nat,Yt2))
  }
  
  argo_ind <- function(truth, argo.state.p, argo.reg.p, argo.nat.p, state_names, state_info, use_reg=FALSE){
    naive.p <- truth
    index(naive.p) <- index(truth) + 7
    common_idx <- index(na.omit(merge(truth, naive.p, argo.state.p, argo.nat.p)))
    
    Y <- truth - naive.p
    Yt2 <- Y
    index(Yt2) <- index(Yt2) + 7
    
    argo.nat.p <- argo.nat.p[common_idx]
    argo.reg.p <- argo.reg.p[common_idx]
    naive.p <- naive.p[common_idx]
    truth <- truth[common_idx]
    argo.state.p <- argo.state.p[common_idx]
    
    X <- argo.state.p - naive.p
    
    X.nat <- as.numeric(argo.nat.p) - naive.p
    X.nat <- X.nat[common_idx]
    
    argo.reg.p.dup <- lapply(state_names, function(each_state){
      region_id_for_state = state_info[state_info$Abbre==strsplit(each_state, "\\.")[[1]][2], "Region"]
      argo.reg.p[,region_id_for_state]
    })
    argo.reg.p.dup <- do.call(merge, argo.reg.p.dup)
    names(argo.reg.p.dup) <- state_names
    X.reg <- argo.reg.p.dup - naive.p
    
    Y.pred <- X
    Y.pred[] <- NA
    
    for(it in 105:length(common_idx)){
      for(state_id in 1:ncol(argo.state.p)){
        training_idx <- common_idx[(it-104):(it-1)]
        t.now <- common_idx[it]
        y <- Y[training_idx,state_id]
        x <- X[training_idx,state_id]
        x.nat <- X.nat[training_idx,state_id]
        x.reg <- X.reg[training_idx,state_id]
        yt2 <- Yt2[training_idx,state_id]
        
        if(use_reg){
          sigma_zw <- cov(y, cbind(x, x.reg, x.nat, yt2))
          vcov.x_xreg_xnats <- var(cbind(x, x.reg, x.nat, yt2))
          
          y.pred <- colMeans(y) +
            sigma_zw %*% 
            solve(vcov.x_xreg_xnats + diag(diag(vcov.x_xreg_xnats)), 
                  c(X[t.now,state_id]-mean(x), X.reg[t.now,state_id]-mean(x.reg), X.nat[t.now,state_id]-mean(x.nat), Yt2[t.now,state_id]-mean(yt2)))
        }else{
          sigma_zw <- cov(y, cbind(x, x.nat, yt2))
          vcov.x_xreg_xnats <- var(cbind(x, x.nat, yt2))
          
          y.pred <- colMeans(y) +
            sigma_zw %*% 
            solve(vcov.x_xreg_xnats + diag(diag(vcov.x_xreg_xnats)), 
                  c(X[t.now,state_id]-mean(x), X.nat[t.now,state_id]-mean(x.nat), Yt2[t.now,state_id]-mean(yt2)))
          
        }
        Y.pred[t.now, state_id] <- y.pred
      }
    }
    
    argo2.p <- Y.pred + naive.p
    err.twostep <- argo2.p - truth
    
    list(onestep=argo.state.p, twostep=argo2.p, naive=naive.p, truth=truth,
         Y.pred=Y.pred, err.twostep=err.twostep)
  }
  
  argo2_result <- argo2(ili_state[,state_names], argo.state.p, argo.reg.p[index(argo.state.p)], argo.nat.p[index(argo.state.p)], state_names, state_info, TRUE)
  state_names_sub <- setdiff(state_names, c("US.HI", "US.AK"))
  state_names_sub <- setdiff(state_names_sub, c("US.ME", "US.MT", "US.ND", "US.SD", "US.VT"))
  argo2_result_sub <- argo2(ili_state[,state_names_sub], 
                            argo.state.p[,state_names_sub], 
                            argo.reg.p[index(argo.state.p)], 
                            argo.nat.p[index(argo.state.p)], 
                            state_names_sub, state_info, TRUE)
  argo_ind_result <- argo_ind(ili_state[,state_names], argo.state.p, argo.reg.p[index(argo.state.p)], argo.nat.p[index(argo.state.p)], state_names, state_info, use_reg = FALSE)
  
  if (!is.null(save.folder)){
    outDir <- paste0(save.folder, "/argo2_mix", round(mix, 2))
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
    save(list=ls(), file=paste0(outDir, "/argo2-state-GT-",tail(strsplit(gt.folder, "/")[[1]], 1),".rda"))
  }
  return(list(individual_pred=argo_ind_result$twostep, pooled_pred=argo2_result_sub$twostep, argo_ind_result=argo_ind_result, argo2_result_sub=argo2_result_sub))()
}
