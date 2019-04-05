#' Parsing of raw data
#'
#' Data related to the PNAS paper. Accessed on Nov 14, 2015.
#'
#' Parse and load CDC's ILI data, Google Flu Trend data,
#' Google Correlate data trained with ILI as of 2010, Google Correlate data trained with ILI as of 2009,
#' Google Trend data with search terms identified from Google Correlate (2010 version).
#'
#' Each week ends on the Saturday indicated in the xts object
#'
#' Google Correlate data is standardized by Google, and we rescale it to 0 -- 100 during parsing.
#' Google Trends data is in the scale of 0 -- 100.
#'
#' @importFrom zoo index
#'
#' @param type the type of the data to be loaded. If \code{type=="extdata"} it loads the data
#' to reproduce the PNAS paper, and if \code{type=="athdata"} it loads the data to reproduce the
#' CID(?) paper.
#' @param ili.weighted logical indicator to specify whether to load weighted ILI or not,
#' if \code{FALSE} unweighted ILI is loaded.
#'
#' @examples
#' system.file("extdata", "correlate-Influenza_like_Illness_h1n1_CDC_.csv", package = "argo")
#' system.file("extdata", "correlate-Influenza_like_Illness_CDC_.csv", package = "argo")
#' system.file("extdata", "GFT.csv", package = "argo")
#' system.file("extdata", "ILINet.csv", package = "argo")
#' load_data()
#'
#' @return A list of following named xts objects if \code{type=="extdata"}
#' \itemize{
#'  \item \code{GC10} Google Correlate trained with ILI available as of 2010.
#'    Available online at \url{https://www.google.com/trends/correlate/search?e=id:20xKcnNqHrk&t=weekly}
#'  \item \code{GC09} Google Correlate trained with ILI available as of 2009.
#'    Not directly available online, you have to manually input ILI time series
#'    at \url{https://www.google.com/trends/correlate}
#'  \item \code{GT} Google Trends data for search queries identified using Google Correlate.
#'    Not directly available online, you have to manually input query terms
#'    at \url{https://www.google.com/trends}
#'  \item \code{CDC} CDC's ILI dataset.
#'    Available online at \url{http://gis.cdc.gov/grasp/fluview/fluportaldashboard.html}
#'  \item \code{GFT} Google Flu Trend (historical predictions).
#'    Available online at \url{https://www.google.org/flutrends}
#' }
#'
#' A list of following named xts objects if \code{type=="athdata"}
#' \itemize{
#'  \item \code{GT} Google Trends data for search queries identified.
#'    Not directly available online, you have to manually input query terms
#'    at \url{https://www.google.com/trends}
#'  \item \code{CDC} CDC's ILI dataset.
#'    Available online at \url{http://gis.cdc.gov/grasp/fluview/fluportaldashboard.html}
#'  \item \code{ili_idx} the indexing information that includes the week number and year
#'  number, the date of ending Saturday, and the season number
#'    Available online at \url{http://www.cdc.gov/flu/weekly/}
#'  \item \code{ATH} Athenahealth data that includes the proportion of "Flu Visit",
#'  "ILI Visit", and "Unspecified Viral or ILI Visit" compared to total number of visit to the
#'  Athenahealth partner healthcare providers.
#'  \item \code{ili_unrevised} Historical unrevised ILI activity level.
#'    The unrevised ILI published on week ZZ of season XXXX-YYYY is available at
#'    \url{www.cdc.gov/flu/weekly/weeklyarchivesXXXX-YYYY/data/senAllregtZZ.html}
#'    or \url{.htm}. For example, original ILI report for week 7 of season 2015-2016 is available at
#'    \url{www.cdc.gov/flu/weekly/weeklyarchives2015-2016/data/senAllregt07.html},
#'    and original ILI report for week 50 of season 2012-2013 is available at
#'    \url{www.cdc.gov/flu/ weekly/weeklyarchives2012-2013/data/senAllregt50.htm}
#' }
#'
#' @references
#' Yang, S., Santillana, M., & Kou, S. C. (2015). Accurate estimation of influenza epidemics using Google search data via ARGO. Proceedings of the National Academy of Sciences, \href{https://dx.doi.org/10.1073/pnas.1515373112}{doi: 10.1073/pnas.1515373112}.
#'
#' @import utils
#'
#' @export
load_data <- function(type = "extdata", ili.weighted=TRUE) {
  if(type=="extdata"){
    resacle_gc <- function(GC_data) {
      GC_data$Date <- as.Date(as.character(GC_data$Date))
      GC_data$Date <- GC_data$Date + 6
      GFT_rescale <- xts::xts(GC_data[, -(1:2)], order.by = as.Date(GC_data$Date))
      for (i in 1:ncol(GFT_rescale)) {
        GFT_rescale[, i] <- GFT_rescale[, i] - min(GFT_rescale[, i])
        GFT_rescale[, i] <- GFT_rescale[, i] * 100/max(GFT_rescale[, i])
      }
      GFT_rescale
    }
    GC09_data <- read.csv(system.file(type, "correlate-Influenza_like_Illness_h1n1_CDC_.csv",
                                      package = "argo"), skip = 11)
    GC10_data <- read.csv(system.file(type, "correlate-Influenza_like_Illness_CDC_.csv",
                                      package = "argo"), skip = 11)

    GC09_data <- resacle_gc(GC09_data)
    GC10_data <- resacle_gc(GC10_data)
    updated_gft <- read.csv(system.file(type, "GFT.csv", package = "argo"),
                            skip = 11)
    updated_gft <- xts::xts(as.numeric(updated_gft$United.States)/1000,
                            order.by = as.Date(updated_gft$Date) + 6)
    names(updated_gft) <- "GFT"
    updated_gft <- na.omit(updated_gft)
  }else if(grepl("athdata", type)){
    ATH_data <- read.csv(system.file(type, "ATHpercent.csv", package = "argo"))
    ATH_data <- xts::xts(ATH_data[,-1], order.by = as.Date(ATH_data[,1]))
    if(ili.weighted){
      ili_unrevised <- read.csv(system.file(type, "ILI_unrevised_weighted.csv", package = "argo"))
    }else{
      ili_unrevised <- read.csv(system.file(type, "ILI_unrevised_unweighted.csv", package = "argo"))
    }
    ili_unrevised <- xts::xts(ili_unrevised[,-1], order.by = as.Date(ili_unrevised[,1]))
    colnames(ili_unrevised) <- as.Date(colnames(ili_unrevised), "X%Y.%m.%d")
  }


  GT_data <- read.csv(system.file(type, "GTdata.csv", package = "argo"))
  GT_data <- xts::xts(GT_data[, -1], order.by = as.Date(GT_data[,1]))

  ili <- read.csv(system.file(type, "ILINet.csv", package = "argo"),
                  skip = 1)

  ili_idx <- ili[,c("YEAR", "WEEK")]
  ili_idx$DATE <- as.Date("1997-10-04") + (1:nrow(ili)-1)*7
  ili_idx$SEASON <- ifelse(ili_idx$WEEK < 40, paste0(ili_idx$YEAR-1,"-",ili_idx$YEAR),
                           paste0(ili_idx$YEAR,"-",ili_idx$YEAR+1))

  ili <- xts::as.xts(ili[, c("X..WEIGHTED.ILI", "X.UNWEIGHTED.ILI")],
                     order.by = as.Date("1997-10-04") + (1:nrow(ili) - 1) * 7)
  if(ili.weighted){
    cdc <- ili[, "X..WEIGHTED.ILI"]
  }else{
    cdc <- ili[, "X.UNWEIGHTED.ILI"]
  }

  cdc[cdc=="X"] <- NA
  cdc <- xts::xts(as.numeric(cdc), order.by = index(cdc))
  names(cdc) <- "CDC.data"



  data_all <- list()
  data_all$GT <- GT_data
  data_all$CDC <- cdc
  data_all$ili_idx <- ili_idx
  if(type=="extdata"){
    data_all$GC09 <- GC09_data
    data_all$GC10 <- GC10_data
    data_all$GFT <- updated_gft
  }else if(grepl("athdata", type)){
    data_all$ATH <- ATH_data
    data_all$ili_unrevised <- ili_unrevised
  }
  data_all
}

#' Parsing of unrevised ili from online source
#' @param type the type of data folder to parse
#' @param ili.weighted indicator to use weighted ILI or not
#' @references
#' Yang, S., Santillana, M., & Kou, S. C. (2015). Accurate estimation of influenza epidemics using Google search data via ARGO. Proceedings of the National Academy of Sciences, \href{https://dx.doi.org/10.1073/pnas.1515373112}{doi: 10.1073/pnas.1515373112}.
#' @import xts XML
#' @examples
#'
#' \donttest{
#' parse_unrevised_ili()
#' }
#'
#' @export
parse_unrevised_ili <- function(type = "extdata", ili.weighted=TRUE){
  all_data <- load_data(type, ili.weighted)
  GFT_xts <- all_data$CDC
  GFT_xts <- GFT_xts["2004/"]
  ili_unrivised <- matrix(rep(GFT_xts$CDC.data, nrow(GFT_xts)), nrow=nrow(GFT_xts))




  colnames(ili_unrivised) <- as.character(index(GFT_xts))
  rownames(ili_unrivised) <- as.character(index(GFT_xts))

  for(i in 1:ncol(ili_unrivised)){
    ili_unrivised[i:nrow(ili_unrivised),i] <- NA
  }

  ili_idx <- all_data$ili_idx

  missing_data_id <- c()
  for(j in 2:ncol(ili_unrivised)){
    i <- which(ili_idx$DATE == index(GFT_xts)[j-1])
    if(ili_idx$DATE[i]>=as.Date("2014-10-04")){
      theurl <- paste0("http://www.cdc.gov/flu/weekly/weeklyarchives",
                       ili_idx$SEASON[i],"/data/senAllregt",
                       sprintf("%.2d", ili_idx$WEEK[i]),".html")
    }else{
      theurl <- paste0("http://www.cdc.gov/flu/weekly/weeklyarchives",
                       ili_idx$SEASON[i],"/data/senAllregt",
                       sprintf("%.2d", ili_idx$WEEK[i]),".htm")
    }
    ili_table <- try(readHTMLTable(theurl)[[1]], silent = TRUE)
    if(class(ili_table)=="try-error"){
      missing_data_id <- c(missing_data_id, j)
      next
    }
    for(colid in 2:ncol(ili_table)){
      ili_table[,colid] <- as.numeric(gsub("x", NA, ili_table[,colid]))
    }
    ili_today <- ili_table[,"Total ILI"]/ili_table[,"Total Patients"]*100

    ili_weektoday <-  as.numeric(as.character(ili_table$Week))

    if(ili_idx$DATE[i]>=as.Date("2014-10-04")){
      end_index <- tail(ili_weektoday,1)%%100-ili_idx$WEEK[i]
    }else{
      end_index <- tail(ili_weektoday,1)-ili_idx$WEEK[i]
    }
    #     if(!end_index%in%c(0,1))
    #       break
    if(ili.weighted)
      ili_today <- as.numeric(as.character(ili_table[,"% Weighted ILI"]))
    if(j-length(ili_today) + end_index<1)
      ili_today <- ili_today[-(1:(1-j+length(ili_today)-end_index))]
    ili_unrivised[(j-length(ili_today)):(j-1) + end_index,j] <- ili_today
  }

  missing_data_ili_idx_id <- sapply(missing_data_id, function(j) which(ili_idx$DATE == index(GFT_xts)[j-1]))
  ili_unrivised <- as.xts(ili_unrivised, order.by = index(GFT_xts))
  result <- list()
  result$ili_unrivised <- ili_unrivised
  result$missing_date <-  ili_idx[missing_data_ili_idx_id,]# ILI revision hisotry is available only during week 40 to week 20 of next year
  result
}

#' Parsing of Google Trends data
#' @param folder folder with weekly Google Trends file
#' @references
#' Yang, S., Santillana, M., & Kou, S. C. (2015). Accurate estimation of influenza epidemics using Google search data via ARGO. Proceedings of the National Academy of Sciences, \href{https://dx.doi.org/10.1073/pnas.1515373112}{doi: 10.1073/pnas.1515373112}.
#' @import xts
#' @import zoo
#'
#' @examples
#'
#' \donttest{
#' download.file("https://scholar.harvard.edu/files/syang/files/gt2016-10-24.zip",
#' file.path(tempdir(), "gt2016-10-24.zip"))
#' unzip(file.path(tempdir(), "gt2016-10-24.zip"), exdir = tempdir())
#' gt.folder <- file.path(tempdir(), "2016-10-19")
#' parsed_data <- parse_gt_weekly(gt.folder)
#' }
#'
#' @export
parse_gt_weekly <- function(folder){
  gtfiles <- list.files(folder)
  gtfiles <- gtfiles[which(sapply(strsplit(gtfiles, "[.]"), function(x) x[2])=="csv")]
  GTdata <- list()
  GTdata_week <- list()
  GTdata_month <- list()
  for(f in gtfiles){
    fscaned <- scan(file.path(folder, f),
                    what=character(), sep="\n", quiet=TRUE)
    if(length(fscaned) < 4)
      next

    endline <- which(!grepl(",", fscaned))[4]-1
    if(is.na(endline))
      endline <- length(fscaned)
    con <- textConnection(paste(fscaned[4:endline], collapse = '\n'))
    GTdata[[f]] <- read.csv(con)
    close(con)
    if(substr(fscaned[4],1,4)=="Week"){
      GTdata_week[[f]] <- xts(GTdata[[f]][,2,drop=FALSE], as.Date(GTdata[[f]]$Week)+6)
    }else if(substr(fscaned[4],1,4)=="Mont"){
      GTdata_month[[f]] <- xts(GTdata[[f]][,2,drop=FALSE], as.yearmon(GTdata[[f]]$Month))
    }
  }
  GTdata_week <- do.call(merge, GTdata_week)

  GTdata <- GTdata_week
  GTdata <- na.omit(GTdata)
  GTdata.print <- as.data.frame(GTdata)
  GTdata.print <- cbind(Week=rownames(GTdata.print), GTdata.print)
  list(GTdata=GTdata, GTdata.print=GTdata.print)
}

#' Parsing each Google Trends file downloaded from website
#' @param gt.folder folder that contains Google Trends file
#' @param f filename for Google Trends file
#' @import utils
#' @import zoo
gt.parser.pub.web <- function(gt.folder, f){
  fscaned <- scan(file.path(gt.folder, f), what=character(), sep="\n", quiet = TRUE)
  if(length(fscaned) < 4 || all(fscaned != "Interest over time"))
    return(list(state=NA, gtdata=NA, gtfreq=NA))
  state <- substr(fscaned[2], 1, gregexpr(" 2004", fscaned[2])[[1]]-1)
  endline <- which(!grepl(",", fscaned))[4]-1
  if(is.na(endline))
    endline <- length(fscaned)
  con <- textConnection(paste(fscaned[4:endline], collapse = '\n'))
  gtdata <- read.csv(con)
  if(substr(fscaned[4],1,4)=="Week"){
    gtdata <- xts(gtdata[,2,drop=FALSE], as.Date(gtdata$Week)+6)
    gtfreq <- "week"
  }else{
    gtdata <- xts(gtdata[,2,drop=FALSE], as.yearmon(gtdata$Month))
    gtfreq <- "month"
  }
  if(!grepl(colnames(gtdata), gsub(" ",".",f))){
    cat(f,colnames(gtdata),"\n")
  }
  close(con)
  return(list(state=state, gtdata=gtdata, gtfreq=gtfreq))
}

#' Parsing each Google Trends file downloaded from Google Trends API
#' @param gt.folder folder that contains Google Trends file
#' @param f filename for Google Trends file
#' @import utils
#' @import zoo
gt.parser.pub.api <- function(gt.folder, f){
  gtdata <- read.csv(file.path(gt.folder, f))
  f_info <- strsplit(f, "_")[[1]]
  state <- f_info[1]
  gtdata <- xts(gtdata[,ncol(gtdata),drop=FALSE], as.Date(gtdata$date))
  if(all(diff(index(gtdata))==7)){
    index(gtdata) <- index(gtdata) + 6
    gtfreq <- "week"
  }else{
    index(gtdata) <- as.yearmon(index(gtdata))
    gtfreq <- "month"
  }
  term <- colnames(gtdata)
  term <- gsub(".*\\.m\\.", "",term)
  if(!grepl(term, gsub(" ",".",f))){
    cat(f,colnames(gtdata),"\n")
  }
  return(list(state=state, gtdata=gtdata, gtfreq=gtfreq))
}


#' Parsing of raw data for regional ILI estimation
#'
#' @importFrom zoo index
#' @param gt.folder folder with all Google Trends data
#' @param ili.folder folder with all ILI data
#' @param population.file csv file path with state population data
#' @param gft.file csv file path for Google Flu Trends
#' @param gt.parser Google Trends data parser function, could be `gt.parser.pub.web` or `gt.parser.pub.api`
#'
#' @examples
#' \donttest{
#' download.file("https://scholar.harvard.edu/files/syang/files/gt2016-10-24.zip",
#' file.path(tempdir(), "gt2016-10-24.zip"))
#' unzip(file.path(tempdir(), "gt2016-10-24.zip"), exdir = tempdir())
#' gt.folder <- file.path(tempdir(), "2016-10-19")
#'
#' data_parsed <- load_reg_data(
#'   gt.folder=gt.folder,
#'   ili.folder=system.file("regiondata", "ili20161121", package = "argo"),
#'   population.file=system.file("regiondata", "Population.csv", package = "argo"),
#'   gft.file=system.file("regiondata", "GFT.txt", package = "argo")
#' )
#' }
#'
#'
#' @references
#' Shaoyang Ning, Shihao Yang, S. C. Kou. Accurate Regional Influenza Epidemics Tracking Using Internet Search Data. Scientific Reports
#'
#' @export
load_reg_data <- function(gt.folder, ili.folder, population.file, gft.file, gt.parser = gt.parser.pub.web) {
  states <- read.csv(population.file, stringsAsFactor=F)
  states$Population <- as.numeric(gsub(",","", states$Population))
  state.abb <- states$Abbre
  state.name <- states$State

  gtfiles <- list.files(gt.folder)
  gtfiles <- gtfiles[which(sapply(strsplit(gtfiles, "[.]"), function(x) x[2])=="csv")]

  GTdata_week <- list()
  GTdata_month <- list()
  state.info <- c()
  for(f in gtfiles){
    fparsed <- gt.parser(gt.folder, f)
    if(is.na(fparsed$state))
      next
    state.info[f] <- fparsed$state
    if(fparsed$gtfreq == "week"){
      GTdata_week[[f]] <- fparsed$gtdata
    }else{
      GTdata_month[[f]] <- fparsed$gtdata
    }
  }

  GTdata_week.state <- tapply(GTdata_week, state.info, function(gt.eachstate){
    tab <- do.call(merge, gt.eachstate)
    tab[,!grepl("1",colnames(tab))]
  })

  if(any(nchar(names(GTdata_week.state)) > 5)){
    names(GTdata_week.state) <- sapply(names(GTdata_week.state), function(x){
      if(x=="United States")
        return("US")
      paste0("US-",state.abb[gsub("\\s+$", "", strsplit(x, "\\(")[[1]][1]) == state.name])
    })
  }

  GTdata_week.state <- GTdata_week.state[order(names(GTdata_week.state))]

  #### read combined data file ####
  GT_national <- GTdata_week.state$US
  GT_national <- na.omit(GT_national)

  GT_data <- list()
  for(j in state.abb){
    empty_GT <- GT_national
    empty_GT[] <- 0
    GT_data[[j]] <- GTdata_week.state[[paste0("US-",j)]]
    empty_GT[,colnames(GT_data[[j]])] <- GT_data[[j]][index(GT_national)]
    GT_data[[j]] <- empty_GT
  }


  ili <- read.csv(file.path(ili.folder,"ILINet_regional.csv"), skip = 1)
  ili_idx <- subset(ili, ili$REGION=="Region 1")[,c("YEAR", "WEEK")]
  ili_idx$DATE <- as.Date("1997-10-04") + (1:nrow(ili_idx)-1)*7
  ili_idx$SEASON <- ifelse(ili_idx$WEEK < 40, paste0(ili_idx$YEAR-1,"-",ili_idx$YEAR),
                           paste0(ili_idx$YEAR,"-",ili_idx$YEAR+1))
  ili_regional <- list()
  for(j in 1:10){
    to_num <- as.character(ili$X..WEIGHTED.ILI)
    to_num[to_num=="X"] <- ""
    ili_regional[[paste0("Region.",j)]] <- as.numeric(to_num)[ili$REGION==paste("Region",j)]
  }
  ili_regional <- do.call(cbind, ili_regional)
  ili_regional <- xts(ili_regional, ili_idx$DATE)

  ili_national <- read.csv(file.path(ili.folder,"ILINet_nat.csv"), skip = 1)
  to_num <- as.character(ili_national$X..WEIGHTED.ILI)
  to_num[to_num=="X"] <- ""
  ili_national <- xts(as.numeric(to_num), as.Date("1997-10-04") + (1:nrow(ili_national)-1)*7)

  if(file.exists(file.path(ili.folder,"ILINet_state.csv"))){
    ili_state_raw <- read.csv(file.path(ili.folder,"ILINet_state.csv"), skip = 1)
    ili_state_list <- list()
    for(j in unique(ili_state_raw$REGION)){
      abbv <- state.abb[state.name == j]
      if(length(abbv) == 0){
        abbv <- j
      }
      ili_state_this <- ili_state_raw[ili_state_raw$REGION==j, c("REGION", "YEAR", "WEEK", "X.UNWEIGHTED.ILI")]
      ili_state_this$X.UNWEIGHTED.ILI <- as.numeric(as.character(ili_state_this$X.UNWEIGHTED.ILI))
      ili_state_this <- merge(ili_state_this, ili_idx, by=c("YEAR", "WEEK"))
      ili_state_xts <- xts(ili_state_this$X.UNWEIGHTED.ILI, ili_state_this$DATE)
      colnames(ili_state_xts) <- paste0("US-",abbv)

      ili_state_list[[paste0("US-",abbv)]] <- ili_state_xts
    }
    ili_state <- do.call(merge, ili_state_list)
    colnames(ili_state) <- names(ili_state_list)
    names(ili_state) <- gsub("-", ".", names(ili_state))
  }else{
    ili_state <- NULL
  }

  GT_regional <- lapply(1:10, function(region_number){
    states_id <- which(states$Region==region_number)
    GT_states_wgted <- lapply(states_id, function(j){
      GT_data[[states$Abbre[j]]]*states$Population[j]/sum(states$Population[states_id])
    })
    Reduce("+", GT_states_wgted)
  })
  names(GT_regional) <- paste0("Region.",1:10)

  gft.pred <- read.csv(gft.file, skip = 10)
  gft.pred <- xts(gft.pred[,-1],as.Date(gft.pred[,1])+6)
  gft.pred <- gft.pred[,grep("HHS.Region",colnames(gft.pred))]
  gft.pred <- gft.pred/1000
  names(gft.pred) <- sapply(strsplit(names(gft.pred), "[.][.]"), function(x) x[1])

  names(GTdata_week.state) <- gsub("-", ".", names(GTdata_week.state))
  names(GT_data) <- paste0("US.", names(GT_data))

  list(ili_national=ili_national, ili_regional=ili_regional,
       GT_national=GT_national, GT_regional=GT_regional, GFT=gft.pred,
       GT_state=GTdata_week.state, ili_state=ili_state, GT_state_filled=GT_data)
}
