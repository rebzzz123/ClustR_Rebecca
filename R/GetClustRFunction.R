#' @title Get ClustR
#' @description This function generates dynamically named CSVs with the IDs and disease statuses of subjects within disease clusters identified by clustr().
#' @param data dataframe: The name of the original dataset clustr() used to produce results.
#' @param id character: The name of the variable in data containing the unique identification number for subjects.
#' @param status character: The name of the variable in data containing the disease status of subjects, where 1 indicates cases and 0 indicates controls.
#' @param date character: The the name of the variable in data containing the date associated with disease status, such as diagnosis date or birth date. Date must be in either Y-m-d or Y/m/d format.
#' @param lat character: The name of the variable in data (as a character string) containing the latitude associated with each subject, such as diagnosis address or birth address.
#' @param long character: The name of the variable in data (as a character string) containing the longitude associated with each subject, such diagnosis address or birth address.
#' @param results dataframe: The name of the clustr() output.
#' @param type binary: Either 0 or 1 to indicate whether you want to map control clusters or case clusters, respectively. Defaults to 1.
#' @param adj logical: Indicates whether you want data from clusters based on chi-square values (unadjusted) or based on the Benjamini-Hochberg p-values (adjusted). Defaults to FALSE.
#' @param pctl.cutoff categorical (1,2,3): If you are not mapping clusters based on the adjusted p-value, you must specify what percentile cut-off to use. 1 indicates whether you want clusters with a chi-square value above the 100th percentile, 2 for above the 99th percentile, 3 for above the 97.5th percentile. Defaults to 1.
#' @param adj.cutoff numeric: If using post-adjustment clusters, indicate what the significance cut-off should be. Defaults to 0.05.
#' @param space_incs (Optional) numeric: A single number or vector of numbers representing what size (in decimal degrees) of clusters you want to pull subjects from.
#' @export
#' @importFrom stats chisq.test complete.cases ecdf p.adjust quantile setNames
#' @importFrom utils write.csv

get.clustr <- function(data=NULL, id=NULL, status=NULL, date=NULL, lat=NULL, long=NULL, results=NULL, type=1, adj=FALSE, pctl.cutoff=1, adj.cutoff=0.05, space_incs=NULL){
  if (is.null(data)){stop("You must specify the name of the dataset was used to generate results.")}
  if (is.null(id)){stop("You must specify the name of the id variable.")}
  if (is.null(status)){stop("You must specify the name of the disease status variable.")}
  if (is.null(lat)){stop("You must specify the name of the latitude variable.")}
  if (is.null(long)){stop("You must specify the name of the longitude variable.")}
  if (is.null(results)){stop("You must specify the name of the clustr() results.")}
  if (!(type %in% c(0,1))){stop("Type must be value 0 or 1.")}

  # Format data again, like we did for clustr(), particularly for the date variable
  overall_data <- data[,c(id,status,date,lat,long)]

  # Rename the columns for clearer use in later code.
  colnames(overall_data) <- c("id","status","date","lat","long")

  # Tell R that the birthdate values are dates.
  overall_data$D <- as.Date(overall_data$date)

  # Determine the minimum date in the dataset and create a new column repeating the minumum date at every row.
  min_date <- min(overall_data$D)
  overall_data$min_date <- min_date

  # Calculate the number of days that have passed for each person between the earliest date in the dataset and each person's personal birthdate and put this in a new column called *time*.
  overall_data$time <- as.numeric(difftime(overall_data$D,overall_data$min_date), units = "days")

  # Drop any person with missing data.
  overall_data <- overall_data[complete.cases(overall_data),]
  # Get rid of unnecessary variables.
  data_true <- overall_data[,c("id","status","time","date","lat","long")]

  #Generating dataframes
  if (adj==FALSE){
    if (pctl.cutoff==3){
      info <- subset(results, (excess==type & Percentile %in% c(">97.5th",">99th",">100th")))
      csvname_addendum <- " - Chi square greater than 97 5th percentile.csv"
    } else if (pctl.cutoff==2){
      info <- subset(results, (excess==type & Percentile %in% c(">99th", ">100th")))
      csvname_addendum <- " - Chi square greater than 99th percentile.csv"
    } else if (pctl.cutoff==1){
      info <- subset(results, (excess==type & Percentile==">100th"))
      csvname_addendum <- " - Chi square greater than 100th percentile.csv"
    }
  } else {
    info <- subset(results, (excess==type & FDR_p < adj.cutoff))
    text.adj.cutoff <- as.character(adj.cutoff)
    text.adj.cutoff <- substring(text.adj.cutoff, 3)
    csvname_addendum <- paste0(" - adjusted p less than point ",text.adj.cutoff,".csv")
  }
  if (!(is.null(space_incs))){
    info$Radius <- as.numeric(info$Radius)
    info <- subset(info, Radius %in% space_incs)
    csvname_middle <- " - radii of point "
    for (r in space_incs){
      text.r <- as.character(r)
      text.r <- substring(text.r, 3)
      csvname_middle <- paste0(csvname_middle, text.r, " ",sep="")
    }
    csvname_addendum <- paste0(csvname_middle, csvname_addendum)
  }
  if(nrow(info)==0){stop("There are no clusters with these specifications.")}
  for (i in 1:nrow(info)){
    lat <- info$lat[i]
    long <- info$long[i]
    time <- info$time[i]
    y <- info$chisq_max_time[i]
    if (type==1){
      csvname <- paste0("IDs for Case Cluster Around ID ",info$id[i], csvname_addendum)
    } else if (type==0){
      csvname <- paste0("IDs for Control Cluster Around ID ",info$id[i], csvname_addendum)
    }
    df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("id","status"))
    #cases
    caseids <- data_true$id[which((((data_true$lat-lat)^2 + (data_true$long-long)^2)^.5)<info$chisq_max_size[i] & (abs(data_true$time-time))<y & data_true$status==1)]
    for (c in caseids){
      index <- which(caseids==c)
      df[index,] <- c(c,1)
    }
    n <- nrow(df)
    #controls
    controlids <- data_true$id[which((((data_true$lat-lat)^2 + (data_true$long-long)^2)^.5)<info$chisq_max_size[i] & (abs(data_true$time-time))<y & data_true$status==0)]
    for (o in controlids){
      index <- n+which(controlids==o)
      df[index,] <- c(o,0)
    }
    df2 <- merge(df, data_true, by=c("id","status"))
    df2 <- subset(df2, select=-c(time))
    write.csv(df2, file=csvname, row.names=FALSE)
  }
}
