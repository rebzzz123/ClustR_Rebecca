#' @title ClustR
#' @description This function evaluates space-time clustering of case-control or cohort data. ClustR simulates a random distribution of disease for the input data and then compares the true data against it by bootstrapping samples and comparing chi-square distributions. This function returns a dataframe summarizing which areas and timeframes of the true data were sampled and which yielded unusually high proportions of cases or controls compared to the simulation. It also returns p-values for each sample (generated_p) and then adjusts for false discovery rate via the Benjamini-Hochberg approach (FDR_p).
#' @param id character: The name of the variable containing the unique identification number for subjects.
#' @param status character: The name of the variable containing the disease status of subjects, where 1 indicates cases and 0 indicates controls.
#' @param date character: The the name of the variable containing the date associated with disease status, such as diagnosis date or birth date. Date must be in either Y-m-d or Y/m/d format.
#' @param lat character: The name of the variable (as a character string) containing the latitude associated with each subject, such as diagnosis address or birth address.
#' @param long character: The name of the variable (as a character string) containing the longitude associated with each subject, such diagnosis address or birth address.
#' @param time_incs (Optional) numeric: A single number or vector of numbers representing the maximum difference in time (in days) between subjects in clusters. Defaults to c(90,180,365). If a vector of numbers is given, ClustR will try each constraint and return the most significant cluster among them.
#' @param space_incs (Optional) numeric: A single number or vector of numbers representing the maximum difference in space (in decimal degrees) between subjects in clusters. Defaults to c(.025, .05, .1, .5). If a vector of numbers is given, ClustR will try each constraint and return the most significant cluster among them.
#' @param n (Optional) numeric: A single number representing the size of samples you want ClustR to bootstrap. Defaults to nrow(data) if nrow(data) < 1001 or else defaults to 1000.
#' @param bootstraps (Optional) numeric: A single number representing how many samples ClustR should bootstrap. Defaults to 5.
#' @param data dataframe: The name of the dataset ClustR should evaluate.
#' @keywords space time clustering
#' @examples
#' \dontrun{
#' clustr_test1 <- clustr(id="blindid",status="case_status",date="birthdate",lat="mom_lat",
#'        long="mom_long",time_incs=c(30,90),space_incs=c(.5,1),n=500,bootstraps=20,data=testdata)
#' }
#' @export
#' @importFrom stats chisq.test complete.cases ecdf p.adjust quantile setNames
#' @importFrom utils write.csv

clustr <- function(id=NULL, status=NULL, date=NULL, lat=NULL, long=NULL, time_incs=c(90,180,365), space_incs=c(.025, .05, .1, .5), n=NULL, bootstraps=5, data){

  # For greater precision, upped digits to 10.
  options(digits=10)
  combined_sample_result <- as.data.frame(matrix())
  combined_total <- matrix(,1,1)

  # Format dataset.

  # Limit the dataset to five variables in the following order: ID, case status, date (in YYYY-MM-DD format), latitude and longitude.
  # Based on what the user told R, ClustR extracts these variables and puts them in the order it wants.
  overall_data <- data[,c(id,status,date,lat,long)]

  # Rename the columns for clearer use in later code.
  colnames(overall_data) <- c("id","status","date","lat","long")

  # Checking if options are correctly specified
  if (is.null(id)) stop("You must specify an identification variable.")
  if (is.null(status)) stop("You must specify an disease status variable.")
  if (is.null(date)) stop("You must specify a date variable.")
  if (is.null(lat)) stop("You must specify a latitude variable.")
  if (is.null(long)) stop("You must specify an longitude variable.")
  if (TRUE %in% (sort(unique(overall_data$status))!=c(0,1))) stop("Status must be coded as 0 and 1 for controls and cases, respectively.")
  if (is.null(data)) stop("You must specify what dataframe to use.")


  # Tell R that the birthdate values are dates.
  overall_data$D <- as.Date(overall_data$date)

  # Determine the minimum date in the dataset and create a new column repeating the minumum date at every row.
  min_date <- min(overall_data$D)
  overall_data$min_date <- min_date

  # Calculate the number of days that have passed for each person between the earliest date in the dataset and each person's personal birthdate and put this in a new column called *time*.
  overall_data$time <- as.numeric(difftime(overall_data$D,overall_data$min_date), units = "days")

  # Drop any person with missing data.
  overall_data <- overall_data[complete.cases(overall_data),]

  # Evaluate User Options

  # Limit how many time ranges the user can look at.
  try(if(length(time_incs)>5) stop("You cannot specify more than five time ranges to look at."))

  # Get information about what the starting and ending dates for the data are.
    date_range = c(as.character(overall_data$date[which.min(overall_data$time)]), as.character(overall_data$date[which.max(overall_data$time)]))
    mintime <- min(overall_data$time)
    maxtime <- max(overall_data$time)
    time_range_vector <- mintime:maxtime

  # Drop the variables for birthdate and minimum date. (Now that we're done with them)
  data_true <- overall_data[,c("id","status","time","lat","long")]

  # If the radii (spatial increments) to be used are inappropriate for the data, the user will get an error.
    radius_incs <- space_incs
    groups <- length(radius_incs)
    try(if(groups>5) stop("You cannot specify more than five radii."))
    try(if((range(data_true$lat)[2]-range(data_true$lat)[1])<(max(radius_incs)*2)) stop("At least one of the radii specified by the user or the default is too large for your data."))
    try(if((range(data_true$long)[2]-range(data_true$long)[1])<(max(radius_incs)*2)) stop("At least one of the radii specified by the user or the default is too large for your data."))


  #The size of samples to bootstrap.

  # If the user did not specify the size of the samples to take: if sample size is small, boostrap as large of samples as possible; if sample size is large, boostrap samples of 1000.
  # If the user did specify the size of the samples to take: use that number.
  if (is.null(n)){
    if(length(data_true$id) < 1001) {
      samples <- nrow(data)
    } else {
      samples <- 1000
    }
  } else {
    samples <- n
  }

  # Permutated Data
  # Create a 1x1 matrix called *total* to hold the null distribution.

  total <- matrix(,1,1)

  ## Repeat the following process [bootstraps] times in order to get a stable distribution:

  for(c in 1:bootstraps){
    # Create a new dataset called *data_perm* with case/control status scrambled randomly across participants (using the same proportion as found in the original dataset). This is your permuted dataset.
    data_perm <- data_true
    cases_perm <- sum(data_true$status)  # generate number of cases (ie 20% of # of controls)
    controls_perm <- nrow(data_true) - cases_perm # generate number of controls
    assignment_perm <- sample(c(rep(1,cases_perm),rep(0,controls_perm)), replace = FALSE)  # scramble case and control ordering
    data_perm$status <- assignment_perm  # set permuted status
    # Rename *data_perm* as *data*.
    data <- data_perm
    # Number each row of the permuted dataset starting at 1 and ending at [number of rows in the dataset] in a new column called *num_id*.
    r <- nrow(data)
    data$num_id <- seq(1,r)
    # Create a new dataframe called *cluster* that includes the permuted dataset plus 5*[number of radii specified] empty columns.
    cluster <- cbind(data,matrix(,nrow = r, ncol = 5*groups*length(time_incs)))
    # Name the empty columns near_cases, near_controls, total, distance, time repeating per number of radii specified and time intervals specified.
    # Created dynamic column names that will reflect what radii and time increments are being used.
    cluster_col_names <- c(colnames(data))
    for (g in radius_incs){
      for (t in time_incs){
        name1 <- paste("near_cases", g, t, sep="_")
        name2 <- paste("near_controls", g, t, sep="_")
        name3 <- paste("total", g, t, sep="_")
        name4 <- paste("distance", g, t, sep="_")
        name5 <- paste("time", g, t, sep="_")
        radius_time_cols <- c(name1, name2, name3, name4, name5)
        cluster_col_names <- append(cluster_col_names, radius_time_cols)
      }
    }
    colnames(cluster) <- cluster_col_names

    # Create an empty [number of samples] x 6+5*[number of radii]*[number of time increments] matrix *new_tab* with the same column names as the new dataframe *cluster* for storing soon-to-be generated spatial and temporal data.
    new_tab <- matrix(,nrow = samples, ncol = (5*groups*length(time_incs)+6)) #CE: made the ncol of this matrix flexible for changing amount of radii
    colnames(new_tab) <- colnames(cluster)
    ### For [number of samples] times (main loops for permutation, uses replacement):
    for(a in 1:samples){
      # Pick a random row of the permuted data and take their latitude, longitude and time.
      i <- sample(r, 1)
      lat <- (data[i,"lat"])
      long <- (data[i,"long"])
      time <- (data[i,"time"])
      # **Then for [your number of radii]:**
      loopnumber <- 0
      for(p in 1:groups){
        # Then for [your number of time increments]:
        for (u in 1:length(time_incs)){
          # Replaced for loops here and replaced with WHICH statements: improved speed of this chunk by 93%
          # ClustR notes how many cases and controls are nearby this point in this radius in this time range.
          near_cases <- length(data$status[which((((data$lat-lat)^2 + (data$long-long)^2)^.5)<radius_incs[p] & (abs(data$time-time))<time_incs[u] & data$status==1)])
          near_controls <- length(data$status[which((((data$lat-lat)^2 + (data$long-long)^2)^.5)<radius_incs[p] & (abs(data$time-time))<time_incs[u] & data$status==0)])
          # For each set of five columns for each radius and time, for each randomly chosen observation, in the following order: input the number of nearby cases, the number of nearby controls, the total number of nearby observations of either status, what radius was being used and what the randomly chosen time interval was.
          loopnumber <- loopnumber + 1
          n = (loopnumber-1)*5
          cluster[i,7+n] <- near_cases
          cluster[i,8+n] <- near_controls
          cluster[i,9+n] <- near_cases + near_controls
          cluster[i,10+n] <- radius_incs[p]
          cluster[i,11+n] <- time_incs[u]
        } #(End of For Loop for [number of time increments])
      }
      # (End of For Loop for [your number of radii].)
      # For each randomly chosen row, take the complete row you just generated above and place that into the matrix called *new_tab*.
      entry <- as.matrix(cluster[i,])
      new_tab[a,] <- entry
    }
    # (End of For Loop for [number of samples] times.)
    # Convert *new_tab* to a dataframe and rename it as *sample_result*.
    sample_result <- as.data.frame(new_tab)
    # Create an empty matrix *chisq* which has the number of rows as *sample_result* but with ten extra columns at the end for holding chi-square information.
    chisq <- matrix(,nrow = nrow(sample_result), ncol = 2*groups*length(time_incs)) # Number of columns is flexible with changing # of radii and time incs
    # Create a vector called *baseline* that contains the total number of cases in the permuted dataset and the total number of rows in the dataset.
    cases <- sum(data$status)
    baseline <- c(cases, r)
    ### For each row in your *sample_result* matrix:
    for(k in 1:nrow(sample_result)){
      # **For each radii and time interval**
      loopnumber <- 0
      for(l in 1:groups){
        for (u in 1:length(time_incs)){
          # Compare the number of cases and total number of observations nearby each row (*clust*) to the number of cases and total number of observations in the whole dataset (*baseline*).
          loopnumber <- loopnumber + 1
          tcol <- (loopnumber*5)+4
          ccol <- (loopnumber*5)+2
          clust <- c(sample_result[k,ccol], (sample_result[k,tcol]))
          compare <- baseline - clust
          test_table <- cbind(compare,clust)
          # Run a chi-square test comparing proportion of cases in sampled result to proportion of cases in permuted dataset.
          suppressWarnings(cs <- chisq.test(test_table)) #CE: suppresses the thousands of warnings that chi-square values may be too small to give reliable p-values
          # Input the chi-square statistic and p-value for each test/row for each radii and for each time into the *chisq* matrix.
          chisq[k,loopnumber] <- cs$statistic
          chisq[k,loopnumber+(groups*length(time_incs))] <- cs$p.value #CE: changed 5 to groups here to make code more flexible for varying # of radii
        }}}
    # (End of For Loop for each radii and each time interval.)
    # (End of For Loop for each row in the *sample_result* matrix.)
    # Put these chi-square statistics in your main *sample_result* matrix for your 1000 sampled points.
    loopnumber = 0
    for (g in 1:groups){
      for (u in 1:length(time_incs)){
        loopnumber = loopnumber + 1
        sample_result$new <- chisq[,loopnumber]
        names(sample_result)[names(sample_result)=='new'] <- paste("chisq", radius_incs[g], time_incs[u], sep="_")
      }
    }
    # For each row for each sampled person, take the maximum chi-square value across the radii and put that into a new column in your *sample_result* matrix.
    sample_result$chisq_max <- apply(sample_result[,grep("chisq_", names(sample_result), value=TRUE)], 1, max)
    # Append the column of the highest chi-square values to your empty 1x1 matrix called *total*.
    total <- sample_result$chisq_max

    # Add the results of this single bootstrap to dataframes for holding info from all bootstraps together.
    # If this is the first loop in the for loop, need to adjust the combined_sample_result matrix to fit the size of the sample_result matrix and to have its column names:
    if (c==1){
      combined_sample_result <- cbind(combined_sample_result,matrix(,nrow = 1, ncol = ncol(sample_result)))
      combined_sample_result <- combined_sample_result[,-1]
      colnames(combined_sample_result) <- colnames(sample_result)
      combined_sample_result <- rbind(combined_sample_result, sample_result)
      combined_sample_result <- combined_sample_result[-1,]
    } else {
      combined_sample_result <- rbind(combined_sample_result, sample_result)
    }
    combined_total <- append(combined_total, total)
  }
  combined_total <- combined_total[-1]

  rm(data_perm, cases_perm, controls_perm, assignment_perm)

  # (End of For Loop for repeating this process [bootstaps] times.)
  # Calculate and save the cut offs for the 97.5th, 99th, and 100th percentile of your chi-square values based on your randomly sampled points from your null distribution.

  # Quantiles are based on the distribution of chi-squares over every bootstrap performed.
  quantile(combined_total, c(.975,.99, 1), na.rm = TRUE)
  nsfpct <- quantile(combined_total, .975, na.rm = TRUE)
  nnpct <- quantile(combined_total, .99, na.rm  = TRUE)
  max <- quantile(combined_total, 1, na.rm = TRUE)


  # Real Data
  # Now repeat all of the above steps for the real data - the dataset where the case status is not scrambled. Differences noted below:

  data <- data_true

  # Reset this object for the Real Data.
  combined_sample_result <- as.data.frame(matrix())

  r <- nrow(data)
  data$num_id <- seq(1,r)

  groups = length(radius_incs)

  for (c in 1:bootstraps){
    cluster <- cbind(data,matrix(,nrow = r, ncol = 5*groups*length(time_incs)))
    cluster_col_names <- c(colnames(data))
    for (g in radius_incs){
      for (t in time_incs){
        name1 <- paste("near_cases", g, t, sep="_")
        name2 <- paste("near_controls", g, t, sep="_")
        name3 <- paste("total", g, t, sep="_")
        name4 <- paste("distance", g, t, sep="_")
        name5 <- paste("time", g, t, sep="_")
        radius_time_cols <- c(name1, name2, name3, name4, name5)
        cluster_col_names <- append(cluster_col_names, radius_time_cols)
      }
    }
    colnames(cluster) <- cluster_col_names

    new_tab <- matrix(,nrow = samples, ncol = (5*groups*length(time_incs)+6))
    colnames(new_tab) <- colnames(cluster)


    # MAIN LOOP

    for(a in 1:samples){
      i <- sample(r, 1)
      lat <- (data[i,"lat"])
      long <- (data[i,"long"])
      time <- (data[i,"time"])

      loopnumber = 0
      for(p in 1:groups){
        for (u in 1:length(time_incs)){
          near_cases <- length(data$status[which((((data$lat-lat)^2 + (data$long-long)^2)^.5)<radius_incs[p] & (abs(data$time-time))<time_incs[u] & data$status==1)])
          near_controls <- length(data$status[which((((data$lat-lat)^2 + (data$long-long)^2)^.5)<radius_incs[p] & (abs(data$time-time))<time_incs[u] & data$status==0)])

          loopnumber = loopnumber + 1
          n = (loopnumber - 1)*5
          cluster[i,7+n] <- near_cases
          cluster[i,8+n] <- near_controls
          cluster[i,9+n] <- near_cases + near_controls
          cluster[i,10+n] <- radius_incs[p]
          cluster[i,11+n] <- time_incs[u]
        }
      }
      entry <- as.matrix(cluster[i,])
      new_tab[a,] <- entry
    }

    sample_result <- as.data.frame(new_tab)
    loopnumber <- 0
    for (g in 1:groups){
      for (u in 1:length(time_incs)){
        loopnumber <- loopnumber + 1
        tcol <- (loopnumber*5)+4
        ccol <- (loopnumber*5)+2
        sample_result$new <- sample_result[,ccol]/sample_result[,tcol]
        name <- paste("case_rate", radius_incs[g], time_incs[u], sep="_")
        names(sample_result)[names(sample_result) == 'new'] <- name
      }
    }
    chisq <- matrix(,nrow = nrow(sample_result), ncol = 2*groups*length(time_incs))
    cases <- sum(data$status)
    baseline <- c(cases, r)

    for(k in 1:nrow(sample_result)){
      loopnumber <- 0
      for(l in 1:groups){
        for (u in 1:length(time_incs)){
          loopnumber <- loopnumber + 1
          tcol <- (loopnumber*5)+4
          ccol <- (loopnumber*5)+2
          clust <- c(sample_result[k,ccol], (sample_result[k,tcol]))
          compare <- baseline - clust
          test_table <- cbind(compare,clust)
          suppressWarnings(cs <- chisq.test(test_table)) # Suppresses the thousands of warnings that p-values may be unreliable due to a small amount of subjects in one or more boxes of your 2x2 table
          chisq[k,loopnumber] <- cs$statistic
          chisq[k,loopnumber+(groups*length(time_incs))] <- cs$p.value
        }}}

    loopnumber <- 0
    for (g in 1:groups){
      for (u in 1:length(time_incs)){
        loopnumber <- loopnumber + 1
        sample_result$new <- chisq[,loopnumber]
        names(sample_result)[names(sample_result)=='new'] <- paste("chisq", radius_incs[g], time_incs[u], sep="_")
      }
    }

    # Grab all chi square values across all the different radii and times.
    rel <- sample_result[,grep("chisq_", names(sample_result), value=TRUE)]

    # Also create matrices of the other characteristics of interest across all the different radii and times.
    rel_size <- sample_result[,grep("distance_", names(sample_result), value=TRUE)]
    rel_time <- sample_result[,grep("time_", names(sample_result), value=TRUE)]
    rel_caserate <- sample_result[,grep("case_rate_", names(sample_result), value=TRUE)]
    rel_casecount <- sample_result[,grep("near_cases_", names(sample_result), value=TRUE)]
    rel_controlcount <- sample_result[,grep("near_controls_", names(sample_result), value=TRUE)]
    sample_result$chisq_max <- apply(sample_result[,grep("chisq_", names(sample_result), value=TRUE)], 1, max)

    # Create a vector that indicates, for each row of sample result, which chi-square column has the largest chi-square value.
    rel_index <- apply(rel,1, which.max)

    # Obtain the characteristics of the clusters with the largest chi-square values from the real data. Add that information to the *sample_result* matrix.
    rez <- matrix(,ncol = 1, nrow = nrow(rel_size))
    rez_time <- matrix(,ncol = 1, nrow = nrow(rel_size))
    rez_rate <- matrix(,ncol = 1, nrow = nrow(rel_size))
    rez_count <- matrix(,ncol = 1, nrow = nrow(rel_size))
    rez_count_cont <- matrix(,ncol = 1, nrow = nrow(rel_size))
    for(i in 1:nrow(rel_size)){
      rez[i] <- rel_size[i,rel_index[i]]
      rez_time[i] <- rel_time[i, rel_index[i]]
      rez_rate[i] <- rel_caserate[i,rel_index[i]]
      rez_count[i] <- rel_casecount[i,rel_index[i]]
      rez_count_cont[i] <- rel_controlcount[i,rel_index[i]]
    }

    sample_result$chisq_max_size <- rez
    sample_result$chisq_max_time <- rez_time
    sample_result$chisq_max_caserate <- rez_rate
    sample_result$chisq_max_casecount <- rez_count
    sample_result$chisq_max_controlcount <- rez_count_cont

    # Calculate the bounds of the time interval covered by each cluster by subtracting and adding the time increment from the randomly chosen point's date.
    sample_result$clust_start_date <- sample_result$time-sample_result$chisq_max_time #CE: Adjusted due to new time sampling
    sample_result$clust_end_date <- sample_result$time+sample_result$chisq_max_time

    # Create a new variable in sample_results that indicates 0 is case rate is less than the overall caserate and 1 if case rate is >= the overall caserate.
    overallcaserate <- sum(overall_data$status)/length(overall_data$id)
    sample_result$excess <- ifelse(sample_result$chisq_max_caserate < overallcaserate, 0, 1)

    # If the maximum chi-square value in each row meets any one of the percentile cut-offs, mark it as so in a new column. If not, mark as "none." Make this a factor variable. (These will be used for plotting purposes.)
    sample_result$category <- ifelse(sample_result$chisq_max > max, ">100th", ifelse(sample_result$chisq_max > nnpct, ">99th", ifelse(sample_result$chisq_max > nsfpct, ">97.5th", "none")))
    sample_result$Percentile <- factor(sample_result$category, levels =  c(">100th", ">99th", ">97.5th", "none"), ordered = TRUE)

    # Create a new factor variable with what radius had the largest chi-square value.
    sample_result$Radius <- as.factor(sample_result$chisq_max_size)

    #Getting rid of unimportant variables
    dropping <- c("category","num_id")
    sample_result <- sample_result[,!(names(sample_result) %in% dropping)]

    # Adding the results of this single bootstrap to dataframes for holding info from all bootstraps together
    # If this is the first loop in the for loop, need to adjust the combined_sample_result matrix to fit the size of the sample_result matrix and to have its column names:
    if (c==1){
      combined_sample_result <- cbind(combined_sample_result,matrix(,nrow = 1, ncol = ncol(sample_result)))
      combined_sample_result <- combined_sample_result[,-1]
      colnames(combined_sample_result) <- colnames(sample_result)
      combined_sample_result <- rbind(combined_sample_result, sample_result)
      combined_sample_result <- combined_sample_result[-1,]
    } else {
      combined_sample_result <- rbind(combined_sample_result, sample_result)
    }

  }
  combined_sample_result$clust_start_date <- as.Date(combined_sample_result$clust_start_date, origin=min_date)
  combined_sample_result$clust_end_date <- as.Date(combined_sample_result$clust_end_date, origin=min_date)

  #(End of bootstrapping real data for loop.)

  #Generate p-values and adjusted-pvalues.
  percentile <- ecdf(combined_total) # Creates cumulative distribution function
  # Creates a new variable with p-values for each chi-square value compared to the null distribution.
  combined_sample_result$generated_p <- sprintf("%.50f",1-percentile(combined_sample_result$chisq_max)) # Uses 50 decimal places for this and below
  combined_sample_result$number_p <- as.numeric(combined_sample_result$generated_p)
  x <- 1-percentile(quantile(combined_total, .999, na.rm  = TRUE)) # Gets p-value for 99.9th percentile
  combined_sample_result$generated_p[which(combined_sample_result$number_p==0)] <- sprintf("%.50f",x/combined_sample_result$chisq_max[which(combined_sample_result$number_p==0)])  # Manually sets p-value for "off the charts" (non-existant in the simulated data) chi-square values
  # Creates a new variable with p-values adjusted for multiple testing.
  combined_sample_result$FDR_p <- sprintf("%.50f",p.adjust(combined_sample_result$generated_p, method = "fdr", n = nrow(combined_sample_result)))
  combined_sample_result <- combined_sample_result[,-which(names(combined_sample_result) %in% c("number_p"))]

  # Returns our dataframe with all results
  return(combined_sample_result)
}
