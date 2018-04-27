#'@title Map ClustR Results
#'@description This function maps the radii containing clusters identified by the clustr() function using ggmap.
#'@param results dataframe: The name of the clustr() output to map from.
#'@param type binary: Either 0 or 1 to indicate whether you want to map control clusters or case clusters, respectively. Defaults to 1.
#'@param adj logical: Indicates whether you want to plot clusters based on chi-square values (unadjusted) or based on the Benjamini-Hochberg p-values (adjusted). Defaults to FALSE.
#'@param adj.cutoff numeric: If mapping post-adjustment clusters, indicate what the significance cut-off should be. Defaults to 0.05.
#'@param zoom numeric: Tells ggmap how far to zoom into your region. Defaults to 6.
#'@param region character: Name of the region (e.g. state, city, country) where your clusters are.
#'@keywords clustr space time clustering ggmap ggplot2
#'@examples
#' \dontrun{
#'map.clustr(results=clustr_test5, type=1, adj=TRUE, adj.cutoff=.10, zoom=6, region='California')
#'}
#'@export
#'@import ggplot2 ggmap
#' @importFrom stats complete.cases setNames
#' @importFrom utils write.csv

map.clustr <- function(results=NULL, type=1, adj=FALSE, adj.cutoff=.05, zoom=6, region=NULL){
  #making sure specified vars
  if (is.null(results)){stop("You must specify what ClustR output to map.")}
  if (is.null(region)){stop("You must specify what region to map the clusters over.")}
  # Get basemap, zoom will need to be fiddled with by the user --> if google has trouble reading your query, will retry up to 1000 times until it gets the map and suppress unnecessary warnings each time
  testmap <- suppressWarnings(try(map <- get_map(location=region, zoom=zoom, maptype="roadmap", scale=2), silent=TRUE))
  suppressWarnings(
  try(
    for (i in 1:1000){
      if("try-error" %in% class(testmap) & i != 1000){
        testmap <- try(map <- get_map(location=region, zoom=zoom, maptype="roadmap", scale=2), silent=TRUE)
      } else {
        if (i==1000){
          stop("Reached query limit for gg map - try again later.")
        } else {
          break
        }
      }
    }
  , silent=TRUE))

  # Code from StackFlow, modified, for generating circles for ggmap
  circularise <- function(d, n=360){
    angle <- seq(-pi, pi, length = n)
    make_circle <- function(x,y,r,id){
      data.frame(x=x+r*cos(angle), y=y+r*sin(angle), id)
    }
    lmat <- mapply(make_circle, id = seq_len(nrow(d)),
                   x = d[,1], y=d[,2], r=d[,3], SIMPLIFY = FALSE)
    do.call(rbind, lmat)
  }
  #--------------------------------------------Graph the radii with significantly high proportion of controls--------------------------------------------#
  if (type==1 & adj==FALSE){
  # Redid all of the plotting to more accurately reflect where ClustR is looking
  casemaptitle <- paste0("Clusters of Cases, ",region)

  # Prepare data
  subdata <- subset(results, Percentile!="none" & excess==1)
  subdata$Radius <- as.numeric(subdata$Radius)
  subdata$RadiusFact <- as.factor(subdata$Radius)

  # Create radius data
  circles <- setNames(data.frame(matrix(ncol=4, nrow=0)),c("x","y","id","Percentile")) #build empty dataframe to add data to
  for (x in unique(subdata$Percentile)){
    subsubdata <- subset(subdata, Percentile==x, select=c(long, lat, Radius))
    colnames(subsubdata) <- c("x", "y", "r")
    groupcircles <- circularise(subsubdata) #creating radius data for each percentile, so we can keep track of what data came from what group for later
    groupcircles$Percentile <- x
    groupcircles$id <- groupcircles$id + nrow(circles) #make sure we don't repeat id numbers
    circles <- rbind(circles, groupcircles)
  }

  # Make sure levels are in the correct order
  circles$Percentile <- factor(circles$Percentile, levels=c(">100th",">99th",">97.5th"), labels=c(">100th",">99th",">97.5th"))

  # Making your map
  plotmap <- ggmap(map)+ geom_polygon(aes(x,y,group=id,colour=Percentile), fill=NA, data=circles, alpha=.5) + coord_fixed() + scale_colour_manual(name="Percentile", values=c(">97.5th"="darkorange", ">99th"="red", ">100th"="darkmagenta")) + ggtitle(casemaptitle) + xlab("Longitude") + ylab("Latitude")
  return(plotmap)
  }

  #--------------------------------------------Graph the radii with significantly high proportion of controls--------------------------------------------#
  if (type==0 & adj==FALSE){
  controlmaptitle <- paste0("Clusters of Controls, ",region)

  # Prepare data
  subdata <- subset(results, Percentile!="none" & excess==0)
  subdata$Radius <- as.numeric(subdata$Radius)
  subdata$RadiusFact <- as.factor(subdata$Radius)

  # Create radius data
  circles <- setNames(data.frame(matrix(ncol=4, nrow=0)),c("x","y","id","Percentile")) #build empty dataframe to add data to
  for (x in unique(subdata$Percentile)){
    subsubdata <- subset(subdata, Percentile==x, select=c(long, lat, Radius))
    colnames(subsubdata) <- c("x", "y", "r")
    groupcircles <- circularise(subsubdata) #creating radius data for each percentile, so we can keep track of what data came from what group for later
    groupcircles$Percentile <- x
    groupcircles$id <- groupcircles$id + nrow(circles) #make sure we don't repeat id numbers
    circles <- rbind(circles, groupcircles)
  }

  # Make sure levels are in the correct order
  circles$Percentile <- factor(circles$Percentile, levels=c(">100th",">99th",">97.5th"), labels=c(">100th", ">99th", ">97.5th"))

  # Making your map
  plotmap <- ggmap(map)+ geom_polygon(aes(x,y,group=id,colour=Percentile), fill=NA, data=circles, alpha=.5) + coord_fixed() + scale_colour_manual(name="Percentile", values=c(">100th"="darkmagenta", ">99th"="red", ">97.5th"="darkorange")) + ggtitle(controlmaptitle) + xlab("Longitude") + ylab("Latitude")
  return(plotmap)
  }

  #--------------------------------------------Case Mapping after B-H Adjustment--------------------------------------------#
  if (type==1 & adj==TRUE){
  # Prepare data
  subdata <- subset(results, FDR_p < adj.cutoff & excess==1)

    if (nrow(subdata)==0){
      stop (paste0("There are no case clusters with an FDR adjusted p-value < ",adj.cutoff," to graph."))
    } else {
    casemaptitle <- paste0("Clusters of Cases after Benjamini-Hochberg Adjustment, ",region,", p<",adj.cutoff)

    subdata$Radius <- as.numeric(subdata$Radius)
    subdata$RadiusFact <- as.factor(subdata$Radius)

    # Creating dataset for this function
    subsubdata <- subset(subdata, select=c(long, lat, Radius))
    colnames(subsubdata) <- c("x", "y", "r")

    # Create radius data
    circles <- circularise(subsubdata)

    # Making your map
    plotmap <- ggmap(map)+ geom_polygon(aes(x,y,group=id), colour="red", fill=NA, data=circles, alpha=.5) + coord_fixed() + ggtitle(casemaptitle) + xlab("Longitude") + ylab("Latitude")
    return(plotmap)
    }
  }
  #--------------------------------------------Control Mapping after B-H Adjustment--------------------------------------------#
  if (type==0 & adj==TRUE){
  #prepare data
  subdata <- subset(results, FDR_p < adj.cutoff & excess==0)

  if (nrow(subdata)==0){
    stop (paste0("There are no control clusters with an FDR adjusted p-value < ",adj.cutoff," to graph."))
  } else {
    controlmaptitle <- paste0("Clusters of Controls after Benjamini-Hochberg Adjustment, ",region,", p<",adj.cutoff)

    subdata$Radius <- as.numeric(subdata$Radius)
    subdata$RadiusFact <- as.factor(subdata$Radius)

    #creating dataset for this function
    subsubdata <- subset(subdata, select=c(long, lat, Radius))
    colnames(subsubdata) <- c("x", "y", "r")

    #create radius data
    circles <- circularise(subsubdata)

    #making your map
    plotmap <- ggmap(map)+ geom_polygon(aes(x,y,group=id), colour="red", fill=NA, data=circles, alpha=.5) + coord_fixed() + ggtitle(controlmaptitle) + xlab("Longitude") + ylab("Latitude")
    return(plotmap)
    }
  }
}
