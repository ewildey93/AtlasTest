# based on A guide to pre-processing high-throughput animal tracking data, Gupte et al. 2022
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13610
library(atlastools)
library(sf)
library(amt)
setwd("C:/Users/eliwi/OneDrive/Documents/R/AtlasTest")

#prepare data for functions
AtlTest <- read.csv("./Est.Locations.Filter.Distance.250_Results.csv")
AtlTest1 <- read.csv("./1minEst.Locations.Filter.Distance.250_Results.csv")
str(AtlTest)
AtlTest$DateTime <- as.POSIXct(paste(AtlTest$Date, AtlTest$Time),
                               format="%Y-%m-%d %H:%M:%S", tz="GMT")
AtlTest1$DateTime <- as.POSIXct(paste(AtlTest1$Date, AtlTest1$Time),
                               format="%Y-%m-%d %H:%M:%S", tz="GMT")
attributes(AtlTest$DateTime)$tzone <- "America/Denver"
attributes(AtlTest1$DateTime)$tzone <- "America/Denver"
AtlTest$UnixTime <- as.numeric(AtlTest$DateTime)
AtlTest1$UnixTime <- as.numeric(AtlTest1$DateTime)
#colnames(AtlTest)[c(2,5,6,7,13)] <- c("TAG", "NBS", "X", "Y", "TIME")
write.csv(AtlTest, "./AtlTest.csv")

#get SD from confidence interval
AtlTest$VarX <- (((AtlTest$x.ci.upper-AtlTest$x.ci.lower)/3.92)*sqrt(AtlTest$No.Nodes))^2
AtlTest$VarY <- (((AtlTest$y.ci.upper-AtlTest$y.ci.lower)/3.92)*sqrt(AtlTest$No.Nodes))^2
AtlTest$CoV <- cov(AtlTest$x.est, AtlTest$y.est)

###############filter for time period of track####################
range(AtlTest$DateTime)
range(AtlTest1$DateTime)
a <- AtlTest%>%filter(DateTime > "2021-08-04 08:30:00" & DateTime < "2021-08-04 13:15:00")
range(a$DateTime)
write.csv(a, "C:/Users/eliwi/OneDrive/Documents/R/AtlasTest/FilterBeepTrackPts.csv")
aList <- split(a, f=a$TagId)
aLines <- lapply(aList, function (x) track(x, x = x$x.est, y = x$y.est, .t=x$DateTime,  
                                               crs=CRS("+init=EPSG:32613"), all_cols=TRUE))
aLines <- lapply(aLines, function (x) as_sf_lines(x))
for (i in 1:length(aLines)){st_write(aLines[[i]], dsn=paste0(names(aLines)[i], "filtertimetrack.shp"))}

#1min data
b <- AtlTest1%>%filter(DateTime > "2021-08-04 08:30:00" & DateTime < "2021-08-04 13:15:00")
write.csv(b, "C:/Users/eliwi/OneDrive/Documents/R/AtlasTest/FilterBeepTrackPts1min.csv")
bList <- split(b, f=b$TagId)
bLines <- lapply(bList, function (x) track(x, x = x$x.est, y = x$y.est, .t=x$DateTime,  
                                           crs=CRS("+init=EPSG:32613"), all_cols=TRUE))
bLines <- lapply(bLines, function (x) as_sf_lines(x))
for (i in 1:length(bLines)){st_write(bLines[[i]], dsn=paste0(names(bLines)[i], "filtertimetrack1min.shp"))}

###################################################################
#make lines for points
str(AtlTest)
AtlList <- split(AtlTest, f=AtlTest$TagId)
AtlList1 <- split(AtlTest1, f=AtlTest1$TagId)

AtlLines <- lapply(AtlList, function (x) track(x, x = x$x.est, y = x$y.est, .t=x$DateTime,  
                                               crs=CRS("+init=EPSG:32613"), all_cols=TRUE))
AtlLines1 <- lapply(AtlList1, function (x) track(x, x = x$x.est, y = x$y.est, .t=x$DateTime,  
                                               crs=CRS("+init=EPSG:32613"), all_cols=TRUE))
AtlLines <- lapply(AtlLines, function (x) as_sf_lines(x))
AtlLines1 <- lapply(AtlLines1, function (x) as_sf_lines(x))

names <- names(AtlLines1)
for (i in 1:length(AtlLines1)){st_write(AtlLines1[[i]], dsn=paste0(names(AtlLines1)[i], "1min.shp"))}


#filter based on error estimates to be figured out later

############look for unrealistic speed##################################
#3min
AtlList <- lapply(AtlList, function (x)
  {x$speed_in <-  atl_get_speed(data=x,x="x.est", y="y.est", 
                                  time="UnixTime", type=c('in'))
  ;x})
AtlList <- lapply(AtlList, function (x)
{x$speed_out <-  atl_get_speed(data=x,x="x.est", y="y.est", 
                              time="UnixTime", type=c('out'))
;x})

par(mfrow=c(1,2))
lapply(AtlList, function(x) hist(x$speed_in))
lapply(AtlList, function(x) hist(x$speed_out))

#1min
AtlList1 <- lapply(AtlList1, function (x)
{x$speed_in <-  atl_get_speed(data=x,x="x.est", y="y.est", 
                              time="UnixTime", type=c('in'))
;x})
AtlList1 <- lapply(AtlList1, function (x)
{x$speed_out <-  atl_get_speed(data=x,x="x.est", y="y.est", 
                               time="UnixTime", type=c('out'))
;x})

par(mfrow=c(1,1))
lapply(AtlList1, function(x) hist(x$speed_in))
lapply(AtlList1, function(x) hist(x$speed_out))

#########################look at angle distribution#################################
#3min
AtlList <- lapply(AtlList, function (x)
{x$angle <-  atl_turning_angle(data=x,x="x.est", y="y.est", 
                               time="UnixTime")
;x})
par(mfrow=c(1,1))
lapply(AtlList, function(x) hist(x$angle))

#1min
AtlList1 <- lapply(AtlList1, function (x)
{x$angle <-  atl_turning_angle(data=x,x="x.est", y="y.est", 
                               time="UnixTime")
;x})
par(mfrow=c(1,1))
lapply(AtlList, function(x) hist(x$angle))


###################################################################################
#filter based on unrealistic speed and distance: keep long straight, short tortuous
###################################################################################
####3min####
AtlList <- lapply(AtlList, function (x) atl_filter_covariates(data=x,
                                                              filters= c(
                                                                "(speed_in < 3 & speed_out < 3) | angle < 70"
                                                              )))
#filter "spikes" or reflections
?atl_remove_reflections
AtlList <- lapply(AtlList, function (x) atl_remove_reflections (data = x, x = "x.est",
                                                                y= "y.est",time = "UnixTime",
                                                                point_angle_cutoff = 45,
                                                                reflection_speed_cutoff = 20))
###1min####
AtlList1 <- lapply(AtlList1, function (x) atl_filter_covariates(data=x,
                                                              filters= c(
                                                                "(speed_in < 12 & speed_out < 12) | angle < 80"
                                                              )))
#filter "spikes" or reflections, do before filter covariates if run
?atl_remove_reflections
AtlList1 <- lapply(AtlList1, function (x) atl_remove_reflections (data = x, x = "x.est",
                                                                y= "y.est",time = "UnixTime",
                                                                point_angle_cutoff = 70,
                                                                reflection_speed_cutoff = 12))                                                               
#######################################################################                                                                                                                                ))
#median smoothing, play with moving window value must be an odd number,
#need not be assigned to new object, modifies in place
########################################################################
#3min
AtlSmooth <- lapply(AtlList, function (x) atl_median_smooth(data = x, x = "x.est",
                                                                y= "y.est",time = "UnixTime",
                                                                moving_window=5))

SmoothLines <- lapply(AtlSmooth, function (x) track(x, x = x$x.est, y = x$y.est, .t=x$DateTime,  
                                           crs=CRS("+init=EPSG:32613"), all_cols=TRUE))
SmoothLines <- lapply(SmoothLines, function (x) as_sf_lines(x))
for (i in 1:length(SmoothLines)){st_write(SmoothLines[[i]], dsn=paste0(names(SmoothLines)[i], "3min.FilterTAS.Smooth.shp"))}

#1min
AtlSmooth1 <- lapply(AtlList1, function (x) atl_median_smooth(data = x, x = "x.est",
                                                            y= "y.est",time = "UnixTime",
                                                            moving_window=5))

SmoothLines1 <- lapply(AtlSmooth1, function (x) track(x, x = x$x.est, y = x$y.est, .t=x$DateTime,  
                                                    crs=CRS("+init=EPSG:32613"), all_cols=TRUE))
SmoothLines1 <- lapply(SmoothLines1, function (x) as_sf_lines(x))
for (i in 1:length(SmoothLines1)){st_write(SmoothLines1[[i]], dsn=paste0(names(SmoothLines1)[i], "1min.FilterTAS.NoSpikes.Smooth.shp"))}

#############################################################
#thin data, do with smoothed data only!#########################
###########################################################
thinned <- lapply(AtlListSmooth, function (x) atl_thin_data(data = x, interval=60,
                                                          id_columns= c("TagId"),
                                                          method= "aggregate or subsampling"))

                                                          
                                                          
#################################################################################                                                                                                                    method = "aggregate or subsample"))
#calculate residence patches
#################################################################################
patches <- lapply(AtlListSmooth, function (x) atl_res_patch(data = x,
                                                            buffer_radius= 10,
                                                            lim_spat_indep = 100,
                                                            lim_time_indep = 30,
                                                            min_fixes = 3,
                                                            summary_variables = c("speed"),
                                                            summary_functions = c('mean', "sd")))

####################################scrap########################
What <- AtlTest[AtlTest$TagId == "072D524C",]
What$speed_out <-  atl_get_speed(data=What,x="x.est", y="y.est", 
                               time="UnixTime", type=c('out'))

