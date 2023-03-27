library(raster)
library(sp)
library(rgdal)
library(sf)
library(ggplot2)
library(geosphere)
library(rgeos)
library(dplyr)
library(tidyr)
library(RPostgres)
library(tidyverse)
library(raster)
library(nlstools)
setwd("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools")
source("functions/data_manager.R")
source("functions/localization.R")
source("C:/Users/eliwi/OneDrive/Documents/R/EcolEvol.Manuscript_Optimizing.Trilateration/4_Functions_RSS.Based.Localizations.R")


#Gather beep data
#Stations: 2300631E77E9 61AA8ADD12C0
db_name <- "Salida"
conn <- dbConnect(RPostgres::Postgres(), dbname=db_name, password= 'arsenal4')
flights_db <- tbl(conn, "raw")
beep <- flights_db %>% 
  filter(time > '2022-01-06' & time < '2022-01-15') %>%
  collect()
beep2 <- flights_db %>% 
  filter(time > '2022-01-06' & time < '2022-01-15' & station_id == "2300631E77E9") %>%
  collect()
beep_data <- rbind(beep, beep2)
colnames(beep_data)[c(3,4,5,6,7,8)] <- c('RadioId','TagId','NodeId','TagRSSI','Validated','Time')
str(beep_data)
attributes(beep_data$Time)$tzone <- "America/Denver"
str(beep_data)
beep_data <- distinct(beep_data, TagId, TagRSSI, NodeId, Time, .keep_all=TRUE )
beep_data$NodeId <- toupper(beep_data$NodeId)

beep_data <- readRDS("C:/Users/eliwi/OneDrive/Documents/R/AtlasTest/BeepMerge8.4.rds")
str(beep_data)
attributes(beep_data$Time)
attributes(beep_data$Time)$tzone <- "America/Denver"

#bring in node file
nodes <- read.csv("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools/Calibration1.20/Nodes1.20.csv", as.is=TRUE, na.strings=c("NA", ""), strip.white=TRUE) 
colnames(nodes)[c(1,3,4)] <- c('NodeId', 'lat','lng')
nodesSF <- st_as_sf(nodes,coords = c(4,3),crs=st_crs(4326))
nodesSF <- st_transform(nodesSF, crs=st_crs(32613))
nodes$NodeUTMx <- st_coordinates(nodesSF)[,1]
nodes$NodeUTMy <- st_coordinates(nodesSF)[,2]
nodes$NodeId <- toupper(nodes$NodeId)


#format calibration DF, check for scientific notation 52551E33 33071E07
calibration <- read.csv("./test8.4/Calibration8.4.csv")
table(grepl(pattern='e+', calibration$TagId))
calibration <- calibration[,1:8]
calibration$start <- as.POSIXct(calibration$start, format="%m/%d/%Y %H:%M", tz="America/Denver")
calibration$end <- as.POSIXct(calibration$end,format="%m/%d/%Y %H:%M", tz="America/Denver")
calibration$start <- calibration$start + 60
calibration$end <- calibration$end - 60

calibrationsf <- st_as_sf(x = calibration,coords = c(7,6), crs=st_crs(4326))
calibrationsf <- st_transform(calibrationsf, crs=CRS("+init=epsg:32613"))
st_coordinates(calibrationsf)[,1]
calibration$TestUTMx <- st_coordinates(calibrationsf)[,1]
calibration$TestUTMy <- st_coordinates(calibrationsf)[,2]
colnames(calibration)[2] <- "TestId"



calibrate <- function(beep_data, calibration, nodes) {
#filter beep_data to those within the calibration time frame
beep_data <- beep_data[beep_data$TagId %in% calibration$TagId,]
dt1 <- data.table(beep_data, start=beep_data$Time, end=beep_data$Time)
dt2 <- data.table(calibration)
setkey(dt2, TagId, start, end)
indx <- foverlaps(dt1, dt2, type='within')
beep_data <- indx[!is.na(indx$start),]
colnames(beep_data)[3] <- "TestId"

## Calculate average RSS value for each unique test and node 
summary.test.tags <- beep_data %>%
  dplyr::group_by(NodeId, TestId) %>%
  dplyr::summarise(avgRSS = mean(TagRSSI),
                   sdRSS = sd(TagRSSI),
                   n.det = n())


## Calculate euclidean distance between each test location and node

# Calculate distance between nodes and test locations
dst <- raster::pointDistance(calibration[,c("TestUTMx", "TestUTMy")], nodes[,c("NodeUTMx", "NodeUTMy")], lonlat = F, allpairs = T)

# Make matrix into a dataframe
dist_df <- data.frame(dst, row.names = calibration$TestId)
colnames(dist_df) <- nodes$NodeId
dist_df$TestId <- rownames(dist_df)
dist_df$TestId <- as.integer(dist_df$TestId)

# rearrange data
dist.gather <- dist_df %>%
  tidyr::gather(key = "NodeId", value = "distance", -TestId)


## Combine distances with summary data 
summary.dist <- summary.test.tags %>%
  dplyr::left_join(dist.gather) 

# Add UTMs of nodes and test locations to dataset
summary.dist <- summary.dist %>%
  dplyr::left_join(nodes[, c("NodeId", "NodeUTMx", "NodeUTMy")]) %>%
  dplyr::left_join(calibration[, c("TestId", "TestUTMx", "TestUTMy")]) 


## save file
write.csv(summary.dist, paste0(outpath, TEST.TYPE, "_Dataset.csv"),
          row.names = F)


## Add to R environment
return(summary.dist)


}
outpath <- './test8.4'
Calib8.4 <- calibrate(beep_data, calibration, nodes)
###########################a,S,K all nodes#################################
CalibDS <- read.csv("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools/Calibration1.20/Calibration_Dataset.csv")
CalibDS <- summary.dist
# Preliminary Model
exp.mod <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc), data = Calib8.4)
a <- coef(exp.mod)[["R0"]]
S <- exp(coef(exp.mod)[["lrc"]])
K <- coef(exp.mod)[["Asym"]]
relation <- relate(a,S,K)

# Final Model
nls.mod <- nls(avgRSS ~ a * exp(-S * distance) + K, start = list(a = a, S = S, K= K), 
               data = summary.dist)
coef(nls.mod)
a <- coef(nls.mod)[["a"]]
S <- coef(nls.mod)[["S"]]
K <- coef(nls.mod)[["K"]]
######################a,s,K by node########################################
nodeQ <- unique(CalibDS$NodeId)
noderelates <- data.frame(NodeId=character(),a=numeric(),S=numeric(),K=numeric(), N=numeric())
for (i in 1:length(nodeQ)){ 
  tryCatch({
    ID <- nodeQ[i]
    print(ID)
    beeps <- filter(CalibDS, NodeId == ID)
    exp.mod <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc), data = beeps)
    intercept <- coef(exp.mod)[["R0"]]
    rate <- exp(coef(exp.mod)[["lrc"]])
    asymp <- coef(exp.mod)[["Asym"]]
    
    nls.mod <- nls(avgRSS ~ a * exp(-S * distance) + K, start = list(a = intercept, S = rate, K= asymp), 
                   data = beeps)
    #a,S,K
    x <- data.frame(NodeId=ID, a=coef(nls.mod)[["a"]], S=(coef(nls.mod)[["S"]]),
                    K=coef(nls.mod)[["K"]], N=nrow(beeps))
    
    #populate dataframe with a,S,K for each node RSS~distance relationship
    noderelates <- rbind(noderelates, x)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

noderelates <- left_join(nodes[,c(1,3,4)], noderelates, by='NodeId')
#taken from a,S,K estimates for all nodes!!!! make sure you're referencing the right parameters
noderelates$a[is.na(noderelates$a)] <- as.numeric(coef(nls.mod)[1])
noderelates$S[is.na(noderelates$S)] <- as.numeric(coef(nls.mod)[2])
noderelates$K[is.na(noderelates$K)] <- as.numeric(coef(nls.mod)[3])
noderelates$N[is.na(noderelates$N)] <- 0
saveRDS(noderelates, "C:/Users/eliwi/OneDrive/Documents/R/AtlasTest/NodeRelates.rds")
############################################################################################################

#For calibration data
# Function to estimate the distance of each test signal from the node based on the RSS values detected
#need to modify for parameters by node
combined.data <- estimate.distance(CalibDS, noderelate = noderelates) #combined.data=CalibDS, summary.dist
combined.data.singleASK <- estimate.distance(CalibDS)
#for individual node a,S,K
estimate.distance <- function(x, noderelate=NULL) {
  
  # supress warnings
  options(warn = -1)
  
  
  if(!is.null(noderelate)) {
  #combine node relate parameters with combined.data frame
  x <- left_join(x, noderelates, by='NodeId')
  # Calculate estimated distance based on RSSI~Distance relationship and indicate simulation round
  combined.data <- x %>%
    dplyr::mutate(e.dist = (log(avgRSS - K) - log(a)) / -S)
  }
  else {
    combined.data <- x %>%
      dplyr::mutate(e.dist = (log(avgRSS - K) - log(a)) / -S)
  }
  
  # Remove rows with NAs
  combined.data <- combined.data[complete.cases(combined.data),]
  
  # Change negative distances to 10 (rss values > intercept of exponential curve and thus negative) - indicates very close to the node
  combined.data <- combined.data %>%
    dplyr::mutate(e.dist = dplyr::case_when(e.dist < 0 ~ 10,
                                            e.dist >=0 ~ e.dist))
  
  # Save File to outpath
  write.csv(combined.data, paste0(outpath, "LocError_Dataset.csv"),
            row.names = F)
  
  return(combined.data)
  
}



########################Trilaterate with distance filters#####################
DIST.FILTER <- c(187.5, 300)

#outpath to dump results of trilateration
outpath <- "C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools/Calibration1.20/"


# Calculate error of location estimates of each test location when Distance filters are applied prior to trilateration 
Dist.filters <- trilateration.TestData.Distance.Filter(combined.data)
Dist.filters.singleASK <- trilateration.TestData.Distance.Filter(combined.data.singleASK)



Dist.filters <- read.csv("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools/Calibration1.20/Trilateration.Test.Data_Filter.Distance.187.5_Results.ASKSbyNode.csv")
Dist.filters.singleASK <- read.csv("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools/Calibration1.20/Trilateration.Test.Data_Filter.Distance.187.5_Results.csv")



######################scrap###################################################
Node3398D3 <- filter(CalibDS, NodeId == "3398D3")
exp.mod3398D3 <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc), data = Node3398D3)
coef(exp.mod3398D3)
a1 <- coef(exp.mod3398D3)[["R0"]]
s1 <- exp(coef(exp.mod3398D3)[["lrc"]])
k1 <- coef(exp.mod3398D3)[["Asym"]]

nls.mod3398D3 <- nls(avgRSS ~ a * exp(-S * distance) + K, start = list(a = a1, S = s1, K= k1), 
               data = Node3398D3)
coef(nls.mod3398D3)
