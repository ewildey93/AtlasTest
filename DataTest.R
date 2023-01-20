library(sf)
library(ggplot2)
library(terra)
library(spatialEco)
setwd("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools")

ErrorData <- read.csv("./Calibration1.20/LocError_Dataset.csv")

xestSF <- st_as_sf(x,coords = c(2,3),crs=st_crs(4326))
xestSF <- st_transform(xestSF, crs=st_crs(32613))
x$EstUTMx <- st_coordinates(xestSF)[,1]
x$EstUTMy <- st_coordinates(xestSF)[,2]

xtagSF <- st_as_sf(x,coords = c(12,11),crs=st_crs(4326))
xtagSF <- st_transform(xtagSF, crs=st_crs(32613))
x$TagUTMx <- st_coordinates(xtagSF)[,1]
x$TagUTMy <- st_coordinates(xtagSF)[,2]

x$diffUTMx <- x$TagUTMx-x$EstUTMx
x$diffUTMy <- x$TagUTMy-x$EstUTMy

ggplot(x, aes(diffUTMx, diffUTMy), scale="globalminmax") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  theme_minimal()




#glmm
elev <- rast("C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers/sljm_usgs_1m.tif")
tpi <- tpi(elev, 250, 'circle')
