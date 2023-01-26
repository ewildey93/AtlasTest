library(sf)
library(ggplot2)
library(terra)
library(spatialEco)
library(tmap)
setwd("C:/Users/eliwi/OneDrive/Documents/R/data_tools/data_tools")

x <- read.csv("./Calibration1.20/Trilateration.Test.Data_Filter.Distance.187.5_Results.csv")
str(x)
ErrorData <- read.csv("./Calibration1.20/LocError_Dataset.csv")

xestSF <- st_as_sf(x,coords = c(3,4),crs=st_crs(4326))
xestSF <- st_transform(xestSF, crs=st_crs(32613))
x$EstUTMx <- st_coordinates(xestSF)[,1]
x$EstUTMy <- st_coordinates(xestSF)[,2]

xtagSF <- st_as_sf(x,coords = c(9,10),crs=st_crs(4326))
xtagSF <- st_transform(xtagSF, crs=st_crs(32613))
x$TagUTMx <- st_coordinates(xtagSF)[,1]
x$TagUTMy <- st_coordinates(xtagSF)[,2]

x$CIX <- x$x.ci.upper-x$x.ci.lower
x$CIY <- x$y.ci.upper-x$y.ci.lower

x$diffUTMx <- x$TagUTMx-x$EstUTMx
x$diffUTMy <- x$TagUTMy-x$EstUTMy

ggplot(x, aes(diffUTMx, diffUTMy), scale="globalminmax") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  theme_minimal()


x$dist <- terra::distance(as.matrix(x[,16:17]), as.matrix(x[,14:15]),pairwise=TRUE, lonlat=FALSE)
str(x)
ggplot(x, aes(CIX,diff.dist )) +
  geom_point()

ggplot(x, aes(CIY,diff.dist )) +
  geom_point()


#glmm
elev <- rast("C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers/sljm_usgs_1m.tif")
elevcrop <- crop(elev, ext)

#https://bookdown.org/hhwagner1/LandGenCourse_book/bonus_2a.html
Map3 <- tm_shape(elev)# + 
  #tm_raster(style="cat", palette=RAT$color, labels=RAT$attribute,
  #         title="Land cover") +
  #tm_layout(legend.outside=TRUE, legend.outside.position="right") +
  #tm_grid(lines=FALSE)
Map3
Map4 <- Map3 + tm_shape(Sites.sf) +
  tm_symbols(size=0.4, col="yellow", border.col="red") +
  tm_compass() + tm_scale_bar(bg.color="lightgray", bg.alpha=0.5)
Map4



tpi <- tpi(elev, 250, 'circle')

plot(st_geometry(xestSF))
plot(elev)
ext <- draw(x='extent') #this is handy but not reproducible unless you save the values
ext
ext_save <- c(ext[1:4])
nlcd_crop2 <- crop(nlcd, ext)
