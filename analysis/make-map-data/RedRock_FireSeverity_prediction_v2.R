library("raster")
library("rgdal")
library("dismo")
library("sp")
library("ggplot2")

#The following code extracts x,y coords, pre-fire NBR, RdNBR_BL, slope, aspect, elevation, heatload values from each cell within the Red Rock fire

setwd(here::here("analysis", "make-map-data"))

#importing RdNBR
RdNBR = raster("RedRockFire_Indices/RR_RdNBR.tif")
#importing shapefile for fire perimeter
perim <- readOGR(dsn = "RedRockFire_Indices/", layer = "RR_firePerimeter")
#physical vars
elev = raster("RR/RR_DEM_v2.tif")
slp = raster("RR/RR_Slope.tif")
asp = raster("RR/RR_aspect.tif")
#prefire NBR
preNBR = raster("RedRockFire_Indices/RR_prefire_NBR.tif")

#masking fire perimeter
RdNBR.m  = mask(RdNBR,perim);plot(RdNBR.m)

#creating a point file at center of each raster cell
RdNBR.pt = rasterToPoints(RdNBR.m,spatial = TRUE);plot(RdNBR.pt)
RdNBR.pt.ll = spTransform(RdNBR.pt, CRS("+proj=longlat +datum=WGS84"))


proj <- CRS("+proj=longlat +datum=WGS84")

#matching all projections, files from MTBS (prefire, rdnbr, perim) all have same projection
elev.p = projectRaster(elev,crs = proj)
slp.p = projectRaster(slp,crs = proj)
asp.p = projectRaster(asp,crs = proj)

#creating dataframe with values from RdNBR.pt
data <- data.frame(coordinates(RdNBR.pt),
                   coordinates(RdNBR.pt.ll),
                   extract(preNBR,RdNBR.pt, method = "simple"),
                   extract(RdNBR.m, RdNBR.pt, method = "bilinear"),
                   extract(slp,RdNBR.pt, method = "simple"),
                   extract(asp,RdNBR.pt, method = "simple"),
                   extract(elev,RdNBR.pt, method = "simple"))
#renaming columns
colnames(data) <- c("x", "y","long","lat","preNBR", "RdNBR_BL","Slope","Aspect","Elevation");head(data)

###computing heat.load
##converting to radians
#lat
data$lat.r = data$lat*pi/180
#slope
data$Slope.r = data$Slope*pi/180
#aspect
data$Aspect.r = data$Aspect*pi/180
#folded aspect
data$fold.aspect = abs(pi-abs(data$Aspect.r-pi*5/4))
#heatload_in
data$heatload.in =
  -1.467+1.582*cos(data$lat.r)*cos(data$Slope.r) -
  1.5*cos(data$fold.aspect)*sin(data$Slope.r)*sin(data$lat.r) -
  0.262*sin(data$lat.r)*sin(data$Slope.r) +
  0.607*sin(data$fold.aspect)*sin(data$Slope.r)
data$heatload = exp(data$heatload.in)

#exporting df
dir.create(here::here("data/generated"), showWarnings = FALSE)
saveRDS(data, here::here("data/generated/RedRocks_map_data.rds"), compress = TRUE)

setwd(here::here())
