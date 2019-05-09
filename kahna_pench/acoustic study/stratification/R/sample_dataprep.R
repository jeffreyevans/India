library(sp)
library(raster)
library(spatialEco)
library(rgdal)
library(rgeos)

setwd("C:/evans/India/kahna_pench/raw_data")

#e <- readOGR(getwd(), "KP_box")
e <- extent(890239.9, 1129212, 2362567, 2501749) 

# Go-NoGo, disturbance, forest 
disturbance <- crop(raster("di3_sum3_pta_mahmp_Q4.tif"),e)
forest <- crop(raster("forest_ANYeq1.tif"),e)
go <- crop(raster("s2_rem1_gt1km2_01.tif"),e)
  go <- crop(resample(go, disturbance, method="ngb"),e)

r <- stack(go, disturbance, forest)

