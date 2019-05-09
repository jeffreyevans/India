###############################################
#### Sample design for Kahna Pench, India  ####
####   bioacoustic study                   ####
###############################################
library(sp)
library(raster)
library(spatialEco)
library(rgdal)

# Working (path) and data directories (dpath)
path = "C:/evans/India/kahna_pench/data"
dpath = "C:/evans/India/kahna_pench/raw_data/" 
setwd(path)

utm.proj <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
out.shape = "stratification_v2"  # name of resulting shapefile
p = 600                       # minimum sample distance
min.dist = 200                # minimum distance for lag sample
max.dist = 600                # maximum distance for lag sample
sn = 4                        # Number of random samples per strata
sr = 4                        # Number of strata replicates

# function for standard error
std.err <- function(x, na.rm = TRUE) {
  if(na.rm) { 
    x <- x[!is.na(x)] 
    return( sd(x) / sqrt(length(x)) )
  } else {
    return( NA )
  }
}

###############################################
# Read, crop and mask stratification data
###############################################
# Stratification (strat.tif) 
# * strata used in stratification
#   strat  disturbance forest go/no
#   1      1           1      0
#  *2      1           0      0
#   3      2           1      0
#   4      3           1      0
#  *5      3           0      0
#   6      4           1      0
#  *7      4           0      0
#  *8      3           0      1
#   9      3           1      1
#   10     4           1      1
#  *11     4           0      1
#  *12     2           0      1
#   13     2           1      1
#  *14     2           0      0
#  *15     1           0      1
#   16     1           1      1

# Read stratification raster and table
r <- raster("strat.tif")
  dat <- read.csv("strat_table.csv")
    dat <- dat[,-2]

# Index "forested" strata and remove 
#   strata raster	
fidx <- dat[,"value"][which(dat$forest == 1)]
  r[which(r[] %in% fidx)] <- NA

# Number of unique strata
cat("Number of unique strata", length(unique(r[])) -1, "\n")

# Read forest/non-forest raster 
f <- raster(paste0(dpath,"forest_ANYeq1.tif"))
 f <- mask(f,crop(f, extent(r), snap="in"),r)

# Read disturbance raster 
d <- raster(paste0(dpath,"di3_sum3_pta_mahmp_1ki.tif"))
  d <- mask(crop(d, extent(r)),r)
    #d <- d / 1000

###############################################
# sample non-forest in go and no-go     
#   using disturbance as sample weights 
###############################################

strata.ids <- as.vector(na.omit(unique(r[]))) 

# Power-based minimal sample size for each strata
n = vector()
  for( i in strata.ids ) {
    n <- append(n, round((( 2 * sd(d[which(r[] == i)], 
                na.rm = TRUE))^2 / 
				std.err(d[which(r[] == i)])) ,0))   
    }
  names(n) <- strata.ids
    n <- round(n / 20, 0)
# Draw stratified random samples
strata <- data.frame(matrix(vector(), 0, 3, 
                     dimnames=list(c(), 
				     c("x", "y", "strata"))))
  for( i in strata.ids ) {
      idx <- sample(which(r[] == i), 
	                size = (n[which(names(n) == i)]), 
	                replace = FALSE)
          cat("Pulling", length(idx), "random samples for strata", i, "\n")
    strata <- rbind(strata, data.frame(xyFromCell(r, idx), 
				    strata = i) )
  }
# Coerce results to sp SpatialPointsDataFrame object
coordinates(strata) <- ~x+y  

# Count number of obs in each strata
  for( i in strata.ids ) {
    cat("Obs in strata: ", i, 
	    length(strata[strata$strata == i,]$strata), "\n")
  }

# N per strata 4; reps per N, per strata 4  
strata.n4r4 <- stratified.random(strata, strata = "strata", 
                                 n = sn, reps = sr)
strata.sub <- strata.n4r4[strata.n4r4$REP == 1,] 
  rad.list <- list()
    for(i in 1:nrow(strata.sub)) {							
      d1 <- sample.annulus(strata.sub[i,], r1=min.dist, r2=(min.dist+min.dist), 
                           n = 1, type = "random",iter=50)
  	  d1$rmin <- min.dist
  	  d1$rmax <- (min.dist + min.dist)
      d2 <- sample.annulus(strata.sub[i,], r1=(min.dist + min.dist), r2=max.dist, 
                           n = 1, type = "random",iter=100)
   	  d2$rmin <- (min.dist + min.dist)
  	  d2$rmax <- max.dist  
  	rad.list[[i]] <- rbind(d1,d2) 
  	}

# Format final stratification SpatialPointsDataFrame
cluster.plots <- do.call("rbind", rad.list)
strata.sub@data <- data.frame(SID=as.numeric(row.names(strata.sub@data)),
                              rmin=0, rmax=0,s=strata.sub$strata)						  
cluster.plots <- merge(cluster.plots, strata.sub@data[,c("SID","s")], by="SID")
  cluster.plots@data <- cluster.plots@data[,-2] 
  						   
strata.final <- rbind(strata.sub,cluster.plots)
  strata.final$strata <- as.numeric(as.character(strata.final$s)) 
  strata.final@data <- strata.final@data[,-which(names(strata.final) == "s")]
  strata.final <- strata.final[order(strata.final$strata, strata.final$SID,strata.final$rmin, 
                               decreasing = c(FALSE,FALSE,FALSE)),] 

# Write stratification shapefile
proj4string(strata.final) <- utm.proj
  writeOGR(strata.final, getwd(), out.shape, driver="ESRI Shapefile",
           check_exists=TRUE, overwrite_layer=TRUE)

###############################################
# Test effect size
###############################################
d <- raster(paste0(dpath,"di3_sum3_pta_mahmp_1ki.tif"))
strata.final$disturb <- extract(d, strata.final)/1000

effect.size <- list()
cbd <- combn(unique(strata.final$strata), 2)
  for(i in 1:ncol(cbd)) {
    es.test <- data.frame(strata = strata.final[strata.final$strata == cbd[,i][1] |
                                                strata.final$strata == cbd[,i][2],]$strata,  
                          y = strata.final[strata.final$strata == cbd[,i][1] |
                                           strata.final$strata == cbd[,i][2],]$disturb)
  effect.size[[ paste0(cbd[,i][1], "_", cbd[,i][2]) ]] <- psych::cohen.d(es.test, "strata")$cohen.d
  }
( effect.size <- data.frame(strata = names(effect.size), 
                          do.call("rbind", effect.size)) )
write.csv(effect.size, "effect_size.csv", row.names=FALSE)

effect.size.plot <- list()
strata.plot <- strata.final[strata.final$rmin ==0,] 
cbd <- combn(unique(strata.plot$strata), 2)
  for(i in 1:ncol(cbd)) {
    es.test <- data.frame(strata = strata.plot[strata.plot$strata == cbd[,i][1] |
                                                strata.plot$strata == cbd[,i][2],]$strata,  
                          y = strata.final[strata.plot$strata == cbd[,i][1] |
                                           strata.plot$strata == cbd[,i][2],]$disturb)
  effect.size.plot[[ paste0(cbd[,i][1], "_", cbd[,i][2]) ]] <- psych::cohen.d(es.test, "strata")$cohen.d
  }
( effect.size.plot <- data.frame(strata = names(effect.size.plot), 
                          do.call("rbind", effect.size.plot)) )
write.csv(effect.size.plot, "effect_size_plot.csv", row.names=FALSE)

# Save R workspace
save.image(paste0(path, "/stratification_v2.RData"))