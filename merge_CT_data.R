install.packages("rgdal", repos = "http://www.gdal.org", type="source")

library(maptools)
library(rgdal)
library(reshape2)

load("CT_array_centroids.RData")

### Below is just an example of how to format a dataset prior to merging it - 
### ArcINFO is picky about the number of characters in the variable names
# Both datasets must be in wide format:
other_data <- read.csv("../gfcanalysis_scripts/summ_stats.csv")
other_data <- other_data[other_data$threshold == "75%" & other_data$aoi == "PA_0", ]
other_data <- dcast(other_data, sitecode ~ variable)
# And ArcGIS requires short variable names (less than 10 characters)
names(other_data) <- abbreviate(names(other_data), 10)
names(other_data) <- make.names(names(other_data), unique=TRUE)
names(other_data) <- gsub("[.]+", "_", names(other_data))
stopifnot(nrow(trap_centroids) == nrow(other_data))

# Merge "other_data" here based on sitecode:
trap_centroids@data <- merge(trap_centroids@data, other_data, by="sitecode")

# Use rgdal to output a shapefile for Kellee
writeOGR(trap_centroids, ".", "CT_array_centroids_merged", driver="ESRI Shapefile")
