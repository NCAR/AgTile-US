library(raster)
library(RNetCDF)

#read tiff file

data<-raster("path to AgTile-US.tif")

# read netcdf

data <- var.get.nc(open.nc("path to AgTile-US.nc"),"AgTile-US")

