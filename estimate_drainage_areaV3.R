# This R code can be used for estimating gridded 30-m 
# tile drainage area from USDA county level tile drainage census data.
# Main inputs are county-level census, SSURGO gridded soil drainage data,
# NLCD cropland mask, and SRTM DEM derived slope.
# Any question?, please contact prasanth@ucar.edu
#---------------------------
rm(list=ls())
library(raster)

#---------------------------
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
#---------------------------
setwd("E:/DRAINAGE/SOIL_DATA/CONUS/FINAL/")
outdir<-"E:/DRAINAGE/CDL/out/"

fc<-list.files("./census/split/",pattern = ".shp")


s<-0

for (s in 0:51){ # States loop
  countyData   <-shapefile(paste("./census/split/",s+1,".shp",sep="")) 
  soilData     <-raster(paste("./soil/split/soil",s,".TIF",sep=""),RAT = TRUE)
  lan_covr     <-raster(paste("E:/DRAINAGE/CDL/split/crop",s,".TIF",sep=""))
  slope        <-raster(paste("./slope/split/slope",s,".TIF",sep=""))
  extent(lan_covr)<-extent(soilData)
  res(lan_covr)<-30
  #-----------------------------------------------------
  
  fid<-countyData$GEOID
  
  catdr<- vector(length = length(fid))
  slpdr<- vector(length = length(fid))
  slopd<- vector(length = length(fid))
  rmse <- vector(length = length(fid))
  pars<-data.frame(fid,catdr,slpdr,slopd,rmse)
  
  i<-5
  j<-1
  k<-1
  
  for(i in 1:length(fid)) { # county loop
    
    countymsk <- subset(countyData,GEOID==fid[i])
    sD_crp    <- crop(soilData,countymsk) 
    sd_msk    <- mask(sD_crp  ,countymsk)
    #plot( countymsk )
    print(paste(i,fid[i]))
    if(countymsk$Area >1){
      
      
      lC_crp  <-crop(lan_covr,countymsk) 
      slp_crp <-crop(slope,countymsk)
      
      
      lc_msk<-mask(lC_crp  ,countymsk)
      sl_msk<-mask(slp_crp ,countymsk)
      
      
      lc_msk[!(lc_msk==1)]<-NA
      
      drc<-5
      for (drc in 7:1){
        drn<-sd_msk
        drn[drn <  drc] <-NA
        drn[drn >= drc] <-1
        
        for(slop in 0:20){
          #slop <- 15
          slp<-sl_msk
          slp[slp >  slop] <-NA
          slp[slp <= slop]<-1
          
          drm <- drn*slp*lc_msk
          #plot(drm)
          A_m <- sum(getValues(drm),na.rm = T)*30*30/4046.86
          A_s <- countymsk$Area
          dA  <- A_s - A_m
          print(paste(slop,A_s,round(A_m),round(dA),drc))
          if(dA<0){break} 
        }
        if(dA<0){break} 
      }
      
      #------------------
      if((drc==1)&&(dA>0)&&(A_s>100)&&(slop>=20)){
        drn3<-sd_msk
        drn3[drn3 <  4] <-NA
        drn3[drn3 >= 4] <-1
        drm3 <- drn3*lc_msk
        soil_sum <- sum(getValues(drm3),na.rm=T)*30*30/4046.86
        
        if(soil_sum<countymsk$Area){
          drm3[is.na(drm3)]<-0
          writeRaster(drm3, filename=paste(outdir,s,"_",fid[i],".tif",sep=""), format="GTiff", overwrite=TRUE)
          print(paste(fid[i],i, "done"))
          pars$catdr[i]<- NA
          pars$slpdr[i]<- NA
          pars$rmse[i] <- NA
          pars$slopd[i]<- 20
        }
        next
      }
      if((drc==1)&&(dA>0)&&(A_s<100)){
        next
      }
      
      if(drc<7){
        drn1<-sd_msk
        drn1[drn1 <  (drc+1)] <-NA
        drn1[drn1 >= (drc+1)] <-1
        
        slp<-sl_msk
        slp[slp >  slop] <-NA
        slp[slp <= slop]<-1
        
        drm1<-drn1*slp*lc_msk
        drm1[is.na(drm1)]<-0
      }else{
        drm1<-sd_msk*0
        drm1[is.na(drm1)]<-0
      }
      #writeRaster(drm1, filename="drm1.tif", format="GTiff", overwrite=TRUE)
      
      slv<-seq(slop-5,slop,by=0.1)
      jv   <- vector(length = length(slv))
      rmsv <- vector(length = length(slv))
      dfrmse<-data.frame(jv,rmsv)
      
      for (j in 1:(length(slv))){
        
        slp<-sl_msk
        slp[slp > slv[j]] <-NA
        slp[slp <= slv[j]]<-1
        
        drn2<-sd_msk
        # drn2[drn2 < drc]<-NA
        # drn2[drn2 > (drc+1)] <-NA
        # drn2[drn2 >= drc ] <-1
        
        drn2[!(drn2==drc)]<-NA
        drn2[!(is.na(drn2))]<-1
        
        drm2 <- drn2*slp*lc_msk
        drm2[is.na(drm2)]<-0
        
        drm <- drm1+drm2
        drm[drm>0]<-1
        A_m <- sum(getValues(drm),na.rm = T)*30*30/4046.86
        A_s <- countymsk$Area
        
        dfrmse$jv[j]   <- j
        dfrmse$rmsv[j] <- RMSE(A_m,A_s)
        print(j)
        
      }
      print("slope done")
      rwn<-which(dfrmse$rmsv==min(dfrmse$rmsv))
      
      j<-dfrmse$jv[rwn]
      
      vl<-max(slv[j])
      slp<-sl_msk
      slp[slp > vl] <-NA
      slp[slp <= vl]<-1
      
      drn2<-sd_msk
      drn2[!(drn2==drc)]<-NA
      drn2[!(is.na(drn2))]<-1
      
      drm2 <- drn2*slp*lc_msk
      drm2[is.na(drm2)]<-0
      #writeRaster(drm2, filename="drm2.tif", format="GTiff", overwrite=TRUE)
      
      drm <- drm1+drm2
      drm[drm>0]<-1
      
      drm[drm<1]<-NA
      
      drm[is.na(drm)]<-0
      
      writeRaster(drm, filename=paste(outdir,s,"_",fid[i],".tif",sep=""), format="GTiff", overwrite=TRUE)
      print(paste(fid[i],i, "done"))
      
      pars$catdr[i]<- drc+1
      pars$slpdr[i]<- vl
      pars$rmse[i] <- dfrmse$rmsv[max(j)]
      pars$slopd[i]<- slop
      
    }else{
      drm <- sd_msk
      drm[drm<10]<-NA
      drm[is.na(drm)]<-0
      writeRaster(drm, filename=paste(outdir,s,"_",fid[i],".tif",sep=""), format="GTiff", overwrite=TRUE)
      print(paste(fid[i],i, "done"))
      pars$catdr[i]<- NA
      pars$slpdr[i]<- NA
      pars$rmse[i] <- NA
      pars$slopd[i]<- NA
    }
    
  }
  write.csv(pars,paste(outdir,s,"_",fid[i],".csv",sep=""))
}




