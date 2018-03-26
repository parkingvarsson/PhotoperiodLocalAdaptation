#! /usr/bin/Rscript --no-save --no-restore


###This script is used to estimate the physical distance between populations
###step1: calculate the mean latitude and longtitude for all 12 populatiosn in this study
library(taRifx)

setwd("/gulo/proj_nobackup/b2011141/local_adaptation_paper/asp201_94/env/raw_data")
#swasp=read.table("SwAsp.Clone.Origin.data.pop.txt",header=T)

#swasp_lat=by(swasp$Tree.latitude,swasp$pop,mean)
#swasp_long=by(swasp$Tree.longitude,swasp$pop,mean)
#lat_mean=as.data.frame(swasp_lat)
#long_mean=as.data.frame(swasp_long)
#names(lat_mean)=c("Pop","Latitude")
#names(long_mean)=c("Pop","Longitude")
#
#distance=merge(lat_mean,long_mean)
#
#write.table(distance,file="SwAsp.94samples.lat.long.txt",sep="\t",quote=F,row.names=F,col.names=T)


###Step2, calculate physical distance for all samples
### Geodist is a function wich calculates the geographic distance between two
### points defined by their coordinates
### As input the function needs a text file with three columns. Names of points
### are in the first column, latitudes in the second and longitudes in the third.

geod <- function(la1, lo1, la2, lo2) {
  geo <- 2 * 6378.2 * asin (sqrt( ((sin((pi / 180) * ((la1-la2)/2)))^2) + (cos((pi / 180) *(la1)) * cos((pi / 180) *(la2)) * (sin((pi / 180) *((lo1-lo2)/2)))^2)))
  return(geo)
}

geod.simul <- function(la1, lo1, la2, lo2) {
  geo <- sqrt((la1-la2)^2 + (lo1-lo2)^2)
  return(geo)
}


distmat <- function (coord){
  dmat <- matrix(NA, dim(coord)[1], dim(coord)[1])
  for(i in 1:dim(coord)[1]){
    for (j in 1:dim(coord)[1]){
      dmat[i, j] <- geod(coord[i, 1], coord[i, 2], coord[j, 1], coord[j, 2])
    }
  }
  return(dmat)
}

distmat.simul <- function (coord){
  dmat <- matrix(NA, dim(coord)[1], dim(coord)[1])
  for(i in 1:dim(coord)[1]){
    for (j in 1:dim(coord)[1]){
      dmat[i, j] <- geod.simul(coord[i, 1], coord[i, 2], coord[j, 1], coord[j, 2])
    }
  }
  return(dmat)
}

##here the input file is the file created in step1 "SwAsp.94samples.lat.long.txt"
Geodist <- function(filename) {
  
  coordlist <- read.table(file=filename,header=T)
  coords <- as.matrix(coordlist[ , 2:3])
  pointnames <- as.vector(coordlist[ , 1])
  
  distances <- distmat(coords)
  
  matr <- cbind(pointnames, distances)
  matr <- rbind(c("x", pointnames), matr)
  write.table(matr, file="geo-distances.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}



