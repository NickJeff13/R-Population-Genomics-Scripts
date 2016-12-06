#Get least cost distances among sampling sites using marmap
#Set working directory
setwd("C:/Users/Nick/Desktop/Postdoc/Atlantic Salmon Project/Paper 1 - Clines/sPCA")

#libraries -------------
library(gdistance)
library(maps)
library(mapdata)
library(maptools)
library(marmap)
library(PBSmapping)
library(ggplot2)
library(dplyr)
library(data.table)
library(gtable)
library(vegan)

#source mapping functions
# source("c:/Users/rystanley/OneDrive/PostDoc/DFO/Cline Analysis/Cartesian Analysis/ScaleBar.R")


#Make colours
blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))


#data ---------

#Load workspace
#Load coordinates as degrees Lat and Long
SpeciesCoordinates <- read.csv("Salmon_Coords.csv")




## Set map limits------
Lat.lim=c(42,62)
Long.lim=c(-72,-52)

#Download bathymetry layer for salmon sampling range - this gets the bathymetry for your area to then calculate lc dists
#Setting the resolution to a lower number (e.g. 1 is lowest) takes longer but grabs more detail
getNOAA.bathy(lon1 = Long.lim[1], lon2 = Long.lim[2], lat1 = Lat.lim[1], lat2 = Lat.lim[2],
              resolution = 1,keep=TRUE) -> SalmonDepths

#plot the bathy map to see how it looks
plot(SalmonDepths,image = TRUE, land = TRUE, lwd = 0.03,
     bpal = list(c(0, max(SalmonDepths), greys),
                 c(min(SalmonDepths), 0, blues)))
plot(SalmonDepths, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
points(SpeciesCoordinates$AdjustLong,  SpeciesCoordinates$AdjustLat,pch=19,cex=1,col="black") #Add points
#Add text labels
text(x=SpeciesCoordinates$AdjustLong,y=SpeciesCoordinates$AdjustLat,labels=SpeciesCoordinates$Code, cex=0.6,pos=3,offset = 0.5,font = 2)
#### Calculate distance (km) between stations for each species -------------

#Transition object - this can take a few minutes and is used as the trans object in the lc.dist function below
trans1 <- trans.mat(SalmonDepths,min.depth = -10) 

sites1<-SpeciesCoordinates[,5:6]
rownames(sites1)<-SpeciesCoordinates$Code

#calculate least cost distances - just use res="dist" to get a distance matrix
dist.salmon <- lc.dist(trans1, 
                       sites1, 
                       res="dist")

#res = "path" will give you a list of all distances between each pairwise site but can take a long time
salmon.lines<- lc.dist(trans1, 
                       sites1, 
                       res="path")