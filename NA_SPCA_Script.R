#Description ------------
## Run a Spatial principal components analysis for salmon using neutral and outlier loci
setwd("/home/ian/Desktop/Nick/AtlanticSalmon/SPCA/")
load("NA_SPCA_Data.RData")
#STEP 1 is to load cartesian data and create a vector which matches the
#genepop file for each species. The order of samples needs to match the coordinates
#so that the coordinates can be assigned to each GENIND object from
#adegenet

#STEP 2 is to create genind objects and assign coordinate data to the other$xy slot

#STEP 3 will run the sPCA for each loci subset. 

#STEP 4 will combine the sPCA results and make a master plot which summarizes each.

#load libraries ----------
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(adegenet)
library(genepopedit)
library(akima)
library(splancs)
library(gridExtra)

inds<-genepop_detective("NA_75OutliersGP_NoAqua.txt",variable = "Inds")
duplicated(inds)
sum(duplicated(inds))
inds3<-genepop_detective("NA_Neutral_1712lociGP.txt",variable = "Inds")
## Step 1 Load coordintes ------------

NACoords <- read.csv("SalmonCartesianCoords.csv",header=T)
head(NACoords)
inds2<-as.vector(NACoords$SampleID)
duplicated(inds2)
sum(duplicated(inds2))
length(intersect(inds3,inds2))
setdiff(inds,inds2)

##Jitter the species coordinates so that a Delaunay triangulation can be applied. Here we use a jitter of 2 km 
NACoords$jCartx <- jitter(NACoords$MDS1,2) 
NACoords$jCarty <- jitter(NACoords$MDS2,2)

##### Step 2 Read in genepop files as GENIND objects and add coordinatse------------

NAoutliers <- adegenet::read.genepop("NA_75OutliersGP_NoAqua.GEN",
                                       ncode=3,quiet = TRUE)
NAoutliers <- as.genind(tab(NAoutliers,NA.method="zero"))

#tt<-row.names(NAoutliers$tab)
#setdiff(inds2,tt)
#length(intersect(tt,inds2))

NAneutral <- adegenet::read.genepop("NA_Neutral_1712lociGP.GEN",
                                     ncode=3,quiet = TRUE)
NAneutral <- as.genind(tab(NAneutral,NA.method="zero"))

##Add the spatial coordinates

#######Outlier loci#####
NAOutlier_Cart <- NAoutliers # cartesian coordinates
NAOutlier_Cart@other$xy <- NACoords[,c("MDS1","MDS2")]

NAOutlier_jCart <- NAoutliers # jittered Cartesian coordinates (Delaunay triangulation)
NAOutlier_jCart@other$xy <- NACoords[,c("jCartx","jCarty")]

NAOutlier_Coords <- NAoutliers1 # Standard lat long
NAOutlier_Coords@other$xy <- NACoords[,c("Long","Lat")]

save.image("NA_SPCA_Data.RData")
gc()

###########Now run the SPCA#########
########outlier loci#################
####################################
#Neighbourhood distance (type=5)
NAOutlier_sPCA_Cart <- spca(NAOutlier_Cart,ask=F,type=5,d1=0,d2=3000,scannf = F)
barplot(NAOutlier_sPCA_Cart$eig)
screeplot(NAOutlier_sPCA_Cart)

x <- other(NAOutlier_Cart)$xy[,1]
y <- other(NAOutlier_Cart)$xy[,2]
png("NA_Outlier_SPCA_SpatialHeatmap.png",width=2400,height=2000,res=300)
temp <- interp(x, y, NAOutlier_sPCA_Cart$li[,1],duplicate = "mean")
image(temp, col=azur(100))
text(x = NACoords$Cartx,y=NACoords$Carty,labels = NACoords$Population,cex = 0.75)
dev.off()

#4. Delaunay (type=1)   
NAOutlier_sPCA_jCart <- spca(NAOutlier_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(NAOutlier_sPCA_jCart$eig,col=spectral(length(NAOutlier_sPCA_jCart$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")


x <- other(NAOutlier_Cart)$xy[,1]
y <- other(NAOutlier_Cart)$xy[,2]
temp2<-interp(x,y,NAOutlier_sPCA_jCart$li[,1],duplicate="mean")
png("NA_Outlier_SPCA_SpatialHeatmap_Delaunay.png",width=4800,height=4000,res=300)
image(temp2, col=azur(100))
text(x = NACoords$MDS1,y=NACoords$MDS2,labels = NACoords$Population,cex = 0.5)
dev.off()

myPal2 <- colorRampPalette(c("firebrick2","white","forestgreen"))

annot <- function(){
  title("sPCA - interpolated map of individual scores")
  text(x = NACoords$MDS1,y=NACoords$MDS2,labels = NACoords$Population,cex = 0.75)
}
filled.contour(temp2, color.pal=rainbow, nlev=30,
               key.title=title("lagged \nscore 1"), plot.title=annot())

save.image("NA_SPCA_Data.RData")
gc()


#################Neutral loci#########################

NAneutral_Cart <- NAneutral # cartesian coordinates
NAneutral_Cart@other$xy <- NACoords[,c("MDS1","MDS2")]
head(NAneutral_Cart@other$xy)
tail(NAneutral_Cart@other$xy)

NAneutral_jCart <- NAneutral # jittered Cartesian coordinates (Delaunay triangulation)
NAneutral_jCart@other$xy <- NACoords[,c("jCartx","jCarty")]
head(NAneutral_jCart@other$xy)
tail(NAneutral_jCart@other$xy)

NAneutral_Coords <- NAneutral # Standard lat long
NAneutral_Coords@other$xy <- NACoords[,c("Long","Lat")]

###Now run the SPCA
#Neighbourhood distance (type=5)
#NAneutral_sPCA_Cart <- spca(NAneutral_Cart,ask=TRUE,type=5,d1=0,d2=4000,scannf = FALSE)
#barplot(NAneutral_sPCA_Cart$eig)
#screeplot(NAneutral_sPCA_Cart)


#Delaunay
NANeutral_sPCA_jCart <- spca(NAneutral_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(NANeutral_sPCA_jCart$eig,col=spectral(length(NANeutral_sPCA_jCart$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

x <- other(NAneutral_Cart)$xy[,1]
y <- other(NAneutral_Cart)$xy[,2]
png("NA_neutral_SPCA_SpatialHeatmap.png",width=2400,height=2000,res=300)
temp3 <- interp(x, y, NANeutral_sPCA_jCart$li[,1],duplicate = "mean")
image(temp3, col=azur(100))
text(x = NACoords$MDS1,y=NACoords$MDS2,labels = NACoords$Population,cex = 0.75)
dev.off()


png("NA_neutral_Delaunay_InterpolatedMap.png",width=2400,height=2000,res=300)

annot <- function(){
  title("sPCA - interpolated map of individual scores")
  text(x = NACoords$MDS1,y=NACoords$MDS2,labels = NACoords$Population,cex = 0.75)
}
filled.contour(temp3, color.pal=rainbow, nlev=30,
               key.title=title("lagged \nscore 1"), plot.title=annot())

dev.off()


#Global and local tests for significance on outliers and neutral SNPs
myoutlierGtest<-global.rtest(NAOutlier_jCart@tab,NAOutlier_sPCA_jCart$lw,nperm=9999)
plot(myoutlierGtest)
myoutlierLtest<-local.rtest(NAOutlier_jCart@tab,NAOutlier_sPCA_jCart$lw,nperm=9999)
plot(myoutlierLtest)
myneutralGtest<-global.rtest(NAneutral_jCart@tab,NAOutlier_sPCA_jCart$lw,nperm=9999)
plot(myneutralGtest)
myneutralLtest<-local.rtest(NAneutral_jCart@tab,NAOutlier_sPCA_jCart$lw,nperm=9999)
plot(myneutralLtest)

save.image("NA_SPCA_Data.RData")
