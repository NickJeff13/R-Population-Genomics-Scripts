#Description ------------
## Run a Spatial principal components analysis for salmon using neutral and outlier loci
setwd("/home/ian/Desktop/Nick/AtlanticSalmon/SPCA/")
load("Europe_SPCA_Data.RData")

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
library(akima) #masks lazyeval interp function
library(splancs)
library(gridExtra)
library(G)
?chooseCN #this will show your options for the SPCA type

inds<-genepop_detective("Europe_19OutliersGP.txt",variable = "Inds")
duplicated(inds)
sum(duplicated(inds))
## Step 1 Load coordintes ------------

euroCoords <- read.csv("EuropeCartesianCoords.csv",header=T)
head(euroCoords)
inds2<-as.vector(euroCoords$SampleID)
duplicated(inds2)
sum(duplicated(inds2))
length(intersect(inds3,inds2))
setdiff(inds,inds2)

##Jitter the species coordieurotes so that a Delaueuroy triangulation can be applied. Here we use a jitter of 2 km 
euroCoords$jCartx <- jitter(euroCoords$MDS1,2) 
euroCoords$jCarty <- jitter(euroCoords$MDS2,2)

##### Step 2 Read in genepop files as GENIND objects and add coordieurotse------------

eurooutliers <- adegenet::read.genepop("Europe_Outliers.GEN",
                                       ncode=3,quiet = T)

eurooutliers <- as.genind(tab(eurooutliers,NA.method="zero"))

euroneutral <- adegenet::read.genepop("Europe_Neutral_withGilbey_GP.GEN",
                                      ncode=3,quiet = TRUE)
euroneutral <- as.genind(tab(euroneutral,NA.method="zero"))

##Add the spatial coordieurotes

#######Outlier loci#####
euroOutlier_Cart <- eurooutliers # cartesian coordieurotes
euroOutlier_Cart@other$xy <- euroCoords[,c("MDS1","MDS2")]

euroOutlier_jCart <- eurooutliers # jittered Cartesian coordieurotes (Delaueuroy triangulation)
euroOutlier_jCart@other$xy <- euroCoords[,c("jCartx","jCarty")]

euroOutlier_Coords <- eurooutliers # Standard lat long
euroOutlier_Coords@other$xy <- euroCoords[,c("Long","Lat")]

save.image("Europe_SPCA_Data.RData")
gc()

#########Now run the SPCA###########
########outlier loci#################
####################################
#Neighbourhood distance (type=5)
euroOutlier_sPCA_Cart <- spca(euroOutlier_Cart,ask=F,type=5,d1=0,d2=5000,scannf = F)
barplot(euroOutlier_sPCA_Cart$eig,col=spectral(length(euroOutlier_sPCA_Cart$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
#screeplot(euroOutlier_sPCA_Cart)
gc()
x <- other(euroOutlier_Cart)$xy[,1]
y <- other(euroOutlier_Cart)$xy[,2]
png("euro_Outlier_NeighbourhoodSPCA_SpatialHeatmap.png",width=2400,height=2000,res=300)
temp <- interp(x, y, euroOutlier_sPCA_Cart$li[,1],duplicate = "mean")
image(temp, col=azur(100))
text(x = euroCoords$MDS1,y=euroCoords$MDS2,labels = euroCoords$Population,cex = 0.75)
dev.off()

png("euro_Outlier_NeighbourhoodSPCA_SpatialHeatmap2.png",width=2400,height=2000,res=300)
annot <- function(){
  title("sPCA - interpolated map of individual scores")
  text(x = euroCoords$MDS1,y=euroCoords$MDS2,labels = euroCoords$Population,cex = 0.75)
}
filled.contour(temp, color.pal=myPal2, nlev=50,
               key.title=title("lagged \nscore 1"), plot.title=annot())
dev.off()
#gc()

#4. Delaunay (type=1)   
euroOutlier_sPCA_jCart <- spca(euroOutlier_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(euroOutlier_sPCA_jCart$eig,col=spectral(length(euroOutlier_sPCA_jCart$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")


temp2<-akima::interp(x,y,euroOutlier_sPCA_jCart$li[,1],duplicate="mean")
png("euro_Outlier_SPCA_SpatialHeatmap_Delaunay.png",width=2400,height=2000,res=300)
image(temp2, col=azur(100))
text(x = euroCoords$MDS1,y=euroCoords$MDS2,labels = euroCoords$Population,cex = 0.5)
dev.off()

png("euro_Outlier_SPCA_SpatialHeatmap_Delaunay.png",width=4800,height=4000,res=300)

myPal2 <- colorRampPalette(c("firebrick2","white","forestgreen"))

annot <- function(){
  title("sPCA - interpolated map of individual scores")
  text(x = euroCoords$MDS1,y=euroCoords$MDS2,labels = euroCoords$Population,cex = 0.75)
}
filled.contour(temp2, color.pal=rainbow, nlev=50,
               key.title=title("lagged \nscore 1"), plot.title=annot())
dev.off()

#Global and local tests for significance on outliers
myGtest<-global.rtest(euroOutlier_jCart@tab,euroOutlier_sPCA_jCart$lw,nperm=9999)
plot(myGtest)
myLtest<-local.rtest(euroOutlier_jCart@tab,euroOutlier_sPCA_jCart$lw,nperm=9999)
plot(myLtest)
#######################################################
#################Neutral loci#########################
#######################################################

euroneutral_Cart <- euroneutral # cartesian coordieurotes
euroneutral_Cart@other$xy <- euroCoords[,c("MDS1","MDS2")]
head(euroneutral_Cart@other$xy)
tail(euroneutral_Cart@other$xy)

euroneutral_jCart <- euroneutral # jittered Cartesian coordieurotes (Delaueuroy triangulation)
euroneutral_jCart@other$xy <- euroCoords[,c("jCartx","jCarty")]
head(euroneutral_jCart@other$xy)
tail(euroneutral_jCart@other$xy)

euroneutral_Coords <- euroneutral # Standard lat long
euroneutral_Coords@other$xy <- euroCoords[,c("Long","Lat")]


save.image("Europe_SPCA_Data.RData")

#####Now run the SPCA######
#Neighbourhood distance (type=5)
euroneutral_sPCA_Cart <- spca(euroneutral_Cart,ask=TRUE,type=5,d1=0,d2=5000,scannf = FALSE)


barplot(euroneutral_sPCA_Cart$eig)
#screeplot(euroneutral_sPCA_Cart)


#Delaunay
euroNeutral_sPCA_jCart <- spca(euroneutral_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(euroNeutral_sPCA_jCart$eig,col=spectral(length(euroNeutral_sPCA_jCart$eig)))
#legend("topright", fill=spectral(2),
 #      leg=c("Global structures", "Local structures"))
#abline(h=0,col="grey")

save.image("Europe_SPCA_Data.RData")


x <- other(euroneutral_Cart)$xy[,1]
y <- other(euroneutral_Cart)$xy[,2]
png("euro_neutral_SPCA_DelaunaySpatialHeatmap.png",width=2400,height=2000,res=300)
temp3 <- interp(x, y, euroNeutral_sPCA_jCart$li[,1],duplicate = "mean")
image(temp3, col=azur(100))
text(x = euroCoords$MDS1,y=euroCoords$MDS2,labels = euroCoords$Population,cex = 0.75)
dev.off()


png("euro_neutral_Delaunay_InterpolatedMap.png",width=2400,height=2000,res=300)


annot <- function(){
  title("sPCA - interpolated map of individual scores")
  text(x = euroCoords$MDS1,y=euroCoords$MDS2,labels = euroCoords$Population,cex = 0.75)
}
filled.contour(temp3, color.pal=rainbow, nlev=50,
               key.title=title("lagged \nscore 1"), plot.title=annot())

dev.off()

#Global and local tests for significance on neutral
myneutralGtest<-global.rtest(euroneutral_jCart@tab,euroNeutral_sPCA_jCart$lw,nperm=9999)
plot(myneutralGtest)
myneutralLtest<-local.rtest(euroneutral_jCart@tab,euroNeutral_sPCA_jCart$lw,nperm=9999)
plot(myneutralLtest)

#Save data again
save.image("Europe_SPCA_Data.RData")
