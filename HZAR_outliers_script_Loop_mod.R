#Nick's Atlantic salmon hzar script to analyze geographic clines in allele frequency

setwd("/home/ian/Desktop/Nick/GreenCrab/HZAR/")
rm(list=ls())

library(hzar)
library(doMC)




#This will make sure your models run in parallel
if(require(doMC)){
  registerDoMC()
} else {
  #use foreach in sequential mode
  registerDoSEQ();
}

## A typical chain length. This value is the default setting in the package.
chainLength=1e5

#Run each model off a different set of seeds
mainSeed=list(A=c(596,528,124,978,544,99),
              B=c(528,124,978,544,99,596),
              C=c(124,978,544,99,596,528),
              D=c(109,824,555,713,82,341),
              E=c(54,216,544,621,672,930),
              G=c(28,64,447,523,754,800))



#Read in the data. This should have a list of locations, the least cost distance to a 'reference' pop 
#(e.g. this reference will have a distance of 0)
#And Allele frequencies for each locus you want to test, plus sample sizes

crabout<-read.table("Crab_HZAROutliersIn.txt",header=TRUE)
#print(samout)

      rm(i)
      rm(outlier)
      if(length(apropos("^outlier$",ignore.case=FALSE)) == 0 ||
         + !is.list(outlier) ) outlier <- list()
      
      
      #outlist=c(1:4,9:45) #what loci you want to use.
      
      for(i in 1:45){
        
        ## Save all plots in a series of png files
       
        
        locusname=colnames(crabout)[i+2]
        locusnumber=colnames(crabout)[i+47]
      
        png(width=900, height=900, res=200, family="Arial", filename=paste0(locusname,"clinePlot.png"),pointsize=8)
       
        ## good to stay organized.
        outlier$locus <- list();
        outlier$locus$obs <- list();
        ## Space to hold the models to fit
        outlier$locus$models <- list();
        ## Space to hold the compiled fit requests
        outlier$locus$fitRs <- list();
        ## Space to hold the output data chains
        outlier$locus$runs <- list();
        ## Space to hold the analysed data
        outlier$locus$analysis <- list();
        
        outlier$locus$obs <- hzar.doMolecularData1DPops(crabout$Distance,
                                       crabout[,locusname],
                                       crabout[,locusnumber])
        #Plot the data (this will plot all your loci as separate pngs)
        hzar.plot.obsData(outlier$locus$obs)
       
        #####Make your models (4 models per locus)#####
        #outlier$locus$models$model1 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
        #                                                                scaling="fixed", tails="none")
        #outlier$locus$models$model2 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
        #                                                                scaling="fixed", tails="both")
        #outlier$locus$models$model3 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
        #                                                                scaling="free", tails="none")
        #outlier$locus$models$model4 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
        #                                                                scaling = "free",tails="both")
        #Graham Derryberry's Suggestion
        outlier$locus$models$model1 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
                                                            scaling="none", tails="none")
        outlier$locus$models$model2 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
                                                            scaling="none", tails="both")
        outlier$locus$models$model3 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
                                                            scaling="free", tails="none")
        outlier$locus$models$model4 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
                                                            scaling = "free",tails="both")
        outlier$locus$models$model5 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
                                                            scaling="fixed", tails="none")
        outlier$locus$models$model6 <- hzar.makeCline1DFreq(data=outlier$locus$obs, 
                                                            scaling="fixed", tails="both")
        
        print(outlier$locus$models)
        
        #data collected between 0 and 5100km
        outlier$locus$models <- sapply(outlier$locus$models, hzar.model.addBoxReq,
                                                   -100, 5200, simplify=FALSE)
        
        #Check parameters
        print(outlier$locus$models)
        
        
        #Compile models to prepare for fitting ---> creates hzar.fitRequest from each clineModel object
        outlier$locus$fitRs$init <- sapply(outlier$locus$models, 
                                                       hzar.first.fitRequest.old.ML,
                                                       obsData = outlier$locus$obs,
                                                       verbose=FALSE,
                                                       simplify=FALSE)
        #update settings for the fitter using chainLength and mainSeed created before
        outlier$locus$fitRs$init$model1$mcmcParam$chainLength <- chainLength
        outlier$locus$fitRs$init$model1$mcmcParam$burnin <- chainLength %/% 10
        outlier$locus$fitRs$init$model1$mcmcParam$seed[[1]] <- mainSeed$A
        
        outlier$locus$fitRs$init$model2$mcmcParam$chainLength <- chainLength
        outlier$locus$fitRs$init$model2$mcmcParam$burnin <- chainLength %/% 10
        outlier$locus$fitRs$init$model2$mcmcParam$seed[[1]] <- mainSeed$B
        
        outlier$locus$fitRs$init$model3$mcmcParam$chainLength <- chainLength
        outlier$locus$fitRs$init$model3$mcmcParam$burnin <- chainLength %/% 10
        outlier$locus$fitRs$init$model3$mcmcParam$seed[[1]] <- mainSeed$C
        
        outlier$locus$fitRs$init$model4$mcmcParam$chainLength <- chainLength
        outlier$locus$fitRs$init$model4$mcmcParam$burnin <- chainLength %/% 10
        outlier$locus$fitRs$init$model4$mcmcParam$seed[[1]] <- mainSeed$D
        
        outlier$locus$fitRs$init$model5$mcmcParam$chainLength <- chainLength
        outlier$locus$fitRs$init$model5$mcmcParam$burnin <- chainLength %/% 10
        outlier$locus$fitRs$init$model5$mcmcParam$seed[[1]] <- mainSeed$E
        
        outlier$locus$fitRs$init$model6$mcmcParam$chainLength <- chainLength
        outlier$locus$fitRs$init$model6$mcmcParam$burnin <- chainLength %/% 10
        outlier$locus$fitRs$init$model6$mcmcParam$seed[[1]] <- mainSeed$G
        
        #check fit request settings 
        print(outlier$locus$fitRs$init)
        
        #Run one model, just to check
        #outlier$locus$runs$init <- list()
        #outlier$locus$runs$init$model1 <-hzar.doFit(outlier$locus$fitRs$init$model1)
        #plot(hzar.mcmc.bindLL(outlier$locus$runs$init$model1))
        
        #Compile a new set of fit requests using the initial chains
        outlier$locus$fitRs$chains <- lapply(outlier$locus$init, 
                                           hzar.next.fitRequest)
        
        #replicate each fit request 3 times, keeping the original seeds 
        #while switching to a new seed channel
        #12 total fit requests - 4 models, 3 times each
        outlier$locus$fitRs$chains <- hzar.multiFitRequest(outlier$locus$fitRs$init,
                                                                       each=3,
                                                                       baseSeed=NULL)
        
        
        #######have 36 fit requests - models 4, each with 3 chain################################################## 
        #running the chain 3 times - 36 total runs? - THIS WILL TAKE A WHILE
        outlier$locus$runs$doSeq <- lapply(outlier$locus$fitRs$chains,
                                                       hzar.chain.doSeq, 
                                                       count = 3)
        
        names(outlier)[i]=locusname
        
        dev.off()
        
        gc() # just in case ram is limited
        
      }
      
      
      ####Save the data before trying analysis  
      save.image("HZAR_Models_Completed_CrabEnviroOutliers.RData")
      #load("HZAR_Models_Completed_CrabEnviroOutliers.RData")
      
      
      ################ANALYSIS Loop#############################
      rm(list=ls())
      load("HZAR_Models_Completed_CrabEnviroOutliers.RData")
      #rm(j)
      #j=1
      
      
      locusnames <- names(outlier)
      
      for(j in 1:45) {
        
        png(width=900, height=900, res=200, family="Arial", filename=paste0(locusnames[j],"clinePlot.png"),pointsize=8)
        
        #did model1 converge?
        outlier[[j]]$runs$doSeq[1:3]
        #YES
        
        #did model2 converge?
        outlier[[j]]$runs$doSeq[4:6]
        #YES
        
        #did model3 converge?
        outlier[[j]]$runs$doSeq[7:9]
        #YES
        
        #did model4 converge?
        outlier[[j]]$runs$doSeq[10:12]
        #YES
        
        #did model5 converge?
        outlier[[j]]$runs$doSeq[13:15]
        #YES
        
        #did model6 converge?
        outlier[[j]]$runs$doSeq[16:18]
        #YES
        
        #drop bad models
        
      #outlier[[j]]$runs$doSeq <- outlier[[j]]$runs$doSeq[c("model1","model1")] 
        
        #####ANALYSIS#####
        #start aggregation of data for analysis
        
        #create a data group for the null model
        outlier[[j]]$analysis$initDGs <- list(nullModel=hzar.dataGroup.null(outlier[[j]]$obs))
        
        #create a model data group for each model from the initial runs
       outlier[[j]]$analysis$initDGs$model1 <- hzar.dataGroup.add(outlier[[j]]$runs$doSeq$model1)
       if(is.null(outlier[[j]]$analysis$initDGs$model1)){rm(outlier[[j]]$analysis$initDGs$model1)}                                                      
      
       outlier[[j]]$analysis$initDGs$model2 <- hzar.dataGroup.add(outlier[[j]]$runs$doSeq$model2)
       if(is.null(outlier[[j]]$analysis$initDGs$model2)){rm(outlier[[j]]$analysis$initDGs$model2)}                                                    
      
        outlier[[j]]$analysis$initDGs$model3 <- hzar.dataGroup.add(outlier[[j]]$runs$doSeq$model3)
        if(is.null(outlier[[j]]$analysis$initDGs$model3)){rm(outlier[[j]]$analysis$initDGs$model3)}                                       
      
       outlier[[j]]$analysis$initDGs$model4 <- hzar.dataGroup.add(outlier[[j]]$runs$doSeq$model4)
         if(is.null(outlier[[j]]$analysis$initDGs$model4)){rm(outlier[[j]]$analysis$initDGs$model4)}                                             
      
        outlier[[j]]$analysis$initDGs$model5 <- hzar.dataGroup.add(outlier[[j]]$runs$doSeq$model5) 
        if(is.null(outlier[[j]]$analysis$initDGs$model5)){rm(outlier[[j]]$analysis$initDGs$model5)}
      
        outlier[[j]]$analysis$initDGs$model6 <- hzar.dataGroup.add(outlier[[j]]$runs$doSeq$model6)
        if(is.null(outlier[[j]]$analysis$initDGs$model6)){rm(outlier[[j]]$analysis$initDGs$model6)}                                                  
        
      
      #identify the models which worked
        workingmodels=names(outlier[[j]]$analysis$initDGs)
        workingmodels=rep(workingmodels[-grep("nullModel",workingmodels)],each=3) #subset out the nullModel from this list
        outlier[[j]]$runs$doSeq=outlier[[j]]$runs$doSeq[workingmodels]#subset out any models which did not converge as specified by the workingmodels list
        
        ##create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme
        
        outlier[[j]]$analysis$oDG<- hzar.make.obsDataGroup(outlier[[j]]$analysis$initDGs)
        
        outlier[[j]]$analysis$oDG <- hzar.copyModelLabels(outlier[[j]]$analysis$initDGs,
                                                                       outlier[[j]]$analysis$oDG)
        
        ##convert all 48 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object
        outlier[[j]]$analysis$oDG <- hzar.make.obsDataGroup(lapply(outlier[[j]]$runs$doSeq, 
                                                                                hzar.dataGroup.add),
                                                                         outlier[[j]]$analysis$oDG)
        
        #check to make sure there are only 7 hzar.dataGroup objects
        print(summary(outlier[[j]]$analysis$oDG$data.groups))
      
        #compare the 6 cline models to the null model graphically
        hzar.plot.cline(outlier[[j]]$analysis$oDG)
        
        
        ## Do model selection based on the AICc scores
        print(outlier[[j]]$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(outlier[[j]]$analysis$oDG))
        write.table(x=outlier[[j]]$analysis$AICcTable,file = paste0(locusnames[j],"_AIC_Table.txt"),quote=F)
        
        #print the model with the minimum AICc score
        print(outlier[[j]]$analysis$model.name <-
                 rownames(outlier[[j]]$analysis$AICcTable)[[which.min(outlier[[j]]$analysis$AICcTable$AICc )]])
        #[1] "model3"
        
        
        #Extract the hzar.dataGroup object for the selected model
        outlier[[j]]$analysis$model.selected <-
          outlier[[j]]$analysis$oDG$data.groups[[outlier[[j]]$analysis$model.name]]
        
        #look at the variation in parameters for the selected model
        #print(hzar.getLLCutParam(outlier[[j]]$analysis$model.selected,
         #                         names(outlier[[j]]$analysis$model.selected$data.param)))
        #center2LLLow center2LLHigh width2LLLow width2LLHigh pMin2LLLow pMin2LLHigh pMax2LLLow pMax2LLHigh
        #1     769.2505      1284.908    12.13201     1369.043  0.2200383   0.3115534  0.5054706   0.7085702
        
        
        ####Print Params####
        
        #print the Max. Likelihood cline width for the selected model
        outlier[[j]]$analysis$modeldetails <- print(hzar.get.ML.cline(outlier[[j]]$analysis$model.selected))
        print(outlier[[j]]$analysis$modeldetails$param.all$width)
        print(outlier[[j]]$analysis$modeldetails$logLike)
        write.table(x = outlier[[j]]$analysis$modeldetails$param.all, file = paste0(locusnames[j],"_cline.txt"),quote = FALSE)
       
        #plot the maximum likelihood cline for the selected model
        hzar.plot.cline(outlier[[j]]$analysis$model.selected)
        
        #plot the 95% credible cline region for the selected model
       if (outlier[[j]]$analysis$model.name!="nullModel"){
         hzar.plot.fzCline(outlier[[j]]$analysis$model.selected)
       }
         ###DONE###
        
        #names(outlier)[j]=locusname
        
        dev.off()
       
        gc() # just in case ram is limited
        
        }    #end of j loop
     
  save.image("Crab_Outliers_Modeled.RData")
      
      
      
dev.off()



####Pick loci that are actually clinal!####

for (k in 1:45){
#print all AICc scores
#print(outlier[[k]]$analysis$AICcTable)

#print the min AICc score
#print(min(outlier[[k]]$analysis$AICcTable$AICc))
print(outlier[[k]]$analysis$modeldetails$logLike) #minimize negative log-likelihood - < -20 is good for crab
  
#make subset of clines that have AICc values under a certain criteria, then plot
}




####Save data image####
save.image("~/Desktop/Nick/GreenCrab/HZAR/Crab_Outliers_LogLikelihood.RData")
load("Crab_Outliers_LogLikelihood.RData")













#
#
#
#
#
#
#
#
#
#########################OLD SCRIPT FOR ONE LOCUS AT A TIME#################################################3

###REPLACE THIS SECTION WITH MY OWN DATASET NAME AND ALLELE NAMES
## Blank out space in memory to hold molecular analysis
if(length(apropos("^outlier$",ignore.case=FALSE)) == 0 ||
   + !is.list(outlier) ) outlier <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
outlier$ESTNV_29129_700 <- list();
## Space to hold the observed data
outlier$ESTNV_29129_700$obs <- list();
## Space to hold the models to fit
outlier$ESTNV_29129_700$models <- list();


## Space to hold the compiled fit requests## Blank out space in memory to hold molecular analysis
if(length(apropos("^outlier$",ignore.case=FALSE)) == 0 ||
   + !is.list(outlier) ) outlier <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
outlier$ESTNV_29129_700 <- list();
## Space to hold the observed data
outlier$ESTNV_29129_700$obs <- list();
## Space to hold the models to fit
outlier$ESTNV_29129_700$models <- list();
## Space to hold the compiled fit requests
outlier$ESTNV_29129_700$fitRs <- list();
## Space to hold the output data chains
outlier$ESTNV_29129_700$runs <- list();
## Space to hold the analysed data
outlier$ESTNV_29129_700$analysis <- list();

# Locus ESTNV_29129_700 from Europe
outlier$ESTNV_29129_700$obs <-hzar.doMolecularData1DPops(samout$Distance,
                                                         samout$ESTNV_29129_700,
                                                         samout$NumberLocus4);
## Look at a graph of the observed data
hzar.plot.obsData(outlier$ESTNV_29129_700$obs)



#set the model I want to look at - check param options using ?hzar.makeCline1DFreq
#EXAMPLE - testing 4 models#
outlier$"ESTNV_29129_700"$models$model1 <- hzar.makeCline1DFreq(data=outlier$"ESTNV_29129_700"$obs, 
                                                                scaling="fixed", tails="none")
outlier$"ESTNV_29129_700"$models$model2 <- hzar.makeCline1DFreq(data=outlier$"ESTNV_29129_700"$obs, 
                                                                scaling="fixed", tails="both")
outlier$"ESTNV_29129_700"$models$model3 <- hzar.makeCline1DFreq(data=outlier$"ESTNV_29129_700"$obs, 
                                                                scaling="free", tails="none")
outlier$"ESTNV_29129_700"$models$model4 <- hzar.makeCline1DFreq(data=outlier$"ESTNV_29129_700"$obs, 
                                                                scaling = "free",tails="both")


#data collected between 0 and 5100km
outlier$"ESTNV_29129_700"$models <- sapply(outlier$"ESTNV_29129_700"$models, hzar.model.addBoxReq,
                                           -30, 5200, simplify=FALSE)

#Check parameters
print(outlier$"ESTNV_29129_700"$models)


#Compile models to prepare for fitting ---> creates hzar.fitRequest from each clineModel object
outlier$"ESTNV_29129_700"$fitRs$init <- sapply(outlier$"ESTNV_29129_700"$models, 
                                               hzar.first.fitRequest.old.ML,
                                               obsData = outlier$"ESTNV_29129_700"$obs,
                                               verbose=FALSE,
                                               simplify=FALSE)
#update settings for the fitter using chainLength and mainSeed created before
outlier$"ESTNV_29129_700"$fitRs$init$model1$mcmcParam$chainLength <- chainLength
outlier$"ESTNV_29129_700"$fitRs$init$model1$mcmcParam$burnin <- chainLength %/% 10
outlier$"ESTNV_29129_700"$fitRs$init$model1$mcmcParam$seed[[1]] <- mainSeed$A

outlier$"ESTNV_29129_700"$fitRs$init$model2$mcmcParam$chainLength <- chainLength
outlier$"ESTNV_29129_700"$fitRs$init$model2$mcmcParam$burnin <- chainLength %/% 10
outlier$"ESTNV_29129_700"$fitRs$init$model2$mcmcParam$seed[[1]] <- mainSeed$B

outlier$"ESTNV_29129_700"$fitRs$init$model3$mcmcParam$chainLength <- chainLength
outlier$"ESTNV_29129_700"$fitRs$init$model3$mcmcParam$burnin <- chainLength %/% 10
outlier$"ESTNV_29129_700"$fitRs$init$model3$mcmcParam$seed[[1]] <- mainSeed$C

outlier$"ESTNV_29129_700"$fitRs$init$model4$mcmcParam$chainLength <- chainLength
outlier$"ESTNV_29129_700"$fitRs$init$model4$mcmcParam$burnin <- chainLength %/% 10
outlier$"ESTNV_29129_700"$fitRs$init$model4$mcmcParam$seed[[1]] <- mainSeed$D

#check fit request settings 
print(outlier$"ESTNV_29129_700"$fitRs$init)


#replicate each fit request 3 times, keeping the original seeds 
#while switching to a new seed channel
#12 total fit requests - 4 models, 3 times each
outlier$"ESTNV_29129_700"$fitRs$chains <- hzar.multiFitRequest(outlier$"ESTNV_29129_700"$fitRs$init,
                                                               each=3,
                                                               baseSeed=NULL)


##have 36 fit requests - models 4, each with 3 chain####################################################### 
#running the chain 3 times - 36 total runs? - THIS WILL TAKE A WHILE
outlier$"ESTNV_29129_700"$runs$doSeq <- lapply(outlier$"ESTNV_29129_700"$fitRs$chains,
                                               hzar.chain.doSeq, 
                                               count = 3)




#did model1 converge?
summary(do.call(mcmc.list, lapply(outlier$"ESTNV_29129_700"$runs$doSeq[1:3],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES

#did model2 converge?
summary(do.call(mcmc.list, lapply(outlier$"ESTNV_29129_700"$runs$doSeq[4:6],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES

#did model3 converge?
summary(do.call(mcmc.list, lapply(outlier$"ESTNV_29129_700"$runs$doSeq[7:9],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES

#did model4 converge?
summary(do.call(mcmc.list, lapply(outlier$"ESTNV_29129_700"$runs$doSeq[10:12],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES


#####ANALYSIS#####
#start aggregation of data for analysis

#create a data group for the null model
outlier$"ESTNV_29129_700"$analysis$initDGs <- list(nullModel=hzar.dataGroup.null(outlier$"ESTNV_29129_700"$obs))

#create a model data group for each model from the initial runs
outlier$"ESTNV_29129_700"$analysis$initDGs$model1 <- hzar.dataGroup.add(outlier$"ESTNV_29129_700"$runs$doSeq$model1)
outlier$"ESTNV_29129_700"$analysis$initDGs$model2 <- hzar.dataGroup.add(outlier$"ESTNV_29129_700"$runs$doSeq$model2)
outlier$"ESTNV_29129_700"$analysis$initDGs$model3 <- hzar.dataGroup.add(outlier$"ESTNV_29129_700"$runs$doSeq$model3)
outlier$"ESTNV_29129_700"$analysis$initDGs$model4 <- hzar.dataGroup.add(outlier$"ESTNV_29129_700"$runs$doSeq$model4)


##create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme

outlier$"ESTNV_29129_700"$analysis$oDG <- hzar.make.obsDataGroup(outlier$"ESTNV_29129_700"$analysis$initDGs)

outlier$"ESTNV_29129_700"$analysis$oDG <- hzar.copyModelLabels(outlier$"ESTNV_29129_700"$analysis$initDGs,
                                                               outlier$"ESTNV_29129_700"$analysis$oDG)

##convert all 36 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object
outlier$"ESTNV_29129_700"$analysis$oDG <- hzar.make.obsDataGroup(lapply(outlier$"ESTNV_29129_700"$runs$doSeq, 
                                                                        hzar.dataGroup.add),
                                                                 outlier$"ESTNV_29129_700"$analysis$oDG)

#check to make sure there are only 5 hzar.dataGroup objects
print(summary(outlier$"ESTNV_29129_700"$analysis$oDG$data.groups))

#compare the 4 cline models to the null model graphically
hzar.plot.cline(outlier$"ESTNV_29129_700"$analysis$oDG)
#hzar.plot.fzCline(outlier$"ESTNV_29129_700"$analysis$model.selected)


print(outlier$"ESTNV_29129_700"$analysis$AICcTable <- 
        hzar.AICc.hzar.obsDataGroup(outlier$"ESTNV_29129_700"$analysis$oDG))



#print the model with the minimum AICc score
print(outlier$"ESTNV_29129_700"$analysis$model.name <-
        rownames(outlier$"ESTNV_29129_700"$analysis$AICcTable
        )[[which.min(outlier$"ESTNV_29129_700"$analysis$AICcTable$AICc )]])
#[1] "model3"



#Extract the hzar.dataGroup object for the selected model
outlier$"ESTNV_29129_700"$analysis$model.selected <-
  outlier$"ESTNV_29129_700"$analysis$oDG$data.groups[[outlier$"ESTNV_29129_700"$analysis$model.name]]

#look at the variation in parameters for the selected model
print(hzar.getLLCutParam(outlier$"ESTNV_29129_700"$analysis$model.selected,
                         names(outlier$"ESTNV_29129_700"$analysis$model.selected$data.param)))
#center2LLLow center2LLHigh width2LLLow width2LLHigh pMin2LLLow pMin2LLHigh pMax2LLLow pMax2LLHigh
#1     769.2505      1284.908    12.13201     1369.043  0.2200383   0.3115534  0.5054706   0.7085702


####Print Params####

#print the cline width for the selected model
outlier$"ESTNV_29129_700"$analysis$modeldetails <- print(hzar.get.ML.cline(outlier$"ESTNV_29129_700"$analysis$model.selected))
print(outlier$"ESTNV_29129_700"$analysis$modeldetails$param.all$width)

print(outlier$"ESTNV_29129_700"$analysis$modeldetails$logLike)


#plot the maximum likelihood cline for the selected model
hzar.plot.cline(outlier$"ESTNV_29129_700"$analysis$model.selected)

#plot the 95% credible cline region for the selected model
hzar.plot.fzCline(outlier$"ESTNV_29129_700"$analysis$model.selected)

write.table(x = outlier$ESTNV_29129_700$analysis$modeldetails$param.all,file = "ESTNV_29129_700cline.txt",quote = FALSE)
###DONE###
