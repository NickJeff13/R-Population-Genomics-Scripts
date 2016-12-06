#Salmon Assignment Script
#Set working directory
setwd("/home/ian/Desktop/Nick/AtlanticSalmon/assigner")
#Load Libraries
rm(list=ls())
ls()
library(devtools)
library(reshape2)
library(ggplot2)
library(stringr)
library(stringi)
library(plyr)
library(dplyr) # load this package after plyr to work properly
library(tidyr)
library(readr)
library(ggvis)
#install.packages(pkgs = "~/Downloads/randomForestSRC", repos = NULL, type = "source")
library(randomForestSRC)
library(purrr)
#install_github("thierrygosselin/assigner")
library(stackr)
library(assigner)
#install_gsi_sim(fromSource=TRUE)
library(adegenet)

#Check gsi sim exists
gsi_sim_binary()

####################################################################################################
###1. Remove duplicate individuals and markers and convert tsv to 'tidy genome' format using stackr#
####################################################################################################
duplicate.markers.removed <- stackr::read_long_tidy_wide(data="Top96_Oct2016_AllPops_Assigner.tsv") %>% distinct(POP_ID, INDIVIDUALS, LOCUS, .keep_all = TRUE)
write_tsv(x = duplicate.markers.removed, path = "Top96_AllPops_DAPCclusters_duplicates_removed.tsv")

#Remove duplicate individuals, try genome=TRUE (takes a long time but more accurate)
duplicates <- find_duplicate_genome(data="Top96_Europe_and_NA_Assigner_duplicates_removed.tsv", 
                                    genome = FALSE, parallel.core = 10)

#This one with genome=TRUE takes forever - not recommended!
#duplicategenomes <- find_duplicate_genome(data="Top285_Europe_and_NA_Assigner_duplicate_markers_removed.tsv", 
 #                                         genome = TRUE, parallel.core = 12)

#Plot the duplicates
duplicates$violin.plot.distance
ggsave("Europe_NA_violin.plot.distance.duplicate.pdf",width=10,height=10,dpi=300,units="cm",useDingbats=F)
plot <- duplicates$jitter.plot.distance + 
  geom_hline(yintercept = 1000, color="red",linetype="longdash") +
  geom_hline(yintercept = 5000, color="red",linetype="longdash") +
  geom_hline(yintercept = 10000, color="red",linetype="longdash")
plot 
ggsave("EuropeNAjitter.plot.distance.duplicate.png",width=10,height=10,dpi=150,units="cm")

#to get the stats
duplicates$distance.stats

# to filter the dataset and get to know the duplicate samples
dup.filtered <- filter(.data = duplicates$distance, DISTANCE < 1000)
write_tsv(x = dup.filtered, path = "dup.filtered1000.tsv")

save.image(file = "DuplicateFilteringData.Rdata")

###############################################################################################
#################Now the actual assignment######################################################
###############################################################################################

#my pop labels
pop.labels = c("Gulf","Aquaculture","Inner","NS","Gaspesie","Upper_Shore",
               "Lower_Shore","SouthQC","Anticosti","Ungava","Lab1","Lab2","Lab3",
               "Lab4","Hunt","Melville","NL1","NL2","NL3",
               "NBT","Avalon")

#1. Assignment time. Sampling method=random, no imputation
test.assignment<-assignment_ngs(data="Top285_Europe_and_NA_Assigner.tsv",
                                assignment.analysis="gsi_sim",
                                common.markers = TRUE,
                                marker.number=285,
                                pop.levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28),
                                #pop.labels = pop.labels,
                                sampling.method = "random",
                                #thl=0.4,
                                iteration.method = 1,
                                subsample = 20,
                                iteration.subsample = 50,
                                #gsi_sim.filename = "Top288FST_ADE_Result.txt",
                                keep.gsi.files = FALSE, 
                                imputation.method=NULL, #changed this to NULL from FALSE as of the last update
                                #impute="allele",
                                #imputations.group="populations",
                                #num.tree=100,
                                #iteration.rf=10,
                                #split.number = 100,
                                folder="test3",
                                parallel.core = 12)



#Get a summary plot for each tested group of markers
assignment.results<-test.assignment$assignment
write.table(assignment.results,file="Top288_FST_IFA_Random_STRUCTUREregions.txt", quote=FALSE)
View(assignment.results)
plot<-test.assignment$plot.assignment + facet_grid(~CURRENT)
plot
ggsave("Top285_NA_and_Europe_Regions.png", height = 8, width = 10,dpi = 600)

#Assignment of separate populations using a whitelist, ranked, GSI, and no imputation
system.time(test.assignment1<-assignment_ngs(data="AllAtlantic_byPOP.tsv",
                                             assignment.analysis="gsi_sim",
                                             common.markers = TRUE,
                                             marker.number=c(50,100,250,500,1000,2000,"all"),
                                             pop.levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
                                                            31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55),
                                             sampling.method = "random",
                                             iteration.method = 10,
                                             gsi_sim.filename = "Test_allPOPs2.txt",
                                             keep.gsi.files = FALSE, 
                                             imputation.method = FALSE,
                                             folder="GSI.random.allPopulationsTest_10Iterations",
                                             parallel.core = 14))

plot1<-test.assignment1$plot.assignment + facet_wrap(nrow = 5,~CURRENT)
plot1
ggsave("GSI_random_allPops_faceted.png", height = 8, width = 12,dpi = 600)

assignment.results<-test.assignment1$assignment
write.table(assignment.results,file="GSI_random_allPOPsNoImput.txt", quote=FALSE)
View(assignment.results2)

#2. Assignment, method=ranked, imputations using random forest

system.time(test.assignment2<-assignment_ngs(data="Top288IFA_AssignerRegions.tsv",
                                             assignment.analysis="gsi_sim",
                                            common.markers = TRUE,
                                            marker.number="all",
                                            pop.levels = c(1,2,3,4,5,6,7,8,9,10,11),
                                            sampling.method = "ranked",
                                            thl=0.3,
                                            iteration.method = 3,
                                            gsi_sim.filename = "salmon.data.ranked.txt",
                                            keep.gsi.files = FALSE, 
                                            imputation.method = "rf",
                                            impute="allele",
                                            imputations.group="populations",
                                            num.tree=100,
                                            iteration.rf=10,
                                            verbose=FALSE,
                                            folder="RankedWithImputation1",
                                            parallel.core = 4))

#Plot
assignment.results2<-test.assignment2$assignment
View(assignment.results2)
plot2<-test.assignment2$plot.assignment + facet_grid(SUBSAMPLE~CURRENT)
plot2
ggsave("salmon.assignment_rankedWithImp.png", height = 8, width = 8,dpi = 600)
write.table(assignment.results,file="ResultsRankedWithImput.txt", quote=FALSE)


#3. Random with no imputations using ADEGENET

system.time(test.assignment3<-assignment_ngs(data="salmon.north.america.regions.tsv",
                                            assignment.analysis="adegenet",
                                            common.markers = TRUE,
                                            marker.number=c(50,100,200,500,1000,2000,"all"),
                                            pop.levels = c(1,2,3,4,5,6,7,8,9,10,11),
                                            sampling.method = "random",
                                            iteration.method = 3,
                                            gsi_sim.filename = "salmon.data.Adeg.txt",
                                            keep.gsi.files = FALSE, 
                                            imputations=FALSE,
                                            folder="RandomNoImputationAde1",
                                            parallel.core = 4)
)

#Get a summary plot for each tested group of markers
assignment.results3<-test.assignment3$assignment
View(assignment.results3)
plot3<-test.assignment3$plot.assignment + facet_grid(SUBSAMPLE~CURRENT)
plot3
ggsave("Salmon_Random_Adegenet_NoImput.png", height = 8, width = 10,dpi = 600)
write.table(assignment.results,file="Salmon_AssignmentResults_Adegenet.txt", quote=FALSE)

#4. Assignment with Adegenet and ranked method

system.time(test.assignment4<-assignment_ngs(data="salmon.north.america.regions.tsv",
                                             assignment.analysis="adegenet",
                                             common.markers = TRUE,
                                             marker.number=c(50,100,200,500,1000,2000,"all"),
                                             pop.levels = c(1,2,3,4,5,6,7,8,9,10,11),
                                             sampling.method = "ranked",
                                             thl=0.3,
                                             iteration.method = 3,
                                             gsi_sim.filename = "Adegenet_ranked_NoImput.txt",
                                             keep.gsi.files = FALSE, 
                                             imputations=FALSE,
                                             folder="RankedNoImputationAde1",
                                             parallel.core = 4)
)

#Get a summary plot for each tested group of markers
assignment.results4<-test.assignment4$assignment
View(assignment.results4)
plot4<-test.assignment4$plot.assignment + facet_grid(~CURRENT)
plot4
ggsave("Salmon_Ranked_Adegenet_NoImput.png", height = 8, width = 10,dpi = 600)
write.table(assignment.results4,file="Salmon_AssignmentResults_Adegenet_Ranked.txt", quote=FALSE)

#5. Assignment using ranked method and gsi_sim

system.time(test.assignment5<-assignment_ngs(data="salmon.north.america.regions.tsv",
                                             assignment.analysis="gsi_sim",
                                             common.markers = TRUE,
                                             marker.number=c(50,100,200,500,1000,2000,"all"),
                                             pop.levels = c(1,2,3,4,5,6,7,8,9,10,11),
                                             sampling.method = "ranked",
                                             thl=0.3,
                                             iteration.method = 5,
                                             gsi_sim.filename = "GSI_ranked_NoImput.txt",
                                             keep.gsi.files = FALSE, 
                                             imputations=FALSE,
                                             folder="RankedNoImputationGSI1",
                                             parallel.core = 4)
)

#Get a summary plot for each tested group of markers
assignment.results5<-test.assignment5$assignment
View(assignment.results5)
plot5<-test.assignment5$plot.assignment + facet_grid(~CURRENT)
plot5
ggsave("Salmon_Ranked_GSI_NoImput.png", height = 8, width = 10,dpi = 600)
write.table(assignment.results5,file="Salmon_AssignmentResults_GSI_Ranked.txt", quote=FALSE)
