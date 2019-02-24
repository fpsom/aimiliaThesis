######################################
#####        REQUIREMENTS      #######
######################################

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(DESeq2)
library(stats)
library(MASS)

######################################
#####      UPLOAD FILES        #######
######################################

cells <- read.csv(file="data/Cells.csv") 
human <- read.csv(file="data/Human.csv") 
mice <- read.csv(file="data/Mice.csv") 
orthologs_human<-read.csv(file="data/orthologs_human.csv") 
names(cells)
names(human)
names(mice)

####################################################
###      CREATE DATA WITH ORTHOLOGS GENES      #####
####################################################

orthologs<-orthologs_human[which(orthologs_human$Rank=='high'), ]
mice_orthologs<-merge(mice,orthologs, by="Mouse_Symbol")
names(mice_orthologs)

###################################################
 ######      FPKM and FC data of cells    #########
###################################################
FPKMcells<-cells[c(6,10,18,22,34)]
names(FPKMcells)
FCcells<-cells[c(34,36,39,42)]
names(FCcells)

#################################################
######      FPKM and FC data of human    ########
#################################################
FPKMhuman<-human[c(2,3,4)]
names(FPKMhuman)
FChuman<-human[c(4,5)]
names(FChuman)

##################################################
######     FPKM and FC of mice data       ########
##################################################
FPKMmice<-mice_orthologs[c(19,32,45,94,107,120,139)]
names(FPKMmice)
FCmice<-mice_orthologs[c(136,137,138,139)] 
names(FCmice)
  
###################################################
#######       FPKM FINAL DATASET       ############
###################################################
  
FPKM1<-merge(FPKMcells,FPKMhuman,by="gene_name")
FPKMFINAL<-merge(FPKM1,FPKMmice, by="gene_name")
names(FPKMFINAL)


###################################################
##########       FC FINAL DATASET        ##########
###################################################
FC1<-merge(FCcells,FChuman, by="gene_name")
FCFINAL<-merge(FC1,FCmice, by="gene_name")
names(FCFINAL)






