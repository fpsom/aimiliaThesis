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

###################################################
#######            DISTRIBUTION           #########
###################################################

#FPKM CELLS#
histFPKMcells<- FPKMFINAL[c(2,3,4,5)] #choose the variabe to visualize the distribution
hist(log10(histFPKMcells))

#FPKM HUMAN# 
histFPKMhuman<- FPKMFINAL[c(6,7)] #choose the variabe to visualize the distribution
hist(log10(histFPKMhuman))

#FPKM MICE#
histFPKMmice<- FPKMFINAL[c(8,9,10,12,13)] #choose the variabe to visualize the distribution
hist(log10(histFPKMmice))

### DISTRIBUTION FC???????#######

################################################
########          CLUSTERRING           ########
################################################
#FPKM CELLS#
clustersFPKMcells <- hclust(dist(histFPKMcells))
plot(clustersFPKMcells) #Hierarhical Clustering

#FPKM HUMAN#
clustersFPKMhuman <- hclust(dist(histFPKMhuman))
plot(clustersFPKMhuman) #Hierarhical Clustering

#FPKM MICE#
clustersFPKMmice <- hclust(dist(histFPKMmice))
plot(clustersFPKMmice) #Hierarhical Clustering

################################################
######                PCA              #########
################################################
#FPKM CELLS#
pcaFPKMcells <- prcomp(histFPKMcells)
summary(pcaFPKMcells)# Prints variance summary for all principal components
biplot(pcaFPKMcells)

#FPKM human# 
pcaFPKMhuman <- prcomp(histFPKMhuman)
summary(pcaFPKMhuman)# Prints variance summary for all principal components
biplot(pcaFPKMhuman)

#FPKM MICE#
pcaFPKMmice <- prcomp(histFPKMmice)
summary(pcaFPKMmice)# Prints variance summary for all principal components
biplot(pcaFPKMmice)

#################################################
#########         REGRESSION             ########
#################################################

linearMod <- lm(FPKM.D0 ~ FPKM.D5, data=FPKMFINAL)  # build linear regression model  data
print(linearMod)
summary(linearMod)
AIC(linearMod)  
BIC(linearMod)


