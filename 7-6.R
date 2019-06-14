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
library(corrplot)
library(igraph)
library(DESeq2)
library(ggpubr)
######################################
#####      UPLOAD FILES        #######
######################################
cells <- read.csv(file="data/Cells.csv")
cells2<-read.csv(file="data/Cells2.csv")
cells3<-read.csv(file="data/Cells3.csv")
human <- read.csv(file="data/Human.csv") 
mice <- read.csv(file="data/Mice.csv") 
mice_orthologs<-read.csv(file="data/mice_orthologs.csv")
names(cells)
names(human)
names(mice)
names(mice_orthologs)
##########################################
#####     DESCRIPTION OF DATA      #######
##########################################
# 1 DATASET OF HUMAN:    csv:"human"-> 8013 GENES
#                        FPKM VALUES OF TWO CONDITIONS  ->"FPKM.IZ_TR_184_S2.Control" (CONTROL)
#                                                       ->"FPKM.IZ_TR_185_S3.SCA1"    (PAΤIENT)
#                                         LO2FC VALUES  ->"LOG2FChuman"
# 1 DATASET OF MICE:     csv:"mice" -> 22783 GENES
#                        FPKM VALUES OF SIX CODITIONS   -> "FPKM_W5_FVB" (CONTROL IN WEEK 5)
#                                                       -> "FPKM_W5_Q82" (PATIENT IN WEEK 5)
#                                                       -> "FPKM_W12_FVB" (CONTROL IN WEEK 12)
#                                                       -> "FPKM_W12_Q82" (PATIENT IN WEEK 12)
#                                                       -> "FPKM_W28_FVB" (CONTROL IN WEEK 28)
#                                                       -> "FPKM_W28_Q82" (PATIENT IN WEEK 28)
#                                         LO2FC VALUES  -> "LOG2FCmice_w5"
#                                                       -> "LOG2FCmice_w12"
#                                                       -> "LOG2FCmice_w28"
# 3 DATASETS OF CELLS:  csv:"cells"-> 3984 GENES
#                       FPKM VALUES OF 4 CONDITIONS     -> "FPKM_D0" (CONTROL)             
#                                                       -> "FPKM_D2" (PATIENT IN DAY 2)
#                                                       -> "FPKM_D5" (PATIENT IN DAY 5)
#                                                       -> "FPKM_D10" (PATIENT IN DAY 10)
#                                         LO2FC VALUES  -> "LOG2FCcellsD2.D0"
#                                                       -> "LOG2FCcellD5.D0"
#                                                       -> "LOG2FCcellsD10.D0"
#                       csv:"cells2"-> 688 GENES (FPKM VALUES AND LOG2FC VALUES AS CSV:"cells" RESPECTIVELY)
#                       csv:"cells3"-> 8591 GENES (FPKM VALUES AND LOG2FC VALUES AS CSV:"cells" RESPECTIVELY)
####################################################
###      CREATE DATA WITH ORTHOLOGS GENES      #####
####################################################
orthologs<-mice_orthologs[which(mice_orthologs$Rank=="high"), ] #in column Rank keep only "high"
mice_orthologs<-merge(mice,orthologs, by="Mouse_Symbol")        #merge by "Mouse Symbol"
names(mice_orthologs)
###################################################
########    FPKM AND LOG2FC OF HUMAN    ###########
###################################################
names(human)
FPKMhuman<-human[c("FPKM.IZ_TR_184_S2.Control","FPKM.IZ_TR_185_S3.SCA1","gene_name")] #crate dataframe FPKM human values and gene name
names(FPKMhuman)
FChuman<-human$"FPKM.IZ_TR_185_S3.SCA1"/human$"FPKM.IZ_TR_184_S2.Control" #create variable FChuman
LOG2FChuman<-log2(FChuman)                                                #create variable LOG2FChuman
LOG2FChuman
LOG2FChuman<-data.frame(human$"gene_name",LOG2FChuman)                    #create dataframe LOG2FChuman and gene name
names(LOG2FChuman)
names(LOG2FChuman)[names(LOG2FChuman) == 'human.gene_name'] <- 'gene_name' #rename the column 'human.gene_name'to "gene name"
names(LOG2FChuman)
###################################################                                 
#######     FPKM AND LOG2FC OF MICE      ##########
###################################################
names(mice_orthologs)
FPKMmice<-mice_orthologs[c("FPKM_W5_Q82","FPKM_W12_Q82","FPKM_W28_Q82","FPKM_W5_FVB","FPKM_W12_FVB","FPKM_W28_FVB","gene_name")] #crate dataframe FPKM mice values and gene name
names(FPKMmice)
FCmice_w5<-mice_orthologs$"FPKM_W5_Q82"/mice_orthologs$"FPKM_W5_FVB"   #create variable FCmice_w5
LOG2FCmice_w5<-log2(FCmice_w5)                                         #create variable LOG2FCmice_w5
FCmice_w12<-mice_orthologs$"FPKM_W12_Q82"/mice_orthologs$"FPKM_W12_FVB"#create variable FCmice_w12
LOG2FCmice_w12<-log2(FCmice_w12)                                       #create variable LOG2FCmice_w12
FCmice_w28<-mice_orthologs$"FPKM_W28_Q82"/mice_orthologs$"FPKM_W28_FVB"#create variable FCmice_w28
LOG2FCmice_w28<-log2(FCmice_w28)                                       #create variable LOG2FCmice_w28
LOG2FCmice<-data.frame(mice_orthologs$"gene_name",LOG2FCmice_w5,LOG2FCmice_w12,LOG2FCmice_w28) #create dataframe LOG2FCmice w5/w12/w28 and gene name
names(LOG2FCmice)
names(LOG2FCmice)[names(LOG2FCmice) == 'mice_orthologs.gene_name'] <- 'gene_name' #rename the column "mice_orthologs.gene_name" to "gene name"
names(LOG2FCmice)
###################################################
###   FPKM AND LOG2FC OF CELLS(3984_genes)     ####
###################################################
names(cells)
FPKMcells<-cells[c("FPKM_D0","FPKM_D2","FPKM_D5","FPKM_D10","gene_name")] #crate dataframe FPKM cells values and gene name
names(FPKMcells)
FCcellsD2.D0<-cells$"FPKM_D2"/cells$"FPKM_D0" #create variable FCcells_D0
LOG2FCcellsD2.D0<-log2(FCcellsD2.D0)          #create variable LO2FCcells_D0
FCcellsD5.D0<-cells$"FPKM_D5"/cells$"FPKM_D0" #create variable FCcells_D5
LOG2FCcellsD5.D0<-log2(FCcellsD5.D0)           #create variable LOG2FCcells_D5
FCcellsD10.D0<-cells$"FPKM_D10"/cells$"FPKM_D0" #create variable FCcells_D10
LOG2FCcellsD10.D0<-log2(FCcellsD10.D0)        #create variable LOG2FCcells_D10
LOG2FCcells<-data.frame(cells$"gene_name",LOG2FCcellsD2.D0,LOG2FCcellsD5.D0,LOG2FCcellsD10.D0) #create dataframe LOG2FCcells  and gene name
names(LOG2FCcells)
names(LOG2FCcells)[names(LOG2FCcells) == 'cells.gene_name'] <- 'gene_name' #rename the column 'cells.gene_name' to "gene name"
names(LOG2FCcells)
###################################################
######### LOG2FC_cells2(688_genes)    #############
###################################################
# there are not FPKM values 
names(cells2)
LOG2FCcells2D2.D0<-log2(cells2$"FC_D2.vs.D0" )
LOG2FCcells2D5.D0<-log2(cells2$"FC_D5.vs.D0" )
LOG2FCcells2D10.D0<-log2(cells2$"FC_D10.vs.D0" )
LOG2FCcells2<-data.frame(cells2$gene_name,LOG2FCcells2D2.D0,LOG2FCcells2D5.D0,LOG2FCcells2D10.D0)
names(LOG2FCcells2)[names(LOG2FCcells2) == 'cells2.gene_name'] <- 'gene_name'
names(LOG2FCcells2)
###################################################
######### LOG2FC_cells3(8591_genes)    ############
###################################################
names(cells3)
FCcells3D2.D0<-cells3$"FPKM_D2"/cells3$"FPKM_DO"
LOG2FCcells3D2.D0<-log2(FCcells3D2.D0)
FCcells3D5.D0<-cells3$"FPKM_D5"/cells3$"FPKM_DO"
LOG2FCcells3D5.D0<-log2(FCcells3D5.D0)
FCcells3D10.D0<-cells3$"FPKM_D10"/cells3$"FPKM_DO"
LOG2FCcells3D10.D0<-log2(FCcells3D10.D0)
LOG2FCcells3<-data.frame(cells3$gene_name,LOG2FCcells3D2.D0,LOG2FCcells3D5.D0,LOG2FCcells3D10.D0)
names(LOG2FCcells3)[names(LOG2FCcells3) == 'cells3.gene_name'] <- 'gene_name'
names(LOG2FCcells3)
###################################################
#     FPKM FINAL DATASET(cells:3984_genes)        #
###################################################
FPKM1<-merge(FPKMcells,FPKMhuman,by="gene_name")  # FPKM FINAL only for cells#
FPKMFINAL<-merge(FPKM1,FPKMmice, by="gene_name")
names(FPKMFINAL)
FPKMFINAL<- do.call(data.frame, lapply(FPKMFINAL, function(x) {replace(x, is.infinite(x) | is.na(x), 0)}))
###################################################
#     LOG2FC FINAL DATASET(cells:3984_genes)      #
###################################################
LOG2FC1<-merge(LOG2FCcells,LOG2FChuman, by="gene_name")
LOG2FCFINALcells<-merge(LOG2FC1,LOG2FCmice, by="gene_name")
names(LOG2FCFINALcells)
LOG2FCFINALcells<- do.call(data.frame, lapply(LOG2FCFINALcells, function(x) {replace(x, is.infinite(x) | is.na(x), 0)}))

###################################################
#     LOG2FC FINAL DATASET(cells:688_genes)       #
###################################################
LOG2FC2<-merge(LOG2FCcells2,LOG2FChuman, by="gene_name")
LOG2FCFINALcells2<-merge(LOG2FC2,LOG2FCmice, by="gene_name")
names(LOG2FCFINALcells2)
LOG2FCFINALcells2<- do.call(data.frame, lapply(LOG2FCFINALcells2, function(x) {replace(x, is.infinite(x) | is.na(x), 0)}))

###################################################
#     LOG2FC FINAL DATASET(cells:8591_genes)      #
###################################################
LOG2FC3<-merge(LOG2FCcells3,LOG2FChuman, by="gene_name")
LOG2FCFINALcells3<-merge(LOG2FC3,LOG2FCmice, by="gene_name")
names(LOG2FCFINALcells3)
LOG2FCFINALcells3
LOG2FCFINALcells3<- do.call(data.frame, lapply(LOG2FCFINALcells3, function(x) {replace(x, is.infinite(x) | is.na(x), 0)}))

###################################################
#######    DISTRIBUTION  OF DATA          #########
###################################################
#FPKM CELLS#
histFPKMcells<- FPKMFINAL[c("FPKM_D0","FPKM_D2","FPKM_D5","FPKM_D10")] 
hist(log10(histFPKMcells))

#FPKM HUMAN# 
histFPKMhuman<- FPKMFINAL[c("FPKM.IZ_TR_184_S2.Control","FPKM.IZ_TR_185_S3.SCA1")] 
hist(log10(histFPKMhuman))

#FPKM MICE#/
histFPKMmice<- FPKMFINAL[c("FPKM_W5_Q82","FPKM_W12_Q82","FPKM_W28_Q82","FPKM_W5_FVB","FPKM_W12_FVB","FPKM_W28_FVB")] 
hist(log10(histFPKMmice))

#LOG2FC CELLS#
histLOG2FCcells<-LOG2FCFINALcells[c("LOG2FCcellsD2.D0","LOG2FCcellsD5.D0","LOG2FCcellsD10.D0")]
hist(histLOG2FCcells)

#LOG2FC CELLS2#
histLOG2FCcells2<-LOG2FCFINALcells2[c("LOG2FCcells2D2.D0","LOG2FCcells2D5.D0","LOG2FCcells2D10.D0")]
hist(histLOG2FCcells2)

#LOG2FC CELLS3#
histLOG2FCcells3<-LOG2FCFINALcells3[c("LOG2FCcells3D2.D0","LOG2FCcells3D5.D0","LOG2FCcells3D10.D0")]
hist(histLOG2FCcells3)

#LOG2FC HUMAN#
histLOG2FChuman<-LOG2FCFINALcells[c("LOG2FChuman")]
hist(histLOG2FChuman)

#LOG2FC MICE#
histLOG2FCmice<-LOG2FCFINALcells[c("LOG2FCmice_w5","LOG2FCmice_w12","LOG2FCmice_w28")]
hist(histLOG2FCmice)

###################################################
###     EXRESSION LEVEL OF RPL-RPS GENES      #####
###################################################
names(LOG2FCFINALcells2)
nrow(LOG2FCFINALcells2)
View(LOG2FCFINALcells2)

names(LOG2FCFINALcells)
nrow(LOG2FCFINALcells)
View(LOG2FCFINALcells)

names(LOG2FCFINALcells3)
nrow(LOG2FCFINALcells3)
View(LOG2FCFINALcells3)

########################################################################################################################
########################################################################################################################

####################               COMRARISON HUMAN/MICE/THRESHOLD:(-1)-1                ###############################

########################################################################################################################
########################################################################################################################
names(LOG2FChuman)
names(LOG2FCmice)

LOG2FChuman<- do.call(data.frame, lapply(LOG2FChuman, function(x) {replace(x, is.infinite(x) | is.na(x), 0)}))
LOG2FCmice<- do.call(data.frame, lapply(LOG2FCmice, function(x) {replace(x, is.infinite(x) | is.na(x), 0)}))

threslog2fchuman<-LOG2FChuman[which(LOG2FChuman$LOG2FChuman<(-0.5) | LOG2FChuman$LOG2FChuman>0.5),]

LOG2FCmice_w28_df<-LOG2FCmice[c("gene_name","LOG2FCmice_w28")]
thresLOG2FCmice_w28_df<-LOG2FCmice_w28_df[which(LOG2FCmice_w28_df$LOG2FCmice_w28<(-0.5) | LOG2FCmice_w28_df$LOG2FCmice_w28>0.5),]

thresLOG2FChuman_mice<-merge(threslog2fchuman,thresLOG2FCmice_w28_df,by="gene_name")
names(thresLOG2FChuman_mice)
nrow(thresLOG2FChuman_mice)
View(thresLOG2FChuman_mice)
thresLOG2FChuman_mice

#correlation
thresLOG2FChuman_mice$gene_name<-NULL
corthresLOG2FChuman_mice<-cor(thresLOG2FChuman_mice)
corrplot(corthresLOG2FChuman_mice,method = "number")

#clustering
clustersthres1_d2.d0_W5 <- hclust(dist(thres1_d2.d0_W5))
plot(dist(thres1_d2.d0_W5) #Hierarhical Clustering

#pca
pcathres1_d2.d0_W5 <- prcomp(thres1_d2.d0_W5)
summary(pcathres1_d2.d0_W5)# Prints variance summary for all principal components
biplot(pcathres1_d2.d0_W5)


########################################################################################################################
########################################################################################################################

####################           COMRARISON HUMAN/MICE/CELLS_THRESHOLD:-1,1                ###############################

########################################################################################################################
########################################################################################################################
##################################
#    cells_dataset_688_genes     #
#comparison_D2.D0(cells)~W5(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells2)
thres1log2fc_d2.d0cells2<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCcells2D2.D0<(-1) | LOG2FCFINALcells2$LOG2FCcells2D2.D0 > 1), ]
thres1log2fc_d2.d0cells2<-thres1log2fc_d2.d0cells2[c("LOG2FCcells2D2.D0","gene_name")]
thres1log2fcmice_W5<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCmice_w5<(-1) | LOG2FCFINALcells2$LOG2FCmice_w5 > 1), ]
thres1log2fcmice_W5<-thres1log2fcmice_W5[c("LOG2FCmice_w5","gene_name")]
thres1_d2.d0_W5<-merge(thres1log2fc_d2.d0cells2,thres1log2fcmice_W5,by="gene_name")
thres1_d2.d0_W5

#correlation#DO NOT RUN BECAUSE THERE IS ONLY ONE GENE#
#thres1_d2.d0_W5$gene_name<-NULL
#corthres1_d2.d0_W5<-cor(thres1_d2.d0_W5)
#corrplot(corthres1_d2.d0_W5,method = "number")

#clustering#DO NOT RUN BECAUSE THERE IS ONLY ONE GENE#
#clustersthres1_d2.d0_W5 <- hclust(dist(thres1_d2.d0_W5))
#plot(dist(thres1_d2.d0_W5) #Hierarhical Clustering
     
#pca#DO NOT RUN BECAUSE THERE IS ONLY ONE GENE#
#pcathres1_d2.d0_W5 <- prcomp(thres1_d2.d0_W5)
#summary(pcathres1_d2.d0_W5)# Prints variance summary for all principal components
#biplot(pcathres1_d2.d0_W5)

##################################
#    cells_dataset_688_genes     #
#comparison_D5.D0(cells)~W12(mice)#
##################################
names(LOG2FCFINALcells2)
thres1log2fc_d5.d0cells2<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCcells2D5.D0<(-1) | LOG2FCFINALcells2$LOG2FCcells2D5.D0 > 1), ]
thres1log2fc_d5.d0cells2<-thres1log2fc_d5.d0cells2[c("LOG2FCcells2D5.D0","gene_name")]
thres1log2fcmice_W12<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCmice_w12<(-1) | LOG2FCFINALcells2$LOG2FCmice_w12 > 1), ]
thres1log2fcmice_W12<-thres1log2fcmice_W12[c("LOG2FCmice_w12","gene_name")]
thres1_d5.d0_W12<-merge(thres1log2fc_d5.d0cells2,thres1log2fcmice_W12,by="gene_name")
thres1_d5.d0_W12

#correlation#
thres1_d5.d0_W12$gene_name<-NULL
corthres1_d5.d0_W12<-cor(thres1_d5.d0_W12)
corrplot(corthres1_d5.d0_W12,method = "number")

#clustering#
#histthres1_d5.d0_W12<-thres1_d5.d0_W12[c("LOG2FCcells2D5.D0","LOG2FCmice_w12")]
#clustersthres1_d5.d0_W12 <- hclust(dist(thres1_d5.d0_W12))
#plot(clustersthres1_d5.d0_W12) #Hierarhical Clustering

#pca#
pcathres1_d5.d0_W12 <- prcomp(histthres1_d5.d0_W12)
summary(pcathres1_d5.d0_W12)# Prints variance summary for all principal components
biplot(pcathres1_d5.d0_W12)


##############################################################
##############    cells_dataset_688_genes     ################
##        comparison_D10.D0(cells)~W28(mice)~human         ###  
##############################################################
names(LOG2FCFINALcells2)
thres1log2fc_d10.d0cells2<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCcells2D10.D0<(-1) | LOG2FCFINALcells2$LOG2FCcells2D10.D0 > 1), ]
thres1log2fc_d10.d0cells2<-thres1log2fc_d10.d0cells2[c("LOG2FCcells2D10.D0","gene_name")]
thres1log2fcmice_W28<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCmice_w28<(-1) | LOG2FCFINALcells2$LOG2FCmice_w28 > 1), ]
thres1log2fcmice_W28<-thres1log2fcmice_W28[c("LOG2FCmice_w28","gene_name")]
thres1_d10.d0_W28<-merge(thres1log2fc_d10.d0cells2,thres1log2fcmice_W28,by="gene_name")
thres1_d10.d0_W28
thres1log2fc_human<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FChuman<(-1) | LOG2FCFINALcells2$LOG2FChuman > 1), ]
thres1log2fc_human<-thres1log2fc_human[c("LOG2FChuman","gene_name")]
thres1_d10.d0_W28_human<-merge(thres1_d10.d0_W28,thres1log2fc_human,by="gene_name")
thres1_d10.d0_W28_human

#correlation#
thres1_d10.d0_W28_human$gene_name<-NULL
corthres1_d10.d0_W28_human<-cor(thres1_d10.d0_W28_human)
corrplot(corthres1_d10.d0_W28_human,method = "number")

#clustering
#histthres1_d10.d0_W28_human<-thres1_d10.d0_W28_human[c("LOG2FCcells2D10.D0","LOG2FCmice_w28","LOG2FChuman")]
#hist(log10(histthres1_d10.d0_W28_human))
#clustersthres1_d10.d0_W28_human <- hclust(dist(histthres1_d10.d0_W28_human))
#plot(clustersthres1_d10.d0_W28_human) #Hierarhical Clustering

#pca#
pcathres1_d10.d0_W28_human <- prcomp(histthres1_d10.d0_W28_human)
summary(pcathres1_d10.d0_W28_human)# Prints variance summary for all principal components
biplot(pcathres1_d10.d0_W28_human)

##################################
#    cells_dataset_3984_genes     #
#comparison_D2.D0(cells)~W5(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells)
thres2log2fc_d2.d0cells<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCcellsD2.D0<(-1) | LOG2FCFINALcells$LOG2FCcellsD2.D0 > 1), ]
thres2log2fc_d2.d0cells<-thres2log2fc_d2.d0cells[c("LOG2FCcellsD2.D0","gene_name")]
thres2log2fcmice_W5<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCmice_w5<(-1) | LOG2FCFINALcells$LOG2FCmice_w5 > 1), ]
thres2log2fcmice_W5<-thres2log2fcmice_W5[c("LOG2FCmice_w5","gene_name")]
thres2_d2.d0_W5<-merge(thres2log2fc_d2.d0cells,thres2log2fcmice_W5,by="gene_name")
thres2_d2.d0_W5

#correlation#
thres2_d2.d0_W5$gene_name<-NULL
corthres2_d2.d0_W5<-cor(thres2_d2.d0_W5)
corrplot(corthres2_d2.d0_W5,method = "number")

#clustering#
histthres2_d2.d0_W5<-thres2_d2.d0_W5[c("LOG2FCcellsD2.D0","LOG2FCmice_w5")]
clustersthres2_d2.d0_W5 <- hclust(dist(histthres2_d2.d0_W5))
plot(clustersthres2_d2.d0_W5) #Hierarhical Clustering

#pca
pcathres2_d2.d0_W5 <- prcomp(thres2_d2.d0_W5)
summary(pcathres2_d2.d0_W5)# Prints variance summary for all principal components
biplot(pcathres2_d2.d0_W5)

##################################
#    cells_dataset_3984_genes     #
#comparison_D5.D0(cells)~W12(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells)
thres2log2fc_d5.d0cells<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCcellsD5.D0<(-1) | LOG2FCFINALcells$LOG2FCcellsD5.D0 > 1), ]
thres2log2fc_d5.d0cells<-thres2log2fc_d5.d0cells[c("LOG2FCcellsD5.D0","gene_name")]
thres2log2fcmice_W12<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCmice_w12<(-1) | LOG2FCFINALcells$LOG2FCmice_w12 > 1), ]
thres2log2fcmice_W12<-thres2log2fcmice_W12[c("LOG2FCmice_w12","gene_name")]
thres2_d5.d0_W12<-merge(thres2log2fc_d5.d0cells,thres2log2fcmice_W12,by="gene_name")
thres2_d5.d0_W12

#correlation#
thres2_d5.d0_W12$gene_name<-NULL
corthres2_d5.d0_W12<-cor(thres2_d5.d0_W12)
corrplot(corthres2_d5.d0_W12,method = "number")

#clustering#
histthres2_d5.d0_W12<-thres2_d5.d0_W12[c("LOG2FCcellsD5.D0","LOG2FCmice_w12")]
clustersthres2_d5.d0_W12 <- hclust(dist(histthres2_d5.d0_W12))
plot(clustersthres2_d5.d0_W12) #Hierarhical Clustering

#pca
pcathres2_d5.d0_W12 <- prcomp(thres2_d5.d0_W12)
summary(pcathres2_d5.d0_W12)# Prints variance summary for all principal components
biplot(pcathres2_d5.d0_W12)


##############################################################
##############    cells_dataset_3984_genes     ################
##        comparison_D10.D0(cells)~W28(mice)~human         ###  
##############################################################
names(LOG2FCFINALcells)
thres2log2fc_d10.d0cells<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCcellsD10.D0<(-1) | LOG2FCFINALcells$LOG2FCcellsD10.D0 > 1), ]
thres2log2fc_d10.d0cells<-thres2log2fc_d10.d0cells[c("LOG2FCcellsD10.D0","gene_name")]
thres2log2fcmice_W28<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCmice_w28<(-1) | LOG2FCFINALcells$LOG2FCmice_w28 > 1), ]
thres2log2fcmice_W28<-thres2log2fcmice_W28[c("LOG2FCmice_w28","gene_name")]
thres2_d10.d0_W28<-merge(thres2log2fc_d10.d0cells,thres2log2fcmice_W28,by="gene_name")
thres2_d10.d0_W28
thres2log2fc_human<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FChuman<(-1) | LOG2FCFINALcells$LOG2FChuman > 1), ]
thres2log2fc_human<-thres2log2fc_human[c("LOG2FChuman","gene_name")]
thres2_d10.d0_W28_human<-merge(thres2_d10.d0_W28,thres2log2fc_human,by="gene_name")
thres2_d10.d0_W28_human

#correlation#
thres2_d10.d0_W28_human$gene_name<-NULL
corthres2_d10.d0_W28_human<-cor(thres2_d10.d0_W28_human)
corrplot(corthres2_d10.d0_W28_human,method = "number")

#clustering#
histthres2_d10.d0_W28_human<-thres2_d10.d0_W28_human[c("LOG2FCcellsD10.D0","LOG2FCmice_w28","LOG2FChuman")]
clustersthres2_d5.d0_W12 <- hclust(dist(histthres2_d10.d0_W28_human))
plot(clustersthres2_d5.d0_W12) #Hierarhical Clustering


#pca#DO NOT RUN BECAUSE THERE IS ONLY ONE GENE#
pcathres2_d10.d0_W28_human <- prcomp(histthres2_d10.d0_W28_human)
summary(pcathres2_d10.d0_W28_human)# Prints variance summary for all principal components
biplot(pcathres2_d10.d0_W28_human)

##################################
#    cells_dataset_8591_genes     #
#comparison_D2.D0(cells)~W5(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells3)
thres3log2fc_d2.d0cells<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCcells3D2.D0<(-1) | LOG2FCFINALcells3$LOG2FCcells3D2.D0 > 1), ]
thres3log2fc_d2.d0cells<-thres3log2fc_d2.d0cells[c("LOG2FCcells3D2.D0","gene_name")]
thres3log2fcmice_W5<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCmice_w5<(-1) | LOG2FCFINALcells3$LOG2FCmice_w5 > 1), ]
thres3log2fcmice_W5<-thres3log2fcmice_W5[c("LOG2FCmice_w5","gene_name")]
thres3_d2.d0_W5<-merge(thres3log2fc_d2.d0cells,thres3log2fcmice_W5,by="gene_name")
thres3_d2.d0_W5

#correlation#
thres3_d2.d0_W5$gene_name<-NULL
corthres3_d2.d0_W5<-cor(thres3_d2.d0_W5)
corrplot(corthres3_d2.d0_W5,method = "number")

#clustering#
histthres3_d2.d0_W5<-thres3_d2.d0_W5[c("LOG2FCcells3D2.D0","LOG2FCmice_w5")]
clustersthres3_d2.d0_W5 <- hclust(dist(histthres3_d2.d0_W5))
plot(clustersthres3_d2.d0_W5) #Hierarhical Clustering

#pca
pcathres3_d2.d0_W5 <- prcomp(thres3_d2.d0_W5)
summary(pcathres3_d2.d0_W5)# Prints variance summary for all principal components
biplot(pcathres3_d2.d0_W5)

##################################
#    cells_dataset_8591_genes     #
#comparison_D5.D0(cells)~W12(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells3)
thres3log2fc_d5.d0cells<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCcells3D5.D0<(-1) | LOG2FCFINALcells3$LOG2FCcells3D5.D0 > 1), ]
thres3log2fc_d5.d0cells<-thres3log2fc_d5.d0cells[c("LOG2FCcells3D5.D0","gene_name")]
thres3log2fcmice_W12<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCmice_w12<(-1) | LOG2FCFINALcells3$LOG2FCmice_w12 > 1), ]
thres3log2fcmice_W12<-thres3log2fcmice_W12[c("LOG2FCmice_w12","gene_name")]
thres3_d5.d0_W12<-merge(thres3log2fc_d5.d0cells,thres3log2fcmice_W12,by="gene_name")
thres3_d5.d0_W12

#correlation#
thres3_d5.d0_W12$gene_name<-NULL
corthres3_d5.d0_W12<-cor(thres3_d5.d0_W12)
corrplot(corthres3_d5.d0_W12,method = "number")

#clustering#
histthres3_d5.d0_W12<-thres3_d5.d0_W12[c("LOG2FCcells3D5.D0","LOG2FCmice_w12")]
clustersthres3_d5.d0_W12 <- hclust(dist(histthres3_d5.d0_W12))
plot(clustersthres3_d5.d0_W12) #Hierarhical Clustering

#pca
pcathres3_d5.d0_W12 <- prcomp(thres3_d5.d0_W12)
summary(pcathres3_d5.d0_W12)# Prints variance summary for all principal components
biplot(pcathres3_d5.d0_W12)


##############################################################
##############    cells_dataset_8591_genes     ################
##        comparison_D10.D0(cells)~W25(mice)~human         ###  
##############################################################
names(LOG2FCFINALcells3)
thres3log2fc_d10.d0cells<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCcells3D10.D0<(-1) | LOG2FCFINALcells3$LOG2FCcells3D10.D0 > 1), ]
thres3log2fc_d10.d0cells<-thres3log2fc_d10.d0cells[c("LOG2FCcells3D10.D0","gene_name")]
thres3log2fcmice_W28<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCmice_w28<(-1) | LOG2FCFINALcells3$LOG2FCmice_w28 > 1), ]
thres3log2fcmice_W28<-thres3log2fcmice_W28[c("LOG2FCmice_w28","gene_name")]
thres3_d10.d0_W28<-merge(thres3log2fc_d10.d0cells,thres3log2fcmice_W28,by="gene_name")
thres3_d10.d0_W28
thres3log2fc_human<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FChuman<(-1) | LOG2FCFINALcells3$LOG2FChuman > 1), ]
thres3log2fc_human<-thres3log2fc_human[c("LOG2FChuman","gene_name")]
thres3_d10.d0_W28_human<-merge(thres3_d10.d0_W28,thres3log2fc_human,by="gene_name")
thres3_d10.d0_W28_human

#correlation#
thres3_d10.d0_W28_human$gene_name<-NULL
corthres3_d10.d0_W28_human<-cor(thres3_d10.d0_W28_human)
corrplot(corthres3_d10.d0_W28_human,method = "number")

#clustering#
histthres3_d10.d0_W28_human<-thres3_d10.d0_W28_human[c("LOG2FCcells3D10.D0","LOG2FCmice_w28","LOG2FChuman")]
clustersthres3_d5.d0_W12 <- hclust(dist(histthres3_d10.d0_W28_human))
plot(clustersthres3_d5.d0_W12) #Hierarhical Clustering


#pca#DO NOT RUN BECAUSE THERE IS ONLY ONE GENE#
pcathres3_d10.d0_W28_human <- prcomp(histthres3_d10.d0_W28_human)
summary(pcathres3_d10.d0_W28_human)# Prints variance summary for all principal components
biplot(pcathres3_d10.d0_W28_human)


########################################################################################################################
########################################################################################################################

################          COMRARISON HUMAN/MICE/CELLS_THRESHOLD:-(0,5),(0,5 )              #############################

########################################################################################################################
########################################################################################################################
##################################
#    cells_dataset_688_genes     #
#comparison_D2.D0(cells)~W5(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells2)
thres1log2fc_d2.d0cells2<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCcells2D2.D0<(-0.5) | LOG2FCFINALcells2$LOG2FCcells2D2.D0 > 0.5), ]
thres1log2fc_d2.d0cells2<-thres1log2fc_d2.d0cells2[c("LOG2FCcells2D2.D0","gene_name")]
thres1log2fcmice_W5<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCmice_w5<(-0.5) | LOG2FCFINALcells2$LOG2FCmice_w5 > 0.5), ]
thres1log2fcmice_W5<-thres1log2fcmice_W5[c("LOG2FCmice_w5","gene_name")]
thres1_d2.d0_W5<-merge(thres1log2fc_d2.d0cells2,thres1log2fcmice_W5,by="gene_name")
thres1_d2.d0_W5

#correlation
thres1_d2.d0_W5$gene_name<-NULL
corthres1_d2.d0_W5<-cor(thres1_d2.d0_W5)
corrplot(corthres1_d2.d0_W5,method = "number")

#clustering#
histthres1_d2.d0_W5<-thres1_d2.d0_W5[c("LOG2FCcells2D2.D0","LOG2FCmice_w5")]
clustersthres1_d2.d0_W5<- hclust(dist(thres1_d2.d0_W5))
plot(clustersthres1_d2.d0_W5) #Hierarhical Clustering

#pca#
pcathres1_d2.d0_W5 <- prcomp(histthres1_d2.d0_W5)
summary(pcathres1_d2.d0_W5)# Prints variance summary for all principal components
biplot(pcathres1_d2.d0_W5)

##################################
#    cells_dataset_688_genes     #
#comparison_D5.D0(cells)~W12(mice)#
##################################
names(LOG2FCFINALcells2)
thres1log2fc_d5.d0cells2<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCcells2D5.D0<(-0.5) | LOG2FCFINALcells2$LOG2FCcells2D5.D0 > 0.5), ]
thres1log2fc_d5.d0cells2<-thres1log2fc_d5.d0cells2[c("LOG2FCcells2D5.D0","gene_name")]
thres1log2fcmice_W12<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCmice_w12<(-0.5) | LOG2FCFINALcells2$LOG2FCmice_w12 > 0.5), ]
thres1log2fcmice_W12<-thres1log2fcmice_W12[c("LOG2FCmice_w12","gene_name")]
thres1_d5.d0_W12<-merge(thres1log2fc_d5.d0cells2,thres1log2fcmice_W12,by="gene_name")
thres1_d5.d0_W12
thres1_d5.d0_W12[-c(24:38), ]

#correlation#
thres1_d5.d0_W12$gene_name<-NULL
corthres1_d5.d0_W12<-cor(thres1_d5.d0_W12)
corrplot(corthres1_d5.d0_W12,method = "number")

#clustering#
histthres1_d5.d0_W12<-thres1_d5.d0_W12[c("LOG2FCcells2D5.D0","LOG2FCmice_w12")]
clustersthres1_d5.d0_W12 <- hclust(dist(thres1_d5.d0_W12))
plot(clustersthres1_d5.d0_W12) #Hierarhical Clustering

#pca#
pcathres1_d5.d0_W12 <- prcomp(histthres1_d5.d0_W12)
summary(pcathres1_d5.d0_W12)# Prints variance summary for all principal components
biplot(pcathres1_d5.d0_W12)


##############################################################
##############    cells_dataset_688_genes     ################
##        comparison_D10.D0(cells)~W28(mice)~human         ###  
##############################################################
names(LOG2FCFINALcells2)
thres1log2fc_d10.d0cells2<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCcells2D10.D0<(-0.5) | LOG2FCFINALcells2$LOG2FCcells2D10.D0 >0.5),]
thres1log2fc_d10.d0cells2<-thres1log2fc_d10.d0cells2[c("LOG2FCcells2D10.D0","gene_name")]
thres1log2fcmice_W28<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FCmice_w28<(-0.5)| LOG2FCFINALcells2$LOG2FCmice_w28 > 0.5), ]
thres1log2fcmice_W28<-thres1log2fcmice_W28[c("LOG2FCmice_w28","gene_name")]
thres1_d10.d0_W28<-merge(thres1log2fc_d10.d0cells2,thres1log2fcmice_W28,by="gene_name")
thres1_d10.d0_W28
thres1log2fc_human<-LOG2FCFINALcells2[ which(LOG2FCFINALcells2$LOG2FChuman<(-0.5) | LOG2FCFINALcells2$LOG2FChuman > 0.5), ]
thres1log2fc_human<-thres1log2fc_human[c("LOG2FChuman","gene_name")]
thres1_d10.d0_W28_human<-merge(thres1_d10.d0_W28,thres1log2fc_human,by="gene_name")
thres1_d10.d0_W28_human
names(thres1_d10.d0_W28_human)

#clustering
histthres1_d10.d0_W28_human<-thres1_d10.d0_W28_human[c("LOG2FCcells2D10.D0","LOG2FCmice_w28","LOG2FChuman")]
clustersthres1_d10.d0_W28_human <- hclust(dist(thres1_d10.d0_W28_human[-1]))
plot(clustersthres1_d10.d0_W28_human, labels=thres1_d10.d0_W28_human$gene_name) #Hierarhical Clustering
rect.hclust(clustersthres1_d10.d0_W28_human, 7)


#correlation#
thres1_d10.d0_W28_human$gene_name<-NULL
corthres1_d10.d0_W28_human<-cor(thres1_d10.d0_W28_human)
corrplot(corthres1_d10.d0_W28_human,method = "number")


#pca#
pcathres1_d10.d0_W28_human <- prcomp(thres1_d10.d0_W28_human)
summary(pcathres1_d10.d0_W28_human)# Prints variance summary for all principal components
biplot(pcathres1_d10.d0_W28_human)


##################################
#    cells_dataset_3984_genes     #
#comparison_D2.D0(cells)~W5(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells)
thres2log2fc_d2.d0cells<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCcellsD2.D0<(-0.5) | LOG2FCFINALcells$LOG2FCcellsD2.D0 > 0.5), ]
thres2log2fc_d2.d0cells<-thres2log2fc_d2.d0cells[c("LOG2FCcellsD2.D0","gene_name")]
thres2log2fcmice_W5<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCmice_w5<(-0.5) | LOG2FCFINALcells$LOG2FCmice_w5 > 0.5), ]
thres2log2fcmice_W5<-thres2log2fcmice_W5[c("LOG2FCmice_w5","gene_name")]
thres2_d2.d0_W5<-merge(thres2log2fc_d2.d0cells,thres2log2fcmice_W5,by="gene_name")
thres2_d2.d0_W5
thres2_d2.d0_W5[-c(5), ]

#correlation#
thres2_d2.d0_W5$gene_name<-NULL
corthres2_d2.d0_W5<-cor(thres2_d2.d0_W5)
corrplot(corthres2_d2.d0_W5,method = "number")

#clustering#
#histthres2_d2.d0_W5<-thres2_d2.d0_W5[c("LOG2FCcellsD2.D0","LOG2FCmice_w5")]
#clustersthres2_d2.d0_W5 <- hclust(dist(histthres2_d2.d0_W5))
#plot(clustersthres2_d2.d0_W5) #Hierarhical Clustering

#pca
#pcathres2_d2.d0_W5 <- prcomp(thres2_d2.d0_W5)
#summary(pcathres2_d2.d0_W5)# Prints variance summary for all principal components
#biplot(pcathres2_d2.d0_W5)

##################################
#    cells_dataset_3984_genes     #
#comparison_D5.D0(cells)~W12(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells)
thres2log2fc_d5.d0cells<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCcellsD5.D0<(-0.5) | LOG2FCFINALcells$LOG2FCcellsD5.D0 > 0.5), ]
thres2log2fc_d5.d0cells<-thres2log2fc_d5.d0cells[c("LOG2FCcellsD5.D0","gene_name")]
thres2log2fcmice_W12<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCmice_w12<(-0.5) | LOG2FCFINALcells$LOG2FCmice_w12 > 0.5), ]
thres2log2fcmice_W12<-thres2log2fcmice_W12[c("LOG2FCmice_w12","gene_name")]
thres2_d5.d0_W12<-merge(thres2log2fc_d5.d0cells,thres2log2fcmice_W12,by="gene_name")
thres2_d5.d0_W12
thres2_d5.d0_W12[-c(18,98:112), ]

#correlation#
thres2_d5.d0_W12$gene_name<-NULL
corthres2_d5.d0_W12<-cor(thres2_d5.d0_W12)
corrplot(corthres2_d5.d0_W12,method = "number")

#clustering#
#histthres2_d5.d0_W12<-thres2_d5.d0_W12[c("LOG2FCcellsD5.D0","LOG2FCmice_w12")]
#clustersthres2_d5.d0_W12 <- hclust(dist(histthres2_d5.d0_W12))
#plot(clustersthres2_d5.d0_W12) #Hierarhical Clustering

#pca
#pcathres2_d5.d0_W12 <- prcomp(thres2_d5.d0_W12)
#summary(pcathres2_d5.d0_W12)# Prints variance summary for all principal components
#biplot(pcathres2_d5.d0_W12)


##############################################################
##############    cells_dataset_3984_genes     ################
##        comparison_D10.D0(cells)~W28(mice)~human         ###  
##############################################################
names(LOG2FCFINALcells)
thres2log2fc_d10.d0cells<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCcellsD10.D0<(-0.5) | LOG2FCFINALcells$LOG2FCcellsD10.D0 > 0.5), ]
thres2log2fc_d10.d0cells<-thres2log2fc_d10.d0cells[c("LOG2FCcellsD10.D0","gene_name")]
thres2log2fcmice_W28<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FCmice_w28<(-0.5) | LOG2FCFINALcells$LOG2FCmice_w28 > 0.5), ]
thres2log2fcmice_W28<-thres2log2fcmice_W28[c("LOG2FCmice_w28","gene_name")]
thres2_d10.d0_W28<-merge(thres2log2fc_d10.d0cells,thres2log2fcmice_W28,by="gene_name")
thres2_d10.d0_W28
thres2log2fc_human<-LOG2FCFINALcells[ which(LOG2FCFINALcells$LOG2FChuman<(-0.5) | LOG2FCFINALcells$LOG2FChuman > 0.5), ]
thres2log2fc_human<-thres2log2fc_human[c("LOG2FChuman","gene_name")]
thres2_d10.d0_W28_human<-merge(thres2_d10.d0_W28,thres2log2fc_human,by="gene_name")
thres2_d10.d0_W28_human

#clustering#
clustersthres2_d10.d0_W28_human <- hclust(dist(thres2_d10.d0_W28_human[-1]))
plot(clustersthres2_d10.d0_W28_human, labels=thres2_d10.d0_W28_human$gene_name) #Hierarhical Clustering
rect.hclust(clustersthres2_d10.d0_W28_human, 10)

#correlation#
thres2_d10.d0_W28_human$gene_name<-NULL
corthres2_d10.d0_W28_human<-cor(thres2_d10.d0_W28_human)
corrplot(corthres2_d10.d0_W28_human,method = "number")


#pca#
pcathres2_d10.d0_W28_human <- prcomp(thres2_d10.d0_W28_human)
summary(pcathres2_d10.d0_W28_human)# Prints variance summary for all principal components
biplot(pcathres2_d10.d0_W28_human)

##################################
#    cells_dataset_8591_genes     #
#comparison_D2.D0(cells)~W5(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells3)
thres3log2fc_d2.d0cells<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCcells3D2.D0<(-0.5) | LOG2FCFINALcells3$LOG2FCcells3D2.D0 > 0.5), ]
thres3log2fc_d2.d0cells<-thres3log2fc_d2.d0cells[c("LOG2FCcells3D2.D0","gene_name")]
thres3log2fcmice_W5<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCmice_w5<(-0.5) | LOG2FCFINALcells3$LOG2FCmice_w5 > 0.5), ]
thres3log2fcmice_W5<-thres3log2fcmice_W5[c("LOG2FCmice_w5","gene_name")]
thres3_d2.d0_W5<-merge(thres3log2fc_d2.d0cells,thres3log2fcmice_W5,by="gene_name")
thres3_d2.d0_W5
thres3_d2.d0_W5[-c(7,14),]

#correlation#
thres3_d2.d0_W5$gene_name<-NULL
corthres3_d2.d0_W5<-cor(thres3_d2.d0_W5)
corrplot(corthres3_d2.d0_W5,method = "number")

#clustering#
histthres3_d2.d0_W5<-thres3_d2.d0_W5[c("LOG2FCcells3D2.D0","LOG2FCmice_w5")]
clustersthres3_d2.d0_W5 <- hclust(dist(histthres3_d2.d0_W5))
plot(clustersthres3_d2.d0_W5) #Hierarhical Clustering

#pca
pcathres3_d2.d0_W5 <- prcomp(thres3_d2.d0_W5)
summary(pcathres3_d2.d0_W5)# Prints variance summary for all principal components
biplot(pcathres3_d2.d0_W5)

##################################
#    cells_dataset_8591_genes     #
#comparison_D5.D0(cells)~W12(mice)#
##################################
#common_up-down_regulated _genes#
names(LOG2FCFINALcells3)
thres3log2fc_d5.d0cells<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCcells3D5.D0<(-0.5) | LOG2FCFINALcells3$LOG2FCcells3D5.D0 > 0.5), ]
thres3log2fc_d5.d0cells<-thres3log2fc_d5.d0cells[c("LOG2FCcells3D5.D0","gene_name")]
thres3log2fcmice_W12<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCmice_w12<(-0.5) | LOG2FCFINALcells3$LOG2FCmice_w12 > 0.5), ]
thres3log2fcmice_W12<-thres3log2fcmice_W12[c("LOG2FCmice_w12","gene_name")]
thres3_d5.d0_W12<-merge(thres3log2fc_d5.d0cells,thres3log2fcmice_W12,by="gene_name")
thres3_d5.d0_W12
thres3_d5.d0_W12[-c(35,71,72,73,103,104,105,106,107,108,109,110,111,112,113,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210), ]
nrow(thres3_d5.d0_W12)

#correlation#
thres3_d5.d0_W12$gene_name<-NULL
corthres3_d5.d0_W12<-cor(thres3_d5.d0_W12)
corrplot(corthres3_d5.d0_W12,method = "number")

#clustering#
histthres3_d5.d0_W12<-thres3_d5.d0_W12[c("LOG2FCcells3D5.D0","LOG2FCmice_w12")]
clustersthres3_d5.d0_W12 <- hclust(dist(histthres3_d5.d0_W12))
plot(clustersthres3_d5.d0_W12) #Hierarhical Clustering

#pca
pcathres3_d5.d0_W12 <- prcomp(thres3_d5.d0_W12)
summary(pcathres3_d5.d0_W12)# Prints variance summary for all principal components
biplot(pcathres3_d5.d0_W12)


##############################################################
##############    cells_dataset_8591_genes     ################
##        comparison_D10.D0(cells)~W25(mice)~human         ###  
##############################################################
names(LOG2FCFINALcells3)
thres3log2fc_d10.d0cells<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCcells3D10.D0<(-0.5) | LOG2FCFINALcells3$LOG2FCcells3D10.D0 > 0.5), ]
thres3log2fc_d10.d0cells<-thres3log2fc_d10.d0cells[c("LOG2FCcells3D10.D0","gene_name")]
thres3log2fcmice_W28<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FCmice_w28<(-0.5) | LOG2FCFINALcells3$LOG2FCmice_w28 > 0.5), ]
thres3log2fcmice_W28<-thres3log2fcmice_W28[c("LOG2FCmice_w28","gene_name")]
thres3_d10.d0_W28<-merge(thres3log2fc_d10.d0cells,thres3log2fcmice_W28,by="gene_name")
thres3_d10.d0_W28
thres3log2fc_human<-LOG2FCFINALcells3[ which(LOG2FCFINALcells3$LOG2FChuman<(-0.5) | LOG2FCFINALcells3$LOG2FChuman > 0.5), ]
thres3log2fc_human<-thres3log2fc_human[c("LOG2FChuman","gene_name")]
thres3_d10.d0_W28_human<-merge(thres3_d10.d0_W28,thres3log2fc_human,by="gene_name")
thres3_d10.d0_W28_human
thres3_d10.d0_W28_human[-c(12,13,14,40,41,42,91),]
nrow(thres3_d10.d0_W28_human)

#correlation#
thres3_d10.d0_W28_human$gene_name<-NULL
corthres3_d10.d0_W28_human<-cor(thres3_d10.d0_W28_human)
corrplot(corthres3_d10.d0_W28_human,method = "number")

#clustering#
histthres3_d10.d0_W28_human<-thres3_d10.d0_W28_human[c("LOG2FCcells3D10.D0","LOG2FCmice_w28","LOG2FChuman")]
clustersthres3_d5.d0_W12 <- hclust(dist(histthres3_d10.d0_W28_human))
plot(clustersthres3_d5.d0_W12) #Hierarhical Clustering


#pca#DO NOT RUN BECAUSE THERE IS ONLY ONE GENE#
pcathres3_d10.d0_W28_human <- prcomp(histthres3_d10.d0_W28_human)
summary(pcathres3_d10.d0_W28_human)# Prints variance summary for all principal components
biplot(pcathres3_d10.d0_W28_human)
