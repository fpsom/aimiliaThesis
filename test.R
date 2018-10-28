library(tidyverse)
library(dplyr)

pdata<-E_MTAB_4222 #name the table ?????
dim(pdata) #length of row names
str(pdata) #display the internal structure
head(pdata) #returns the first parts of the table
tail(pdata) #returns the last parts of the table
summary(pdata) #summary of data
names(pdata) #names of rows of the table
mydata<-select(pdata,gene_id)#create subset
mydata
mydata<-select(pdata,gene_id,FPKM.NIL38_F12_R1_E-MTAB-4222)#????????????????