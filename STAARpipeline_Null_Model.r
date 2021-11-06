##########################################################
# fit STAAR null model
# Xihao Li, Zilin Li
# 11/04/2021
##########################################################

rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################

## Phenotype data
phenotype <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/UKB_WES_lipids/tmp.LDL.20211014.Rdata"))
## GRM data
sgrm <- get(load("/n/holystore01/LABS/xlin/Lab/rdey/UKB_Analysis/FastSparseGRM/grmout/grmout_nRandomSNPs_0.sGRM.RData"))
## file directory for the output file 
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/UKB_WES_lipids/" 
## output file name
output_name <- "obj.STAAR.UKB.LDL.20211014.Rdata"

###########################################################
#           Main Function 
###########################################################

### fit null model
obj.nullmodel <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype, kins = sgrm, kins_cutoff = 0.022, id = "userId", use_sparse = TRUE,family = gaussian(link = "identity"), verbose=T)

save(obj.nullmodel,file = paste0(output_path,output_name))








