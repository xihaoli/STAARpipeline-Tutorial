###########################################################
# fit STAAR null model
# Xihao Li, Zilin Li
# 11/04/2021
###########################################################

rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################

## Phenotype file
phenotype <- read.csv("/path_to_the_file/pheno.csv")
## (sparse) GRM file
sgrm <- get(load("/path_to_the_file/sGRM.Rdata"))
## file directory for the output file 
output_path <- "/path_to_the_output_file/" 
## output file name
output_name <- "obj_STAAR_nullmodel.Rdata"

###########################################################
#           Main Function 
###########################################################

### fit null model
obj.nullmodel <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype, kins = sgrm, kins_cutoff = 0.022, id = "userId", use_sparse = TRUE, family = gaussian(link = "identity"), verbose=TRUE)

save(obj.nullmodel,file = paste0(output_path,output_name))

