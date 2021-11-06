##########################################################
# fit STAAR null model using GENESIS package
# Xihao Li, Zilin Li
# 11/04/2021
##########################################################

rm(list=ls())
gc()

library(GENESIS)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################

## Phenotype data
phenotype <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/LDL_20210915_group_size_30.Rdata"))
## GRM data
sgrm <- get(load("/n/holystore01/LABS/xlin/Lab/zilinli/TopMed/Lipid-Phenotype/pcrelate_kinshipMatrix_sparseDeg4_v2.RData"))
## file directory for the output file 
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/" 
## output file name
output_name <- "obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"


###########################################################
#           Main Function 
###########################################################

### fit null model using GENESIS 
data_GENESIS <- as(phenotype,"AnnotatedDataFrame") # Make AnnotatedDataFrame (specifically required by GENESIS)
obj.nullmodel.GENESIS  <- fitNullModel(data_GENESIS, outcome="LDL.norm", covars=c("age", "age2", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "het_resid_var_group_rolled"), cov.mat=sgrm, group.var="het_resid_var_group_rolled", AIREML.tol = 1e-4, verbose=TRUE)

### convert GENESIS null model to STAAR null model
obj.nullmodel <- genesis2staar_nullmodel(obj.nullmodel.GENESIS)
save(obj.nullmodel,file = paste0(output_path,output_name))






