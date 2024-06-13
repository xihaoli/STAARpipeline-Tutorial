###########################################################
# fit STAAR null model using GENESIS package
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
rm(list=ls())
gc()

library(GENESIS)
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
output_name <- "obj_nullmodel_GENESIS.Rdata"

###########################################################
#           Main Function 
###########################################################
## fit null model using GENESIS 
data_GENESIS <- as(phenotype,"AnnotatedDataFrame") # Make AnnotatedDataFrame (specifically required by GENESIS)
obj_nullmodel_GENESIS <- fitNullModel(data_GENESIS,outcome="LDLadj.norm",
                                      covars=c("age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","study_ethnicity"),
                                      cov.mat=sgrm,group.var="study_ethnicity",AIREML.tol=1e-4,verbose=TRUE)

## convert GENESIS null model to STAAR null model
obj_nullmodel <- genesis2staar_nullmodel(obj_nullmodel_GENESIS)
save(obj_nullmodel,file=paste0(output_path,output_name))

