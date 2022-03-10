###########################################################
# fit SCANG-STAAR null model
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################

rm(list=ls())
gc()

## load required packages
library(STAAR)
library(SCANG)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
staar_nullmodel_path <- "/path_to_the_file/obj_nullmodel.Rdata"
scang_staar_nullmodel_path <- "/path_to_the_output_file/obj_nullmodel_SCANG_STAAR.Rdata"

###########################################################
#           Main Function 
###########################################################
## load STAAR null model
obj_nullmodel <- get(load(staar_nullmodel_path))

obj_nullmodel_SCANG_STAAR <- staar2scang_nullmodel(obj_nullmodel)

save(obj_nullmodel_SCANG_STAAR,file=scang_staar_nullmodel_path)

