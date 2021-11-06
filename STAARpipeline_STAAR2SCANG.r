##########################################################
# fit SCANG-STAAR null model
# Xihao Li, Zilin Li
# 11/04/2021
##########################################################

rm(list=ls())
gc()

### load required package
library(STAAR)
library(SCANG)
library(STAARpipeline)
library(Matrix)

###########################################################
#           User Input
###########################################################
staar_nullmodel_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"
scang_staar_nullmodel_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.SCANG.LDL.fulladj.group.size.30.20210920.Rdata"

###########################################################
#           Main Function 
###########################################################
### load STAAR null model
obj_nullmodel <- get(load(staar_nullmodel_path))

obj_nullmodel <- staar2scang_nullmodel(obj_nullmodel)

save(obj_nullmodel,file=scang_staar_nullmodel_path)

