##########################################################
# Individual Analysis using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
##########################################################

rm(list=ls())
gc()

## load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## job nums
jobs_num <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load("//n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Individual_Analysis/Results/"
## output file name
output_file_name <- "TOPMed_F5_LDL_results_individual_analysis"
## input array id from batch file (Harvard FAS cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################

chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
   groupid <- arrayid
}else{
   groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc < end_loc)
{
	results_individual_analysis <- Individual_Analysis(chr=chr, start_loc=start_loc, end_loc=end_loc, genofile=genofile, obj_nullmodel=obj_nullmodel, variant_type=variant_type, QC_label=QC_label, geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)
