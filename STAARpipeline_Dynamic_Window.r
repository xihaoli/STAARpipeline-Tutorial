#####################################################################
# Dynamic window analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
#####################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(SCANG)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Number of jobs for each chromosome
jobs_num <- get(load("/path_to_the_file/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Null model
obj_nullmodel_SCANG_STAAR <- get(load("/path_to_the_file/obj_nullmodel_SCANG_STAAR.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F5_LDL_SCANG"
## input array id from batch file (Harvard FAS RC cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## Number of jobs for SCANG
sum(jobs_num$scang_num)

chr <- which.max(arrayid <= cumsum(jobs_num$scang_num))
group.num <- jobs_num$scang_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$scang_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*1.5e6 + jobs_num$start_loc[chr]
end_loc <- min(start_loc + 1.5e6 - 1, jobs_num$end_loc[chr])

results_scang <- Dynamic_Window_SCANG(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel_SCANG_STAAR,
                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

save(results_scang,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

