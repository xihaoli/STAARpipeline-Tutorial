#####################################################################
# Dynamic window analysis using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
#####################################################################

rm(list=ls())
gc()

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(SCANG)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## job nums
jobs_num <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.SCANG.LDL.fulladj.group.size.30.20210920.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"


## Annotation_dir
Annotation_dir <- "annotation/info/TOPMedAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_Height/Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein","aPC.Liver")


## output path
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/SCANG/Results/"
## output file name
output_file_name <- "results_scang"
## input array id from batch file (Harvard FAS cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])


###########################################################
#           Main Function 
###########################################################
## job number of scang
sum(jobs_num$scang_num)

chr <- which.max(arrayid <= cumsum(jobs_num$scang_num))
group.num <- jobs_num$scang_num[chr]

if (chr == 1){
   groupid <- arrayid
}else{
   groupid <- arrayid - cumsum(jobs_num$scang_num)[chr-1]
}

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

start_loc <- (groupid-1)*1.5e6 + jobs_num$start_loc[chr]
end_loc <- min(start_loc + 1.5e6 - 1, jobs_num$end_loc[chr])

a <- Sys.time()
results_scang <- Dynamic_Window_SCANG(chr=chr, start_loc=start_loc, end_loc=end_loc, genofile=genofile, obj_nullmodel=obj_nullmodel,
					QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
					Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
					Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
b <- Sys.time()

b - a

save(results_scang,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)
