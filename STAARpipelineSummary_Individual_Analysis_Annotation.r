###########################################################
# Annotate a list of variants using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

###########################################################
#           User Input
###########################################################
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))

## QC_label
QC_label <- "annotation/filter"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS",
                     "CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## single variant analysis results
individual_results_sig <- get(load("/path_to_the_file/results_sig_cond.Rdata"))

###########################################################
#           Main Function 
###########################################################
results_individual_analysis_sig_anno <- Annotate_Single_Variants(agds_dir=agds_dir,single_variants_list=individual_results_sig,
                                                                 QC_label=QC_label,Annotation_dir=Annotation_dir,
                                                                 Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)

save(results_individual_analysis_sig_anno,file="/path_to_the_output_file/results_individual_analysis_sig_anno_cond.Rdata")
write.csv(results_individual_analysis_sig_anno,file="/path_to_the_output_file/results_individual_analysis_sig_anno_cond.csv",row.names=FALSE)

