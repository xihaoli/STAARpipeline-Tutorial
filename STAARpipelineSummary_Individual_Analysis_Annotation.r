##########################################################
# Annotate a List of Variants using STAARpipelineSummary
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
library(STAARpipelineSummary)

###########################################################
#           User Input
###########################################################
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## QC_label
QC_label <- "annotation/filter"

## Annotation_dir
Annotation_dir <- "annotation/info/TOPMedAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## single variant analysis results
individual_results_sig <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Individual_Analysis/Results/results_sig_cond.Rdata"))

###########################################################
#           Main Function 
###########################################################

results_individual_analysis_sig_anno <- Annotate_Single_Variants(agds_dir=agds_dir,single_variants_list=individual_results_sig,
QC_label=QC_label,Annotation_dir=Annotation_dir,Annotation_name=Annotation_name,Annotation_name_catalog=Annotation_name_catalog)

save(results_individual_analysis_sig_anno,file="/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Individual_Analysis/Results/results_individual_analysis_sig_anno_cond.Rdata")
write.csv(results_individual_analysis_sig_anno,file="/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Individual_Analysis/Results/results_individual_analysis_sig_anno_cond.csv")

