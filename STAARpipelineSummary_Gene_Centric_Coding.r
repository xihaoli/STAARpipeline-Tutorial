##########################################################
# Summarization and Visualization of Gene-centric 
# Coding Analysis Results using STAARpipelineSummary
# Xihao Li, Zilin Li
# 11/04/2021
##########################################################

rm(list=ls())
gc()

## load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(dplyr)
library(RColorBrewer)
palette("default") 

###########################################################
#           User Input
###########################################################
## job nums
jobs_num <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))

## results path
input_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Gene_Centric_Coding/Results/"
output_path <- input_path

## number of jobs
gene_centric_coding_jobs_num <- 381
## results name
gene_centric_results_name <- "TOPMed_F5_LDL_Coding"

### known loci
# known_loci <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/TOPMed_F8_LDL_known_loci_genome_LD_pruning.Rdata"))
known_loci <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/TOPMed_F8_LDL_known_loci_individual_analysis_genome_LD_pruning.Rdata"))
## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"
## variant_type
variant_type <- "SNV"
## alpha level
alpha <- 5E-07

## Annotation_dir
Annotation_dir <- "annotation/info/TOPMedAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein","aPC.Liver")

###########################################################
#           Main Function 
###########################################################
col_c <- rainbow(7)

Gene_Centric_Coding_Results_Summary(agds_dir=agds_dir,gene_centric_coding_jobs_num=gene_centric_coding_jobs_num,
input_path=input_path,output_path=output_path,gene_centric_results_name=gene_centric_results_name,
obj_nullmodel=obj_nullmodel,known_loci=known_loci,method_cond=method_cond,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
QC_label=QC_label,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)



