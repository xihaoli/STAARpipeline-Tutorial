##########################################################
# Summarization and Visualization of 
# Individual Analysis Results using STAARpipelineSummary
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

###########################################################
#           User Input
###########################################################
## job nums
jobs_num <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))

## results path
input_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Individual_Analysis/Results/"
output_path <- input_path
## results name
individual_results_name <- "TOPMed_F5_LDL_results_individual_analysis"

### known loci
known_loci <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/TOPMed_F5_LDL_known_loci_genome_LD_pruning.Rdata"))
## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "minor"

## method_cond
method_cond <- "optimal"

###########################################################
#           Main Function 
###########################################################

Individual_Analysis_Results_Summary(agds_dir=agds_dir,jobs_num=jobs_num,input_path=input_path,output_path=output_path,
individual_results_name=individual_results_name,obj_nullmodel=obj_nullmodel,known_loci=known_loci,method_cond=method_cond,
QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,
alpha=5E-08,manhattan_plot=TRUE,QQ_plot=TRUE)

