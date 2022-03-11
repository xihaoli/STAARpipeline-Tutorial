###########################################################
# Summarization and visualization of individual analysis
# results using STAARpipelineSummary
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
## Number of jobs for each chromosome
jobs_num <- get(load("/path_to_the_file/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Known loci
known_loci <- get(load("/path_to_the_file/TOPMed_F5_LDL_known_loci_genome_LD_pruning.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## results path
input_path <- "/path_to_the_results_file/"
output_path <- input_path
## results name
individual_results_name <- "TOPMed_F5_LDL_results_individual_analysis"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"
## alpha level
alpha <- 5E-08

###########################################################
#           Main Function 
###########################################################
Individual_Analysis_Results_Summary(agds_dir=agds_dir,jobs_num=jobs_num,input_path=input_path,output_path=output_path,
                                    individual_results_name=individual_results_name,
                                    obj_nullmodel=obj_nullmodel,known_loci=known_loci,
                                    method_cond=method_cond,
                                    QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,
                                    alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)

