###########################################################
# Summarization and visualization of sliding window
# analysis results using STAARpipelineSummary
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
### Known loci
known_loci <- get(load("/path_to_the_file/TOPMed_F5_LDL_known_loci_individual_analysis_genome_LD_pruning.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## results path
input_path <- "/path_to_the_results_file/"
output_path <- input_path
## results name
sliding_window_results_name <- "TOPMed_F5_LDL_Sliding_Window"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"
## variant_type
variant_type <- "SNV"
## alpha level
alpha <- 0.05

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

###########################################################
#           Main Function 
###########################################################
Sliding_Window_Results_Summary(agds_dir=agds_dir,jobs_num=jobs_num,
                               input_path=input_path,output_path=output_path,sliding_window_results_name=sliding_window_results_name,
                               obj_nullmodel=obj_nullmodel,known_loci=known_loci,
                               method_cond=method_cond,
                               QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                               alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)

