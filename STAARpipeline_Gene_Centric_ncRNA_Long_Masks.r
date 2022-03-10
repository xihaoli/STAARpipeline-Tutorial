#####################################################################
# Gene-centric analysis for noncoding rare variants in long masks 
# of ncRNA genes using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
#####################################################################

rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Number of jobs for each chromosome
jobs_num <- get(load("/path_to_the_file/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

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
output_file_name <- "TOPMed_F5_LDL_ncRNA"

###########################################################
#           Main Function 
###########################################################
## gene number in job
arrayid <- c(117,218,220,220,221,156,219)
sub_seq_id <- c(53,19,208,274,311,41,103)

region_spec <- data.frame(arrayid,sub_seq_id)

gene_num_in_array <- 100 
group.num.allchr <- ceiling(table(ncRNA_gene[,1])/gene_num_in_array)
sum(group.num.allchr)

results_ncRNA <- c()
for(kk in 1:dim(region_spec)[1])
{
	arrayid <- region_spec$arrayid[kk]
	sub_seq_id <- region_spec$sub_seq_id[kk]
	
	chr <- which.max(arrayid <= cumsum(group.num.allchr))
	ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
	
	gds.path <- agds_dir[chr]
	genofile <- seqOpen(gds.path)
	
	gene_name <- ncRNA_gene_chr[sub_seq_id,2]
	results <- c()
	results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
	                     rare_maf_cutoff=0.01,rv_num_cutoff=2,
	                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
	                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
	                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
	
	results_ncRNA <- rbind(results_ncRNA,results)
	
	seqClose(genofile)
}

save(results_ncRNA,file=paste0(output_path,output_file_name,"_",223,".Rdata"))

