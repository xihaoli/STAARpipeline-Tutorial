###########################################################
# Independent variants selection from a list of variants
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
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
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Info of known variants
known_loci_info <- read.csv("/path_to_the_file/TOPMed_F5_LDL_Known_Loci_info.csv")
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F5_LDL_known_loci_genome_LD_pruning"

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"

## maf_cutoff
samplesize <- length(obj_nullmodel$id_include)
maf_cutoff <- 20.5/samplesize

## input chr number from batch file
chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## Variant list
known_loci_info <- known_loci_info[!is.na(known_loci_info[,1]),]
variants_list <- known_loci_info

samplesize <- length(obj_nullmodel$id_include)

### gds file
known_loci_chr <- c()

print(chr)
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

known_loci <- LD_pruning(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,variants_list=variants_list,maf_cutoff=maf_cutoff,
                         method_cond=method_cond,QC_label=QC_label,
                         variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
	
known_loci_chr <- rbind(known_loci_chr,known_loci)

seqClose(genofile)

save(known_loci_chr,file=paste0(output_path,output_file_name,"_chr",chr,".Rdata"))

