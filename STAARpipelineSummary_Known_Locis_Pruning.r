##########################################################
# Independent variants selection from a list of variants
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

###########################################################
#           User Input
###########################################################
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## Info of known variants
known_loci_info <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/TOPMed_F5_LDL_Known_Loci_info.Rdata"))
## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"))

## output path
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/"
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

## input chr number from batch file (Harvard FAS cluster)
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

known_loci <- LD_pruning(chr=chr, genofile=genofile, obj_nullmodel=obj_nullmodel, variants_list=variants_list, 
					variant_type=variant_type, geno_missing_imputation=geno_missing_imputation, method_cond=method_cond, QC_label=QC_label,
					maf_cutoff=maf_cutoff)
	
known_loci_chr <- rbind(known_loci_chr,known_loci)

seqClose(genofile)

save(known_loci_chr,file=paste0(output_path,output_file_name,"_chr",chr,".Rdata"))
