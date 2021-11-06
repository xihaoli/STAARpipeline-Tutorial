##################################################################################
# Gene-centric analysis for coding rare variants in long masks using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
##################################################################################

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
## job nums
jobs_num <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/TOPMedAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein","aPC.Liver")

## output path
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Gene_Centric_Coding/Results/"
## output file name
output_file_name <- "TOPMed_F5_LDL_Coding"

## input array id from batch file (Harvard FAS cluster)
arrayid_longmask <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50 
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

### exclude large genes
arrayid <- c(57,112,112,113,113,113,113,113,113,113)
sub_seq_id <- c(840,543,544,575,576,577,578,579,580,582)

region_spec <- data.frame(arrayid,sub_seq_id) 
sub_seq_id <- ((arrayid_longmask-1)*5+1):min(arrayid_longmask*5,length(arrayid))

### gds file
genes <- genes_info

results_coding <- c()
for(kk in sub_seq_id)
{
	print(kk)
	arrayid <- region_spec$arrayid[kk]
	sub_id <- region_spec$sub_seq_id[kk]
	
	chr <- which.max(arrayid <= cumsum(group.num.allchr))
	gds.path <- agds_dir[chr]
	genofile <- seqOpen(gds.path)
	
	genes_info_chr <- genes_info[genes_info[,2]==chr,]
	gene_name <- genes_info_chr[sub_id,1]

	results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
								QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	
	results_coding <- append(results_coding,results)

	seqClose(genofile)
}

save(results_coding,file=paste0(output_path,output_file_name,"_",arrayid_longmask+379,".Rdata"))


