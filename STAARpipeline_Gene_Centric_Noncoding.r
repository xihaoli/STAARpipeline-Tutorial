#####################################################################
# Gene-centric analysis for noncoding rare variants of protein-coding 
# genes using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
#####################################################################

rm(list=ls())
gc()

## load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

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
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Gene_Centric_Noncoding/Results/"
## output file name
output_file_name <- "TOPMed_F5_LDL_Noncoding"

## input array id from batch file (Harvard FAS cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50 
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
group.num <- group.num.allchr[chr]

if (chr == 1){
   groupid <- arrayid
}else{
   groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num)
{ 
	sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
	sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}	

### exclude large genes
jobid_exclude <- c(21,39,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
sub_seq_id_exclude <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)
	
for(i in 1:length(jobid_exclude))
{
	if(arrayid==jobid_exclude[i])
	{
		sub_seq_id <- setdiff(sub_seq_id,sub_seq_id_exclude[i])
	}
}


### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

genes <- genes_info

results_noncoding <- c()
for(kk in sub_seq_id)
{
	print(kk)
	gene_name <- genes_info_chr[kk,1]
	results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
								QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	
	results_noncoding <- append(results_noncoding,results)
}

save(results_noncoding,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)
