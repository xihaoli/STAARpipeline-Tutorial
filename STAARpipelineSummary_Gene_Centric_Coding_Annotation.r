##########################################################
# Annotate Rare Variants in Coding Masks
# using STAARpipelineSummary
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
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

###########################################################
#           User Input
###########################################################
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## results path
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Specific_Region/"

### known loci
# known_loci <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/TOPMed_F5_LDL_known_loci_genome_LD_pruning.Rdata"))
known_loci <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/TOPMed_F5_LDL_known_loci_individual_analysis_genome_LD_pruning.Rdata"))
## Null Model
obj_nullmodel <- get(load("//n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/obj.GENESIS.STAAR.LDL.fulladj.group.size.30.20210915.Rdata"))

## Annotation_dir
Annotation_dir <- "annotation/info/TOPMedAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM",
					"GeneHancer","CAGE","DHS",
					"CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## Parameter
QC_label <- "annotation/filter"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	
method_cond <- "optimal"				

## Chr
chr_seq <- c(15,17,17,17)
## Gene name
gene_name_seq <- c("ACAN","GH1","CSHL1","CD79B")
## Coding mask
category_seq <- c("missense","disruptive_missense","disruptive_missense","synonymous")

###########################################################
#           Main Function 
###########################################################
for(kk in 1:length(chr_seq))
{

	chr <- chr_seq[kk]
	gene_name <- gene_name_seq[kk]
	category <- category_seq[kk]
	
	print(gene_name)

	### gds file
	gds.path <- agds_dir[chr]
	genofile <- seqOpen(gds.path)

	results_info <- Gene_Centric_Coding_Info(category=category,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,
	QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
	Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)

	seqClose(genofile)

	save(results_info,file=paste0(output_path,gene_name,"_",category,".Rdata"))
	write.csv(results_info,paste0(output_path,gene_name,"_",category,".csv"))
}


