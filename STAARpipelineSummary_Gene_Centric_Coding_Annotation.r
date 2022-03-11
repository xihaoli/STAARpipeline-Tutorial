###########################################################
# Annotate rare variants in coding masks
# using STAARpipelineSummary
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
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Known loci
known_loci <- get(load("/path_to_the_file/TOPMed_F5_LDL_known_loci_individual_analysis_genome_LD_pruning.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## output path
output_path <- "/path_to_the_output_file/"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## variant_type
variant_type <- "SNV"
# method_cond
method_cond <- "optimal"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS",
                     "CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## Chr
chr_seq <- c(1,7,19,19,9)
## Gene name
gene_name_seq <- c("PCSK9","NPC1L1","LDLR","APOE","RNF20")
## Coding mask
category_seq <- c("plof","missense","missense","missense","synonymous")

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
	write.csv(results_info,paste0(output_path,gene_name,"_",category,".csv"),row.names=FALSE)
}

