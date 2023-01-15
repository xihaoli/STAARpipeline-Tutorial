###########################################################
# Pre-step for running STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/21/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

###########################################################
#           User Input
###########################################################
## file directory of aGDS file (genotype and annotation data) 
dir.geno <- "/path_to_the_aGDS_file/"
## file name of aGDS, seperate by chr number 
agds_file_name_1 <- "freeze.5.chr"
agds_file_name_2 <- ".pass_and_fail.gtonly.minDP0.gds"
## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/filter"
## file directory for the output files
output_path <- "/path_to_the_output_file/" 
## annotation name. The first eight names are used to define masks in gene-centric analysis, do not change them!! 
## The others are the annotation you want to use in the STAAR procedure, and they are flexible to change.
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
          "aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
## channel name of the annotations. Make sure they are matched with the name, especially for the first eight one!! 
dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info",
         "/genecode_comprehensive_exonic_category","/metasvm_pred",
         "/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")

###########################################################
#           Main Function 
###########################################################
## aGDS directory
agds_dir <- paste0(dir.geno,agds_file_name_1,seq(1,22),agds_file_name_2) 
save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

## Annotation name catalog (alternatively, can skip this part by providing Annotation_name_catalog.csv with the same information)
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))

## Number of jobs for each chromosome
jobs_num <- matrix(rep(0,66),nrow=22)
for(chr in 1:22)
{
	print(chr)
	gds.path <- agds_dir[chr] 
	genofile <- seqOpen(gds.path)
	
	filter <- seqGetData(genofile, QC_label)
	SNVlist <- filter == "PASS" 

	position <- as.numeric(seqGetData(genofile, "position"))

	jobs_num[chr,1] <- chr
	jobs_num[chr,2] <- min(position[SNVlist])
	jobs_num[chr,3] <- max(position[SNVlist])

	seqClose(genofile)
}

## Individual Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/10e6))
## Sliding Window Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/5e6))
## Dynamic Window Analysis (SCANG-STAAR)
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/1.5e6))

colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num","sliding_window_num","scang_num")
jobs_num <- as.data.frame(jobs_num)

save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))

