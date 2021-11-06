##########################################################
# Pre-step for running STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
##########################################################

rm(list=ls())
gc()

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

###########################################################
#           User Input
###########################################################

## file directory of aGDS file (genotype and annotation data) 
dir.geno <- "/n/holystore01/LABS/xlin/Lab/xihaoli/TOPMed_Freeze_5/TOPMed.Anno-All-In-One-GDS-v1.1.2/"
## file name of aGDS, seperate by chr number 
adgs_file_name_1 <- "freeze.5.chr"
agds_file_name_2 <- ".pass_and_fail.gtonly.minDP0.gds"
## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/filter"
## file directory for the output files
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/" 
## annotation name. The first eight names are used to define masks in gene-centric analysis, do not change them!! 
## The others are the annotation you want to use in the STAAR procedure, and they are flexible to change.
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein","aPC.Liver")
## channel name of the annotations. Make sure they are matched with the name, especially for the first eight one!! 
dir <- c("/dbSNP_rs_num","/GENCODE.Category","/GENCODE.Info","/GENCODE.EXONIC.Category","/dbNSFP/MetaSVM_pred","/GeneHancer","/CAGE.tc","/rOCRs","/CADD.FULL/PHRED","/LINSIGHT.PHRED.rounded","/FATHMM.XF.PHRED.rounded",
					"/APC.PHRED.rounded/aPC.EpigeneticActive","/APC.PHRED.rounded/aPC.EpigeneticRepressed","/APC.PHRED.rounded/aPC.EpigeneticTranscription",
					"/APC.PHRED.rounded/aPC.Conservation.v2","/APC.PHRED.rounded/aPC.LocalDiversity.v2","/APC.PHRED.rounded/aPC.Mappability",
					"/APC.PHRED.rounded/aPC.TF","/APC.PHRED.rounded/aPC.Protein","/APC.PHRED.rounded/aPC.Liver")

###########################################################
#           Main Function 
###########################################################

#### aGDS directory
agds_dir <- paste0(dir.geno,adgs_file_name_1,seq(1,22),agds_file_name_2) 
save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

#### Annotation dir
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))

#### jobs_num
jobs_num <- matrix(rep(0,66),nrow=22)
for(chr in 1:22)
{
	print(chr)
	gds.path <- agds_dir[chr] 
	genofile <- seqOpen(gds.path)
	
	filter <- seqGetData(genofile, QC_label)
	SNVlist <- filter == "PASS" 

	position <- as.numeric(seqGetData(genofile, "position"))
	position_SNV <- position[SNVlist]
  
	jobs_num[chr,1] <- chr
	jobs_num[chr,2] <- min(position[SNVlist])
	jobs_num[chr,3] <- max(position[SNVlist])

	seqClose(genofile)
}

# Individual Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/10e6))
# Sliding Window
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/5e6))
# SCANG
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/1.5e6))

colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num","sliding_window_num","scang_num")
jobs_num <- as.data.frame(jobs_num)

save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))



