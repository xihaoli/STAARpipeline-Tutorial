##########################################################
# Get Chr, Pos, REF and ALT from #rs 
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

###########################################################
#           User Input
###########################################################
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/agds_dir.Rdata"))
## Input GWASCatalog variants (#rs)
known_loci <- read.csv("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/LDL_GWASCatalog_rsid.csv",header=TRUE)
## rs channel name in aGDS
rs_channel <- "annotation/info/TOPMedAnnotation/dbSNP_rs_num"
## output path
output_path <- "/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_LDL/Known_Loci/"
## output file name
output_file_name <- "TOPMed_F5_LDL_Known_Loci_info"

###########################################################
#           Main Function 
###########################################################
known_loci_info <- c()

for(chr in 1:22)
{
	print(chr)
	
	print(paste("Chromosome:",chr))
	gds.path <- agds_dir[chr]
	genofile <- seqOpen(gds.path)
		
	variant.id <- seqGetData(genofile, "variant.id")
	rs_num <- seqGetData(genofile,rs_channel)

	rs_num_in <- rs_num%in%known_loci$rs
	
	if(sum(rs_num_in) > 0)
	{
		variant.id.in <- variant.id[rs_num_in]
		
		rm(rs_num)
		gc()
	
		rm(variant.id)
		gc()

		seqSetFilter(genofile,variant.id=variant.id.in)
	
		### Basic Info of Significant Loci
		position <- as.numeric(seqGetData(genofile, "position"))
		REF <- unlist(lapply(strsplit(seqGetData(genofile, "allele"),","),`[[`,1))
		ALT <- unlist(lapply(strsplit(seqGetData(genofile, "allele"),","),`[[`,2))
	
		known_loci_info_chr <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT)
		known_loci_info <- rbind(known_loci_info,known_loci_info_chr)	
	}
	seqClose(genofile)
}

# Output Info of GWASCatalog SNVs
save(known_loci_info,file=paste0(output_path,output_file_name,".Rdata"))
