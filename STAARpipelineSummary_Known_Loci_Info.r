###########################################################
# Extract CHR, POS, REF and ALT from #rs 
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

###########################################################
#           User Input
###########################################################
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Input GWASCatalog variants (#rs)
known_loci <- read.csv("/path_to_the_file/LDL_GWASCatalog_rsid.csv")
## rs channel name in aGDS
rs_channel <- "annotation/info/FunctionalAnnotation/rsid"
## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F5_LDL_known_loci_info"

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
		REF <- as.character(seqGetData(genofile, "$ref"))
		ALT <- as.character(seqGetData(genofile, "$alt"))
	
		known_loci_info_chr <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT)
		known_loci_info <- rbind(known_loci_info,known_loci_info_chr)	
	}
	seqClose(genofile)
}

# Output Info of GWASCatalog SNVs
save(known_loci_info,file=paste0(output_path,output_file_name,".Rdata"))

