# Takes a VCF file as input and outputs it to GDS format.
# A wrapper for the seqVCF2GDS function in the SeqArray package.
# options(repos = c(CRAN = "http://cran.revolutionanalytics.com"))
# VCF Input: can be one VCF or several, as a character vector
# Author(s): Michael R. Brown, Jennifer A. Brody

args <- commandArgs(trailingOnly = TRUE)
print(args)
import.format <- args[1] # characters, the variable name(s) in the FORMAT field for import; or NULL for all variables
fileformat <- args[2] # either 'vcf' or 'bcf'
base.filename <- args[3] # output GDS file name
num.files <- as.numeric(args[4]) # number of input VCF files
vcf.file <- args[5:(5+num.files-1)] # input VCF file names
print(vcf.file)

# GDS Ouput
gds.file <- paste0(base.filename,".gds")

library(SeqArray)

# Conversion -- This could take some time (if the input files are large)
nSlots <- system("nproc",intern=T)
nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
if (nThreads == 0) nThreads <- 1
message(paste("Running with", nThreads,"thread(s)."))
Sys.setenv(MKL_NUM_THREADS=nThreads)


if(fileformat == 'bcf'){
    ## is this a bcf file?
    ## use bcftools to read text
    message("converting BCF")
    message("BCF always single thread")
	seqBCF2GDS(vcf.file,gds.file,fmt.import=import.format, storage.option="LZMA_RA",bcftools="bcftools")
}else{
    message("converting VCF")
	seqVCF2GDS(vcf.file,gds.file,fmt.import=import.format, storage.option="LZMA_RA",parallel=nThreads)
}
# Quick summary for log's sake
g <- seqOpen(gds.file)
seqSummary(g)

