##########################################################
# Adds QC_label with all "PASS" to a post-QC GDS file
# Authors: Xihao Li, Zilin Li
##########################################################

library(gdsfmt)
library(SeqArray)

dir_geno <- "/path_to_the_GDS_file/"
gds_file_name_1 <- "freeze.5.chr"
gds_file_name_2 <- ".pass_and_fail.gtonly.minDP0.gds"

for (chr in 1:22){
  print(paste("Chromosome:",chr))
  gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
  genofile<-seqOpen(gds.path, readonly = FALSE)
  #genofile
  
  position <- as.integer(seqGetData(genofile, "position"))
  length(position)
  Anno.folder <- index.gdsn(genofile, "annotation/info")
  add.gdsn(Anno.folder, "QC_label", val=factor(rep("PASS", length(position))), compress="LZMA_ra", closezip=TRUE)
  #genofile
  
  seqClose(genofile)
}

