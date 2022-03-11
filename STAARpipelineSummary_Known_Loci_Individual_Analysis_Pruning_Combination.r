###########################################################
# Combine chromosome-wide results into genome-wide   
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
rm(list=ls())
gc()

###########################################################
#           User Input
###########################################################
## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F5_LDL_known_loci_individual_analysis_genome_LD_pruning"

###########################################################
#           Main Function 
###########################################################
known_loci_genome <- c()
for(chr in 1:22)
{
	load(paste0(output_path,output_file_name,"_chr",chr,".Rdata"))
	known_loci_genome <- rbind(known_loci_genome,known_loci_chr)
}

save(known_loci_genome,file=paste0(output_path,output_file_name,".Rdata"))

