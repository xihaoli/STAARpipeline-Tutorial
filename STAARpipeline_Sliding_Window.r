#####################################################################
# Sliding window analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
#####################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Number of jobs for each chromosome
jobs_num <- get(load("/path_to_the_file/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
					
## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F5_LDL_Sliding_Window"
## input array id from batch file (Harvard FAS RC cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$sliding_window_num))
group.num <- jobs_num$sliding_window_num[chr]

if (chr == 1){
   groupid <- arrayid
}else{
   groupid <- arrayid - cumsum(jobs_num$sliding_window_num)[chr-1]
}

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

start_loc <- (groupid-1)*5e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 1000*25 - 1

results_sliding_window <- c()
for(kk in 1:200)
{
	print(kk)
	start_loc_sub <- start_loc + 1000*25*(kk-1)
	end_loc_sub <- end_loc + 1000*25*(kk-1) + 1000
	
	end_loc_sub <- min(end_loc_sub,jobs_num$end_loc[chr])
	
	results <- c()
	if(start_loc_sub < end_loc_sub)
	{
		results <- try(Sliding_Window(chr=chr,start_loc=start_loc_sub,end_loc=end_loc_sub,genofile=genofile,obj_nullmodel=obj_nullmodel,type="multiple",
		                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
		
		if(class(results)!="try-error")
		{
			results_sliding_window <- rbind(results_sliding_window,results)
		}

	}
}

save(results_sliding_window,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

