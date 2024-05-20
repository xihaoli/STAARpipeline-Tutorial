# Load necessary libraries
#library(SeqArray)
#library(data.table)
#library(tidyr)
## Make sure that bcftools is installed on your system !

###########################################################
#           processVariantFiles Function
###########################################################
# To address memory constraints, this function converts and/or combines files into an R data object.
# (Can take in a vcf/vcf.gz, csv, txt, or Rdata files)
# It selectively retrieves data that:
# - Corresponds to the chromosome specified in the aGDS file,
# - Has matching genomic positions within the aGDS file,
# - Meets the criteria based on the specified variant type.
#
# Note: This function does not find exact matches within the aGDS file but is used
# for reading, filtering, and merging variant files while lowering overall memory consumption

########################################################################################################################
# Define a function to process a file
processVariantFiles <- function(files, chr, genofile, variant_type = 'variant') {
  # Retrieve unique genomic positions from the SeqArray file to optimize processing.
  unique_positions <- unique(as.numeric(SeqArray::seqGetData(genofile, "position")))
  
  # Initialize an empty list to store data frames from each file.
  list_data<- list()
  
  # Loop through each file in the list of files.
  for (file in files) {
    
    # Check the file type and construct appropriate command for reading data.
    if (endsWith(tolower(file), "vcf") || endsWith(tolower(file), "vcf.gz")) {
      # Read VCF files using bcftools, focusing on the specified chromosome.
      fread_cmd <- paste0("bcftools view -H -r ", chr, " ", file)
      data<- suppressWarnings(data.table::fread(cmd = fread_cmd, header = FALSE, select = c(1, 2, 4, 5)))
      
      # Retry reading the data with a modified chromosome label if the initial attempt is unsuccessful.
      if(length(data) == 0){
        fread_cmd <- paste0("bcftools view -H -r chr", chr, " ", file)
        data<- suppressWarnings(data.table::fread(cmd = fread_cmd, header = FALSE, select = c(1, 2, 4, 5)))
      }
    } else if (endsWith(tolower(file), "csv") || endsWith(tolower(file), "txt")) {
      # Handle CSV or text files.
      fread_cmd <- paste0("cat ", file)
      data<- suppressWarnings(data.table::fread(cmd = fread_cmd, header = FALSE, select = c(1:4)))
      data<- data[tolower(data[[1]]) %in% c(chr, paste0('chr', chr)), c(1:4)]
    } else if (endsWith(tolower(file), "rdata")) {
      # Load RData files and filter data based on the chromosome.
      data<- get(load(file))
      data<- data[tolower(data[[1]]) %in% c(chr, paste0('chr', chr)), c(1:4)]
    } else {
      # Stop execution and return an error message if the file format is not supported.
      stop("Error: File format is not supported!")
    }
    
    # Filter rows to include only those with positions found in the unique positions list.
    data<- data[data[[2]] %in% unique_positions,]
    
    # Standardize chromosome labeling and separate rows by alternate alleles.
    data[[1]] <- sub("^chr", "", data[[1]])
    colnames(data) <- c("CHR", "POS", "REF", "ALT")
    data<- tidyr::separate_rows(data, ALT, sep = ",")
    
    # Conditionally filter data based on the type of genetic variant.
    if (tolower(variant_type) == 'snv') {
      data<- data[nchar(data$REF) == 1 & nchar(data$ALT) == 1,]
    } else if (tolower(variant_type) == 'indel') {
      data<- data[nchar(data$REF) > 1 | nchar(data$ALT) > 1,]
    } else if (tolower(variant_type) != 'variant') {
      stop("Error: Please choose variant_type from list: 'SNV', 'Indel', or 'variant'")
    }
    
    # Add the processed data frame to the list.
    list_data[[file]] <- data
  }
  
  # Combine all data frames into a single data frame and remove duplicates.
  combined_data<- unique(do.call(rbind, list_data))
  
  # Return the data frame, sorted by the position column.
  return(combined_data[order(as.numeric(combined_data$POS)), ])
}

## Create Variant list Example
# Call the function to process the list of files with specified parameters
# files: Full path to a file or a list of full paths to each file. (Can take in a vcf/vcf.gz, csv, txt, or Rdata files)
#       non-VCF files must have first 4 columns in order of CHR, POS, REF, and ALT
# chr: Chromosome value corresponding to the aGDS file.
# genofile: Opened aGDS file (opened with seqOpen()).
# variant_type: Type of variant to filter within the file(s). Options include 'SNV', 'Indel', or 'variant'. Default is "variant".
# variants_list <- processVariantFiles(files, chr=chr, genofile=genofile, variant_type=variant_type)
