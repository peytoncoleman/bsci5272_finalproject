### Script to format individual counts files from CPTAC into a large matrix
### Author: Peyton Coleman
### Last updated: 4/25/2023

library(data.table)
library(dplyr)
library(stringr)
setwd("/home/colemapl/peypeyjonjon_bsci5272/cptac")
files <- fread("./gene_counts_ind/sample.tsv")

## gdc sample sheet has the submitter ID (our file name) and the sample ID which is in the exposure.tsv file which is then in case_set stuff
sample <- fread("./gdc_sample_sheet.2023-03-29.tsv")
names(sample)<-str_replace_all(names(sample), c(" " = "_"))

make_combined_df <- function(pheno) {
    clin <- fread(paste("./clinical/", pheno, "/clinical.tsv", sep = "")) # read in clinical data
    clin_sample <- dplyr::filter(sample, Case_ID %in% clin$case_submitter_id) # subset the sample file to only those who are in the pheno's clinical file
    filelist <- clin_sample$File_Name # get the filenames for our samples
    parent_path <- "./gene_counts_ind" # set the path to the parent folder
    pattern <- paste0("^", paste(filelist, collapse = "|"), "$") # create a regular expression pattern that matches any of the filenames in file_list
    all_files <- list.files(parent_path, recursive = TRUE, full.names = TRUE) # use list.files() to get a list of all files in parent and its subdirectories
    file_list <- grep(pattern, all_files, value = TRUE) # filter the list to include only the files that match the pattern
    file_data <- lapply(file_list, function(file) {
        # Extract the file name without the extension
        name <- sub("\\.rna_seq\\.augmented_star_gene_counts.tsv", "", basename(file))
        # Read the file using fread and select only the 1st and 7th columns
        fread(file, select = c(1, 7), col.names = c("gene_id", name))
        })
    
    # use Reduce() and left_join() to merge the dataframes into one dataframe
    merged_data <- file_data %>%
        Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="gene_id"), .)

    # remove the first four rows of metadata
    merged_data <- merged_data[-c(1:4),]
    # save the data
    write.table(merged_data, paste("./gene_counts_combined/", pheno, "_gene_counts.txt", sep = ""), quote = F, row.names = F)
}


make_combined_df("brain")
make_combined_df("lung")
make_combined_df("pancreas")


## combine all three cancer types
brain <- fread("./gene_counts_combined/brain_gene_counts.txt")
lung <- fread("./gene_counts_combined/lung_gene_counts.txt")
pancreas <- fread("./gene_counts_combined/pancreas_gene_counts.txt")

brain_lung <- merge(brain, lung, by = "gene_id")
all_cancer <- merge(brain_lung, pancreas, by = "gene_id")

write.table(all_cancer, "./gene_counts_combined/all_cancers_gene_counts.txt", row.names = F, quote = F)