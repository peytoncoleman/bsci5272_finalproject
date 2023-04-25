## Script to format CPTAC and GTEX data for differential expression analysis
## Author: Peyton Coleman
## Last updated: 4/25/2023

# set working directory
setwd("/home/colemapl/peypeyjonjon_bsci5272")

# load necessary libraries
library(data.table)
library(dplyr)
library(stringr)

# read in raw counts
gtex <- fread("./gtex/quantification_gencode.counts.txt")
cptac <- fread("./cptac/gene_counts_combined/all_cancers_gene_counts.txt")

# map cptac ENSG names to ENST to compare to GTEX
map <- fread("ENSG_to_ENST.txt")
colnames(map) <- gsub(" ", "_", colnames(map)) # remove spaces from map column names

map_versions <- select(map, c("Gene_stable_ID_version", "Transcript_stable_ID_version")) # make df with gene and transcript names for merge
cptac_g <- merge(map_versions, cptac, by.x = "Gene_stable_ID_version", by.y = "gene_id") # merge CPTAC
cptac_g <- select(cptac_g, -c("Gene_stable_ID_version")) # remove gene name from final DF

# see how many genes overlap between samples
gtex_in_cptac <- filter(gtex, transcript %in% cptac_g$Transcript_stable_ID_version)
## 92098 overlapping transcripts! 

# merge
alldat <- merge(cptac_g, gtex, by.x = "Transcript_stable_ID_version", by.y = "transcript")
write.table(alldat, "./counts/gtex_cptac_counts.txt", row.names = F, quote = F)
write.table(cptac_g, "/nobackup/user/colemapl/counts/cptac_counts.txt", row.names = F, quote = F)
write.table(gtex, "/nobackup/user/colemapl/counts/gtex_counts.txt", row.names = F, quote = F)

# make metadata

### CPTAC
cptac_dems_brain <- fread("cptac/clinical/brain/clinical.tsv")
cptac_dems_lung <- fread("cptac/clinical/lung/clinical.tsv")
cptac_dems_pancreas <- fread("cptac/clinical/pancreas/clinical.tsv")

# merge cptac dems together
cptac_dems <- rbind(cptac_dems_brain, cptac_dems_lung) 
cptac_dems <- rbind(cptac_dems, cptac_dems_pancreas)

# calculate age
    for(i in 1:nrow(cptac_dems)){
        if(cptac_dems$year_of_death[i] == "'--"){
            cptac_dems$age[i] <- 2023-cptac_dems$year_of_birth[i]
        } else {
            cptac_dems$age[i] <- as.numeric(cptac_dems$year_of_death[i])-cptac_dems$year_of_birth[i]
        }
    }

cptac_dems_merge <- select(cptac_dems, c("case_submitter_id", "age", "gender"))

# read in sample file to map filenames to IDs
sample <- fread("cptac/gdc_sample_sheet.2023-03-29.tsv")
names(sample)<-str_replace_all(names(sample), c(" " = "_"))
sample <- select(sample, "Case_ID", "File_Name")

# merge IDs with dems
cptac_dems_merge_names <- merge(sample, cptac_dems_merge, by.x = "Case_ID", by.y = "case_submitter_id")
cptac_dems_merge_names$ids <- sub("\\.rna_seq\\.augmented_star_gene_counts.tsv", "", basename(cptac_dems_merge_names$File_Name)) # get the basename of the file (corresponds to ID in the dems)

# get tissue
brain <- merge(sample, cptac_dems_brain, by.x = "Case_ID", by.y = "case_submitter_id")
pancreas <- merge(sample, cptac_dems_pancreas, by.x = "Case_ID", by.y = "case_submitter_id")
lung <- merge(sample, cptac_dems_lung, by.x = "Case_ID", by.y = "case_submitter_id")

brain$filestem <- str_extract(brain$File_Name, "^[^.]+")
cptac_dems_merge_names$tissue[cptac_dems_merge_names$ids %in% brain$filestem] <- "brain"

pancreas$filestem <- str_extract(pancreas$File_Name, "^[^.]+")
cptac_dems_merge_names$tissue[cptac_dems_merge_names$ids %in% pancreas$filestem] <- "pancreas"

lung$filestem <- str_extract(lung$File_Name, "^[^.]+")
cptac_dems_merge_names$tissue[cptac_dems_merge_names$ids %in% lung$filestem] <- "lung"

# get final metadata 
cptac_dems_merge_out <- select(cptac_dems_merge_names, c("ids", "age", "gender", "tissue"))
cptac_dems_merge_out$cohort <- "cptac"


### GTEX
gtex_dems <- fread("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
gtex_full_ids <- data.frame(ids = colnames(gtex)[-1])

# get partial GTEX IDs
gtex_full_ids$partial <- sub("^([^-]+-[^-]+)-.*$", "\\1", gtex_full_ids$ids)
gtex_full_ids_all <- merge(gtex_full_ids, gtex_dems, by.x = "partial", by.y = "SUBJID")

# merge our samples with dems
gtex_dems_merge <- gtex_full_ids_all %>%
    select("ids", "AGE", "SEX") %>%
    rename(c(
        "gender" = "SEX",
        "age" = "AGE"
        )
    )
gtex_dems_merge$cohort = "gtex"

## recode gender
gtex_dems_merge$gender[gtex_dems_merge$gender == 1] <- "male"
gtex_dems_merge$gender[gtex_dems_merge$gender == 2] <- "female"

## recode age
gtex_dems_merge$age[gtex_dems_merge$age == "20-29"] <- 25
gtex_dems_merge$age[gtex_dems_merge$age == "30-39"] <- 35
gtex_dems_merge$age[gtex_dems_merge$age == "40-49"] <- 45
gtex_dems_merge$age[gtex_dems_merge$age == "50-59"] <- 55
gtex_dems_merge$age[gtex_dems_merge$age == "60-69"] <- 65
gtex_dems_merge$age[gtex_dems_merge$age == "70-79"] <- 75

# rbind
meta <- rbind(cptac_dems_merge_out, gtex_dems_merge)
write.table(meta, "counts/gtex_cptac_meta.txt", row.names=F, quote=F)


