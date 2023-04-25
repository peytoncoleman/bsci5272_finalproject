### Run differential expression with DESeq2
### Author: Peyton Coleman
### Last updated: 4/25/2023

library(DESeq2)
library(data.table)
library(dplyr)

countData <- fread("counts/gtex_cptac_counts.txt")
metaData <- fread("counts/gtex_cptac_meta.txt")

# make sure they have exact same numbers
metaData <- filter(metaData, ids %in% colnames(countData))
countData <- select(countData, all_of(metaData$ids), Transcript_stable_ID_version)

# fix formatting
geneId <- countData$Transcript_stable_ID_version # make gene names into rownames
rawCounts <- as.matrix(select(countData, -c("Transcript_stable_ID_version")))  # remove gene names col
rownames(rawCounts) <- geneId  # move ids to roawnames for counts
rownames(metaData) <- metaData$ids # move ids to rownames for metadata
metaData$ids <- NULL # remove ID col
metaData$cohort <- factor(metaData$cohort, levels = c("gtex", "cptac")) # turn cohort into factor (cancer as higher level)
metaData$gender <- as.factor(metaData$gender) # gender into factor
metaData$age_group <- cut(metaData$age, breaks = seq(0, 100, 10), include.lowest = TRUE) # make age groups into factor 
metaData$age <- NULL


### normalize separately
counts_control <- fread("counts/gtex_counts.txt")
counts_cancer <- fread("counts/cptac_counts.txt")

# make separate metadatas
meta_control <- filter(metaData, cohort == "gtex")
meta_cancer <- filter(metaData, cohort == "cptac")

# extract gene IDs, make the counts into a matrix, and then add the gene IDs as rownames
geneId_control <- counts_control$transcript
rawCounts_control <- as.matrix(select(counts_control, -c("transcript")))
rownames(rawCounts_control) <- geneId_control 

geneId_cancer <- counts_cancer$Transcript_stable_ID_version
rawCounts_cancer <- as.matrix(select(counts_cancer, -c("Transcript_stable_ID_version")))
rownames(rawCounts_cancer) <- geneId_cancer 

# make sure metadata and counts have the same IDs
counts_control <- select(counts_control, all_of(meta_control$ids))
counts_cancer <- select(counts_cancer, all_of(meta_cancer$ids))

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts_control <- rawCounts_control[,unique(rownames(meta_control))]
all(colnames(rawCounts_control) == rownames(meta_control))

rawCounts_cancer <- rawCounts_cancer[,unique(rownames(meta_cancer))]
all(colnames(rawCounts_cancer) == rownames(meta_cancer))

# create DESeq dataset
dds_control <- DESeqDataSetFromMatrix(countData = round(rawCounts_control), colData = meta_control, design = ~1)
dds_cancer <- DESeqDataSetFromMatrix(countData = round(rawCounts_cancer), colData = meta_cancer, design = ~1)

# Estimate normalization factors for each cohort
dds_control <- estimateSizeFactors(dds_control)
dds_cancer <- estimateSizeFactors(dds_cancer)

# Normalize count matrices within each cohort
counts_normalized_control <- counts(dds_control, normalized = TRUE)
counts_normalized_cancer <- counts(dds_cancer, normalized = TRUE)

control <- as.data.frame(counts_normalized_control)
control$genes <- rownames(control)
control <- select(control, c(genes, everything()))

cancer <- as.data.frame(counts_normalized_cancer)
cancer$genes <- rownames(cancer)
cancer <- select(cancer, c(genes, everything()))

# MERGE IN BASH FOR SPEED
 write.table(cancer, "counts/cancer_normalized.txt", row.names = F, quote = F)
 write.table(control, "counts/control_normalized.txt", row.names = F, quote = F)

### BASH
### 
###tail -n +2 control_normalized.txt | sort -k1 > sorted_control_normalized.txt
###echo "$(head -n 1 control_normalized.txt)" | cat - sorted_control_normalized.txt > temp && mv temp sorted_control_normalized.txt

###tail -n +2 cancer_normalized.txt | sort -k1 > sorted_cancer_normalized.txt
###echo "$(head -n 1 cancer_normalized.txt)" | cat - sorted_cancer_normalized.txt > temp && mv temp sorted_cancer_normalized.txt

### join -1 1 -2 1 sorted_control_normalized.txt sorted_cancer_normalized.txt > counts_normalized_separately.txt

# read in the merged counts
sep <- fread("counts/counts_normalized_separately.txt")
rownames(sep) <- sep$genes
sep$genes <- NULL

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
sep <- sep[,unique(rownames(metaData))]
all(colnames(sep) == rownames(metaData))
sep <- as.matrix(sep)

# create DESeq object
dds <- DESeqDataSetFromMatrix(countData=round(sep), 
                              colData=metaData, 
                              design= (~ gender + age_group + cohort))

# run DESeq
data <- DESeq(dds)

print(results(data))

# save 
save(dds, data, counts_normalized, file = "deseq.RDat")
write.table(rownames(counts_normalized), "gene_order.txt")
