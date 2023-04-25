This repository contains the scripts used to perform differential gene expression analysis for Peyton Coleman and Jonathan Lifferth's final project for BSCI 5272 (Genome Science).

Data were downloaded from the following portals:
CPTAC (cancer cohort): https://gdc.cancer.gov/ (metadata to download the exact files we used can be found at /data/cptac/metadata.card.2023-03-29.json)
GTEx (cancer cohort): https://gtexportal.org/home/datasets (we used the V8 dataset)

Cancer vs. Non-cancer differential expression
	Workflow:
		/cancer_noncancer_DESeq/making_genecounts_files.R
		/cancer_noncancer_DESeq/diff_express_pancancer_format.R
		/cancer_noncancer_DESeq/diff_express_pancancer_deseq.R


