This repository contains the scripts used to perform differential gene expression analysis for Peyton Coleman and Jonathan Lifferth's final project for BSCI 5272 (Genome Science).

Data were downloaded from the following portals:
CPTAC (cancer cohort): https://gdc.cancer.gov/ (metadata to download the exact files we used can be found at /data/cptac/metadata.card.2023-03-29.json)
GTEx (cancer cohort): https://gtexportal.org/home/datasets (we used the V8 dataset)

Cancer vs. Non-cancer differential expression
	Workflow:
		/cancer_noncancer_DESeq/making_genecounts_files.R 
			Makes raw gene count matrices for GTEx and CPTAC samples
		/cancer_noncancer_DESeq/diff_express_pancancer_format.R
			Formats data correctly for differential expression analysis in DESeq2
		/cancer_noncancer_DESeq/diff_express_pancancer_deseq.R
			Runs DESeq2
		/cancer_noncancer_DESeq/final_plots.Rmd
			Creates volcano plot (Figure 1)


Lung vs. Brain vs. Pancreatic Cancer differential expression
	Workflow:
		/cancer_type_ANOVA/normalize_cols.ipynb
			Normalizes data (input is output of /cancer_noncancer_DESeq/making_genecounts_files.R)
		/cancer_type_ANOVA/
			Runs ANOVA and creates manhattan plot (Figure 2)
