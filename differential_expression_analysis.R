# Differential expression of genes between risk and protective PLCG2 variants in Alzheimer's disease

#_ Set up of environment ---------- 
install.packages('R.utils')
install.packages('tidyverse')
library(R.utils)
library(tidyverse)


#_ Load data file ---------- 
data_raw <- read_csv('../data/GSE308813_20240829_LB_rna_counts_only.csv')

#_ Data Exploration and Modification---------- 

print(names(data_raw))
# Verify data structure and sample names.

new_names <- c('gene','WT1','WT2','WT3','WT4',
               'KO1','KO2','KO3','KO4',
               'risk1','risk2','risk3','risk4',
               'prot1','prot2','prot3','prot4')

colnames(data_raw) <- new_names
# Rename columns for clarity

gene_names <- data_raw$gene
# Preserve gene names before filtering low-expressed genes.

data_subset_numeric <- data_raw %>% select(-gene)
# Separate count data from gene identifiers for numerical processing.

row_sums_vector <- rowSums(data_subset_numeric)
# Compute total expression per gene to identify lowly expressed genes.

keep_rows <- row_sums_vector > 100
# Define filter threshold to remove low-expression genes to reduce noise.

data_filtered <- data_subset_numeric[keep_rows, ]
filtered_gene_names <- gene_names[keep_rows]
# Apply filtering to both counts and gene names for consistency.

filtered_gene_names_unique <- make.unique(filtered_gene_names)
# Ensure unique gene identifiers to prevent downstream errors.

data_filtered <- as.data.frame(data_filtered)
rownames(data_filtered) <- filtered_gene_names_unique
# Prepare data frame with gene names as row names for DESeq2 compatibility.

metadata <- data.frame(
  sample = colnames(data_filtered),
  condition = c(rep('WT',4), rep('KO',4), rep('risk',4), rep('prot',4))
)
# Define experimental groups for differential expression analysis.

print(metadata) #data frame
print(head(data_filtered))  #matrix
# Verify metadata and data integrity post-filtering.

print(head(rownames(data_filtered)))
# Confirm gene naming correctness.

#_ Statistical analysis using DESeq2 ---------- 

#Load necessary packages
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

library(DESeq2)
# Load DESeq2 for differential expression statistical testing.

metadata$condition<-factor (metadata$condition, levels=c('prot','risk')
# Define conditions as factors with proper reference levels for analysis.

ddc <- DESeqDataSetFromMatrix(countData=data_filtered,
                            colData=metadata,
                            design=~condition)
# Create dataset object combining counts and metadata for DESeq2.
                            
ddc <-DESeq(ddc)
# Perform normalization and differential expression testing.

res <- results(ddc, contrast = c("condition", "risk", "prot"))
# Extract differential expression results comparing risk vs protective variants.

summary(res)
# Review summary statistics of differential expression results.

plotMA(res, main = "MA plot: Risk vs Protective", ylim = c(-5, 5))
# Visualize log fold changes vs mean expression highlighting significant genes.

legend('topright', 
       legend=c('Significant genes', 'Non-significant genes'),
       col= c('blue', 'grey'),
       pch=20,
       pt.cex=2)
# Add legend to distinguish significant from non-significant genes in the plot.

#_ Identity of upregulated genes ---------- 

upregulated_genes<-subset(res, padj<0.05 & log2FoldChange >0)
up_gene_names<-rownames(upregulated_genes)
head(up_gene_names)
# Identify genes significantly upregulated in risk variant group.


#_ Identity of downregulated genes ---------- 

downregulated_genes<-subset(res, padj<0.05 & log2FoldChange <0)
down_gene_names<-rownames(downregulated_genes)
head(down_gene_names)
# Identify genes significantly downregulated in risk variant group.

