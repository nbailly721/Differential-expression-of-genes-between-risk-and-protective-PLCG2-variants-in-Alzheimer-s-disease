# Differential expression of genes between risk and protective PLCG2 variants in Alzheimer's disease

#_ Set up of environment ---------- 
install.packages('R.utils')
install.packages('tidyverse')

# Load libraries
library(R.utils)
library(tidyverse)


#_ Load data file ---------- 
data_raw <- read_csv('/Users/nicolasbailly/Documents/side_projects/Alzheimers/GSE308813_20240829_LB_rna_counts_only.csv')

#_ Data Exploration and modification---------- 

# Check the column names
print(names(data_raw))

# Rename columns with meaningful names
new_names <- c('gene','WT1','WT2','WT3','WT4',
               'KO1','KO2','KO3','KO4',
               'risk1','risk2','risk3','risk4',
               'prot1','prot2','prot3','prot4')

colnames(data_raw) <- new_names

# Extract gene names before filtering
gene_names <- data_raw$gene

# Extract numeric data (excluding gene column)
data_subset_numeric <- data_raw %>% select(-gene)

# Calculate row sums for filtering
row_sums_vector <- rowSums(data_subset_numeric)

# Create logical filter for genes with total counts > 100
keep_rows <- row_sums_vector > 100

# Filter numeric data and gene names together
data_filtered <- data_subset_numeric[keep_rows, ]
filtered_gene_names <- gene_names[keep_rows]

# Make gene names unique to avoid duplicate rowname error
filtered_gene_names_unique <- make.unique(filtered_gene_names)

# Convert to data frame and assign unique gene names as row names
data_filtered <- as.data.frame(data_filtered)
rownames(data_filtered) <- filtered_gene_names_unique

# Create metadata matching samples (columns)
metadata <- data.frame(
  sample = colnames(data_filtered),
  condition = c(rep('WT',4), rep('KO',4), rep('risk',4), rep('prot',4))
)

# Print metadata and check filtered data
print(metadata) #data frame
print(head(data_filtered))  #matrix

# Check first few gene names as rownames
print(head(rownames(data_filtered)))

#_ Statistical analysis using DESeq2 ---------- 

#Load necessary packages
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

library(DESeq2)

#Ensure that the conditions in the metadata are factors
metadata$condition<-factor (metadata$condition, levels=c('prot','risk')

#Create DESeq2 dataset object
ddc <- DESeqDataSetFromMatrix(countData=data_filtered,
                            colData=metadata,
                            design=~condition)

#Normalizes data and performs statistcal tests
ddc <-DESeq(ddc)

#Compares expression results of both groups
res <- results(ddc, contrast = c("condition", "risk", "prot"))

#Summarize results
summary(res)

#Visualization of results-MA plot
plotMA(res, main = "MA plot: Risk vs Protective", ylim = c(-5, 5))

#_ Identity of upregulated genes ---------- 

upregulated_genes<-subset(res, padj<0.05 & log2FoldChange >0)
up_gene_names<-rownames(upregulated_genes)
head(up_gene_names)

#_ Identity of downregulated genes ---------- 

downregulated_genes<-subset(res, padj<0.05 & log2FoldChange <0)
down_gene_names<-rownames(downregulated_genes)
head(down_gene_names)
