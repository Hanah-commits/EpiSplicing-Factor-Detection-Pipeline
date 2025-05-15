# Load necessary library
library(edgeR)

# Read in the featureCounts output
file = commandArgs(trailingOnly = TRUE)[1]

featureCounts_data <- read.delim(file, comment.char = "#")

# Extract the counts matrix and gene lengths
counts <- featureCounts_data[, -c(1, 2, 3, 4, 5, 6, ncol(featureCounts_data))]
gene_lengths <- featureCounts_data$Length
gene_names <- featureCounts_data$gene_name

# Create a DGEList object
dge <- DGEList(counts=counts)

# Calculate RPKM
rpkm_values <- rpkm(dge, gene.length=gene_lengths)

# Convert RPKM values to a data frame and add gene names
rpkm_df <- data.frame(Geneid=gene_names, rpkm_values)

# Perform the log2(rpkm + 1) transformation
rpkm_df[,-1] <- log2(rpkm_df[,-1] + 1)

# Save the transformed RPKM values to a TSV file
write.table(rpkm_df, file="rpkm_values_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Read the text file with RBP names to include
gene_list <- readLines("proteins.txt")

# Filter the RPKM values data frame to include only genes in the gene list
filtered_rpkm_df <- rpkm_df[rpkm_df$Geneid %in% gene_list, ]

# Save the filtered RPKM values to a TSV file
write.table(filtered_rpkm_df, file="rpkm_values_rbps.tsv", sep="\t", quote=FALSE, row.names=FALSE)