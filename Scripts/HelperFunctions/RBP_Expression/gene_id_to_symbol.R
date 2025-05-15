library(biomaRt)


file = commandArgs(trailingOnly = TRUE)[1]
# Read the TSV file into a DataFrame
data <- read.table(file, skip=1, header=TRUE, sep='\t')

# Extract Ensembl IDs
data$ensembl_ids = sub("\\..*$", "", data$Geneid) # ENSG00000183878.15 -> ENSG00000183878

ensembl_ids <- unique(data$ensembl_ids)
ensembl <- useMart("ensembl","hsapiens_gene_ensembl")

# Get gene symbols using BioMart
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = ensembl)

# Identify duplicated gene names
duplicated_gene_names <- duplicated(gene_symbols$external_gene_name) | duplicated(gene_symbols$external_gene_name, fromLast = TRUE)

# Replace duplicate gene names with NA
gene_symbols$external_gene_name[duplicated_gene_names] <- NA

# Merge gene symbols back into the DataFrame
merged_data <- merge(data, gene_symbols, by.x = "ensembl_ids", by.y = "ensembl_gene_id", all.x = TRUE)

# Identify rows with missing or whitespace external_gene_name
missing_or_whitespace <- is.na(merged_data$external_gene_name) | merged_data$external_gene_name == " "

if (any(missing_or_whitespace)) {
    # Replace missing or whitespace values with corresponding geneSymbol values
    merged_data$external_gene_name[missing_or_whitespace] <- merged_data$Geneid[missing_or_whitespace]
}

# Remove ensembl_ids and geneSymbol columns
merged_data <- merged_data[, !names(merged_data) %in% c("ensembl_ids")]

# Rename external_gene_name to gene_name
names(merged_data)[names(merged_data) == "external_gene_name"] <- "gene_name"

colnames(merged_data) <- gsub('X.home.user.data.bam_dir.validation.', "", colnames(merged_data), ignore.case = TRUE)


# Save the dataframe as a TSV file
write.table(merged_data, file = "./counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)