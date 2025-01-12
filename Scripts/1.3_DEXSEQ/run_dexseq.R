library('dplyr')
suppressPackageStartupMessages( library( "DEXSeq" ) )
library(jsonlite)

# Ensure that the correct number of command-line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 4) {
  stop("Usage: Rscript script.R directory_path tissue1 tissue2 processname")
}

proc = args[4]

json_data <- fromJSON("paths.json") 
json_data <- lapply(json_data, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

do.call("rbind", json_data)
json_data = json_data$proc
tissue1 <- json_data$tissue1
tissue2 <- json_data$tissue2
dexseq_dir <- json_data[["DEXSEQ directory"]]
output_dir <- json_data[["Output directory"]]

source(paste0(dexseq_dir,"python_scripts/load_SubreadOutput.R"))

# Retrieve command-line arguments
directory_path <- commandArgs(trailingOnly = TRUE)[1]
tissue1 <- commandArgs(trailingOnly = TRUE)[2]
tissue2 <- commandArgs(trailingOnly = TRUE)[3]

# List files based on the provided tissue types
files_tissue1 <- list.files(directory_path, pattern = paste0(tissue1, "_.*\\.bam$"), full.names = TRUE)
files_tissue2 <- list.files(directory_path, pattern = paste0(tissue2, "_.*\\.bam$"), full.names = TRUE)

all_files <- c(files_tissue1, files_tissue2)
flattenedFile = paste0(output_dir, 'DEXSEQ/DEXSEQ_reference.gtf')

# make a dataframe of metadata
sampleTable <- data.frame(
  row.names = sub("\\.bam$", "", basename(all_files)),
  condition = sub("_.*", "", basename(all_files)),
  libType = rep("paired-end", length(all_files))
)

#  construct an DEXSeqDataSet object
dxd <- DEXSeqDataSetFromFeatureCounts(paste0(output_dir,"DEXSEQ/",tissue1, "_", tissue2, "_count.out"), flattenedfile = flattenedFile, sampleData = sampleTable)

# Run DEU analysis
dxr1 = DEXSeq(dxd)

dxr2 = as.data.frame(dxr1)

# Identify non-list and non-matrix columns
non_list_matrix_cols <- sapply(dxr2, function(x) !is.list(x) && !is.matrix(x))

# Select columns of type "character" or "double"
selected_columns <- dxr2[, non_list_matrix_cols]

# Write to TSV
write.table(selected_columns, file = paste0(output_dir, "DEXSEQ/DEXSEQ_",tissue1,"_",tissue2,".tsv"), sep = "\t", row.names = FALSE)