suppressPackageStartupMessages( library( "DEXSeq" ) )

## read output files from previous steps into R
inDir <- getwd()  
countFiles <- list.files(inDir, pattern="\\.txt$", full.names=TRUE)

flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)

# make a dataframe of metadata
sampleTable <- data.frame(
  row.names = sub("\\.txt$", "", basename(countFiles)),
  condition = sub("_.*", "", basename(countFiles)),
  libType = rep("paired-end", length(countFiles))
)

#  construct an DEXSeqDataSet object
dxd = DEXSeqDataSetFromHTSeq(
countFiles,
sampleData=sampleTable,
design= ~ sample + exon + condition:exon,
flattenedfile=flattenedFile )


# Run DEU analysis
dxr1 = DEXSeq(dxd)


dxr2 = as.data.frame(dxr1)
# Identify non-list and non-matrix columns
non_list_matrix_cols <- sapply(dxr2, function(x) !is.list(x) && !is.matrix(x))
# Select columns of type "character" or "double"
selected_columns <- dxr2[, non_list_matrix_cols]
# Write to TSV
write.table(selected_columns, file = "DEXSEQ_H1_neuonral.tsv", sep = "\t", row.names = FALSE)