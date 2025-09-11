#!/bin/bash

## Check paths. Hardcoded ##

# STEP 1: Get count matrices
featureCounts -O -s 2 -p --countReadPairs -B -t exon -g gene_name -T 8 -F GTF -a ~/data/ref_genome/gencode.v24.PRI.gtf -o counts.txt ~/data/bam_dir/*.bam ~/data/bam_dir/validation/*.bam

# -s 2 reverse
# -O overlapping reads
# -p specify reads are paired ended
# --countReadPairs count reads as pairs
# -B consider read-pairs of good quality
# -t count reads mapping to exons
# -g summarize counts at gene level

wait $!
# STEP 2: Convert geneID to geneSymbol in count matrices
Rscript gene_id_tosymbol.R

wait $!
# STEP 3: Run edgeR - get TPM values
Rscript get_TPM.R

wait $!
mkdir  -p 0_Files/Post-processing/Analyses/expression/counts
mv rpkm_values_rbps.tsv 0_Files/Post-processing/Analyses/expression/counts/tpm_values_rbps.tsv