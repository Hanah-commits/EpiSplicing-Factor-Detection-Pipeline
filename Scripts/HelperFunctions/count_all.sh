#!/bin/bash

# Get the paths to script and output dirs
json_file="paths.json"
field1="\"Output directory\""
field2="\"DEXSEQ directory\""
output_dir=$(jq -r ".[$field1]" "$json_file")
dexseq_dir=$(jq -r ".[$field2]" "$json_file")

# Get arguments
tissue1="${1}"
tissue2="${2}"
dir="${3}"

# Count using Subread
featureCounts -f -O -s 2 -p -T 40 \
-F GTF -a "$output_dir"DEXSEQ/DEXSEQ_reference.gtf -t exon \
-o "${output_dir}DEXSEQ/${tissue1}_${tissue2}_count.out" "${dir}/${tissue1}"_*.bam "${dir}/${tissue2}"_*.bam