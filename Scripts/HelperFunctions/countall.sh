#!/bin/bash

# Iterate over each file passed as an argument
for file in "$@"; do
    # Extract the filename without the extension
    filename=$(basename -- "$file")
    filename_noext="${filename%.*}"

    # Run the Python command with the current file
    python /home/hanah/miniconda3/envs/majiq/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -s reverse ./gencode.v24.DEXSEQ.gff "$file" "${filename_noext}.txt"

    # Remove HTSEQ count info (last five lines)
    head -n -5 "${filename_noext}.txt" > temp.txt && mv temp.txt "${filename_noext}.txt"
    
    # Remove quotes in count file
    cp "${filename_noext}.txt" temp.txt  
    sed 's/\"//g' temp.txt > temp1.txt && mv temp1.txt temp.txt 
    mv temp.txt "${filename_noext}.txt"

done