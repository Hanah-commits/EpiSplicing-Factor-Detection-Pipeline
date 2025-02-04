#!/bin/bash

# Check if the user provided enough arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file="$1"  # Input file containing the data

# Process the input file
while read -r line; do
    # Split the line into fields
    IFS=$'\t' read -r chr col1 col2 feature score strand gene HM col3 col4 col5 <<< "$line"

    # Split both col4 and col5 by comma
    IFS=',' read -r -a col4 <<< "$col4"
    IFS=',' read -r -a col5 <<< "$col5"

    # Loop through the values in col4 and col5 to create combinations
    for i in "${!col4[@]}"; do
        for j in "${col5[@]}"; do
            # Print the relevant fields with the current values from col5 and last_col
            echo -e "$chr\t${col4[i]}\t$j\t$gene\t$HM\t$strand"
        done
    done

done < "$input_file"
