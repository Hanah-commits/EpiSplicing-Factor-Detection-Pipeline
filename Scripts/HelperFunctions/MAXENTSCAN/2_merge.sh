#!/bin/bash

# Check if the user provided enough arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <number_of_recursions>"
    exit 1
fi

input_file="$1"               # Input file containing the data
num_recursions="$2"           # Number of recursions to determine how many columns to merge

# Calculate the number of columns to merge (last n*2 columns)
cols_to_merge=$((num_recursions * 2))

# Process the input file
while read -r line; do
    # Remove the last n*2 columns and keep the rest of the line
    remaining_columns=$(echo "$line" | awk -v cols="$cols_to_merge" -F'\t' '{for(i=1; i<=NF-cols; i++) printf "%s\t", $i; print ""}')
    
    # Extract the last n*2 columns (split using tabs, not spaces)
    last_columns=$(echo "$line" | awk -v cols="$cols_to_merge" -F'\t' '{for(i=NF-cols+1; i<=NF; i++) printf "%s\t", $i; print ""}')
    
    # Split the last_columns into an array
    IFS=$'\t' read -ra columns <<< "$last_columns"
    
    # Arrays to store unique first and second column values
    first_values=()
    second_values=()

    # Loop through the array, process the values two by two
    for (( i=0; i<${#columns[@]}; i+=2 )); do
        first_col="${columns[i]}"
        second_col="${columns[i+1]}"

        # Only add non-dash unique values to respective arrays
        if [[ "$first_col" != "-" && ! " ${first_values[@]} " =~ " ${first_col} " ]]; then
            first_values+=("$first_col")
        fi
        if [[ "$second_col" != "-" && ! " ${second_values[@]} " =~ " ${second_col} " ]]; then
            second_values+=("$second_col")
        fi
    done
    
    # Join the unique values using commas for both first and second column arrays
    first_output=$(IFS=, ; echo "${first_values[*]}")
    second_output=$(IFS=, ; echo "${second_values[*]}")

    # Print the remaining columns followed by the merged output
    echo -e "${remaining_columns}\t${first_output}\t${second_output}"
    
done < "$input_file"
