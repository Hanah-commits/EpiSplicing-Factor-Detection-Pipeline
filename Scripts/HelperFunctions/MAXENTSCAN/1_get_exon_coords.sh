#!/bin/bash

# Define the path to your directories that contain rmats_exons_coords.bed files
parent_dir="/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output"
expression="*/0_Files/RMATS"
dir_list=($(find "$parent_dir" -type d -wholename "$expression"))
echo $dir_list


# Define the initial 'flanks.bed' file
flanks_bed="flanks.bed"

# Temporary file to store intermediate results
temp_bed="temp.bed"

# Loop through each directory and process the rmats_exons_coords.bed file
for dir in "${dir_list[@]}"; do
    rmats_file="$dir/rmats_exons_coords.bed"
    
    if [ -f "$rmats_file" ]; then
        echo "Processing $rmats_file with $flanks_bed"
        
        # Perform the AWK operations
        awk -F'\t' 'BEGIN { OFS="\t" }
        NR==FNR {
            key[$1":"$2] = $2":"$3;
            key[$1":"$3] = $2":"$3;
            next
        }
        {
            coord = $1":"$9;
            if (coord in key) {
                split(key[coord], vals, ":");
                print $0 "\t" vals[1] "\t" vals[2];
            } else {
                print $0 "\t-\t-";
            }
        }' "$rmats_file" "$flanks_bed" | awk -F'\t' '{
            sub(/\t$/, "");        # Remove the trailing tab at the end
            sub(/\t+/, "\t");     # Replace multiple tabs with a single tab
            print;                 # Print the modified line
        }' OFS='\t' > "$temp_bed" && mv "$temp_bed" "$flanks_bed"
    else
        echo "File $rmats_file does not exist. Skipping."
    fi
done

echo "All directories processed. Final output is in $flanks_bed"
