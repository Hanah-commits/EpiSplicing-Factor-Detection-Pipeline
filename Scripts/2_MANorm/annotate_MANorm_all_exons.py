import os
import sys
import json
import pandas as pd
from pathlib import Path

# STEP 0: Create directories to store MANorm files
output_dir = str(Path(os.getcwd())) + "/0_Files/MANorm/"
Path(output_dir).mkdir(parents=True, exist_ok=True)


with open('paths.json') as f:
    d = json.load(f)

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
hms = d["Histone modifications"]

prefix = sys.argv[1]+ 'MANorm/'

for hm in hms:
    
    input = prefix + hm + '_' + tissue1 + '_peak_vs_' + hm + '_' + tissue2 +  '_peak_all_MAvalues.xls'
    output1 = output_dir+ hm + '_all_exons.bed'
    output2 = output_dir+ hm + '_exons.bed'

    ## STEP 1: annotate all exons
    os.system('bedtools intersect -loj -a 0_Files/all_exons.bed -b ' + input + ' | sort | uniq > ' + output1)
    ## STEP 12: annotate non-tss-overlap-exons
    os.system('bedtools intersect -loj -a 0_Files/exon_coords.bed -b ' + input + ' | sort | uniq > ' + output2)

    ## STEP 3: mark exons overlapping with TSS

    # read output into dataframe
    exons = pd.read_csv(output1, delimiter='\t', header=None)
    exons.drop([12, 14, 15, 16], axis=1, inplace=True)
    exons.columns = ['chr', "exon_start", "exon_end", "feature", "score", "strand", "geneSymbol", 'chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']

    non_tss_exons = pd.read_csv( output2, delimiter='\t', header=None)
    non_tss_exons.drop([12, 14, 15, 16], axis=1, inplace=True)
    non_tss_exons.columns = ['chr', "exon_start", "exon_end", "feature", "score", "strand", "geneSymbol", 'chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']

    # # Create a new column 'TSS_exon' in exons with initial value 'True'
    exons['TSS_exon'] = True

    # Use merge to check for matches and update 'TSS_exon' accordingly
    merged = pd.merge(exons, non_tss_exons, how='inner', suffixes=('', '_non_tss'))

    # Update 'TSS_exon' to 'False' for rows that have matches
    exons.loc[exons.index.isin(merged.index), 'TSS_exon'] = False

    ## STEP 4: get length of overlap
    exons['overlap_bp'] = exons.apply(lambda row: max(0, min(row['exon_end'], row['peak_end']) - max(row['exon_start'], row['peak_start'])) 
                            if row['peak_start'] != 0 and row['peak_end'] != 0 else 0, axis=1)
    
    #'overlap_bp_norm' based on the normalized overlap
    exons['overlap_bp_norm'] = exons.apply(lambda row: row['overlap_bp'] / (row['exon_end'] - row['exon_start']) if row['peak_start'] != 0 and row['peak_end'] != 0 else 0, axis=1)

    # assign 0 to flanks that have no peak coords
    exons.replace([-1, '.'], [0, 0], inplace=True)

    ## STEP 5: find if same histone peak is annotated to multiple exons
    duplicates = exons[exons.duplicated(subset=['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value'], keep=False) &
                    (exons['peak_start'] != 0) & (exons['peak_end'] != 0) & (exons['summit'] != 0) & (exons['M_value'] != 0)]

    ## STEP 5.1: Assign peaks to exons where TSS_exon is true
    mask = (duplicates.groupby(['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value'])['TSS_exon'].transform('any'))
    duplicates.loc[~mask, 'TSS_exon'] = False
    # Mark the rest with 0
    duplicates.loc[~mask, ['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']] = 0
    exons.update(duplicates)

    ## STEP 5.2: Assign peaks to exons with at least 50% exon overlap
    exons.loc[exons['overlap_bp_norm'] < 0.5, ['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']] = 0

    # print(len(exons[exons.chr_2 !=0][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates()))



