import os
import sys
import json
import pandas as pd
from pathlib import Path
import numpy as np

def adjust_pvalue(df, col):

    # get indices of null values
    na_idx = df[df[col].isnull()].index.tolist()

    # adjust non-null p values
    pvals = df[col].values.tolist()
    pvals = [x for x in pvals if str(x) != 'nan']
    adj_pval = p_adjust_bh(pvals).tolist()

    # insert null at original indices
    for idx in na_idx:
        adj_pval.insert(idx, None)

    # adjusted p values as new df
    df['adj_pval'] = adj_pval
    return df


def p_adjust_bh(p):
    ## multiple hypothesis testing
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]



# STEP 0: Create directories to store MANorm files
output_dir = str(Path(os.getcwd())) + "/0_Files/MANorm/"
Path(output_dir).mkdir(parents=True, exist_ok=True)


with open('paths.json') as f:
    d = json.load(f)

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
hms = d["Histone modifications"]
ref = d['Reference genome']
fasta = d['Reference fasta']
ref_genome= fasta+".fai"

prefix = sys.argv[1]+ 'MANorm/'

for hm in hms:
    print(hm)

    input = prefix + hm + '_' + tissue1 + '_peak_vs_' + hm + '_' + tissue2 +  '_peak_all_MAvalues.xls'
    output1 = output_dir+ hm + '_all_exons.bed'
    output2 = output_dir+ hm + '_exons.bed'

    # STEP 1: annotate all exons
    os.system('bedtools intersect -loj -a 0_Files/all_exons.bed -b ' + input + ' | sort | uniq > ' + output1)
    ## STEP 2: annotate non-tss-overlap-exons
    os.system('bedtools intersect -loj -a 0_Files/exon_coords.bed -b ' + input + ' | sort | uniq > ' + output2)

    ## STEP 3.a: mark exons overlapping with TSS

    # read output into dataframe
    exons = pd.read_csv(output1, delimiter='\t', header=None)
    exons.drop([12, 14, 15, 16], axis=1, inplace=True)
    exons.columns = ['chr', "exon_start", "exon_end", "feature", "score", "strand", "geneSymbol", 'chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']

    non_tss_exons = pd.read_csv( output2, delimiter='\t', header=None)
    non_tss_exons.drop([12, 14, 15, 16], axis=1, inplace=True)
    non_tss_exons.columns = ['chr', "exon_start", "exon_end", "feature", "score", "strand", "geneSymbol", 'chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']

    # # Create a new column 'TSS_exon' in exons with initial value 'False'
    exons['TSS_exon'] = False

    exons.reset_index(drop=True, inplace=True)
    non_tss_exons.reset_index(drop=True, inplace=True)

    # Use merge to check for matches and update 'TSS_exon' accordingly
    merged = exons.merge(non_tss_exons, on=['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'geneSymbol'], how='left', suffixes=('', '_non_tss'))
    # Drop duplicate rows while keeping the first occurrence
    merged.drop_duplicates(subset=['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'geneSymbol', 'chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value', 'TSS_exon'], keep='first', inplace=True)

    # Realign indices of both DataFrames
    exons.reset_index(drop=True, inplace=True)
    merged.reset_index(drop=True, inplace=True)

    # Update TSS_exon based on the merge (exons not in non_tss are TSS exons)
    exons.loc[merged.index, 'TSS_exon'] = merged['M_value_non_tss'].isna()

    ## STEP 3.b: keep M-values of peaks with adj P-value <= 0.05
    # assign 0 to exons that have no peak coords
    exons.replace([-1, '.'], [0, 0], inplace=True)
    exons['p_value'] = pd.to_numeric(exons['p_value'])
    exons = adjust_pvalue(exons, col='p_value')
    exons.loc[pd.to_numeric(exons['adj_pval']) > 0.05, ['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value']] = 0 

    print('Post FDR-filtering')
    print('All annotated exons:       ', len(exons[exons.chr_2 !=0][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates()))
    print('Annotated non-TSS exons:   ', len(exons[(exons.chr_2 !=0) & (exons.TSS_exon == False)][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates())) # final num non-TSS exons
    print('Annotated TSS exons:       ', len(exons[(exons.chr_2 !=0) & (exons.TSS_exon == True)][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates())) # final num non-TSS exons
    print('All annotated Genes:       ', len(set(exons[(exons.chr_2 !=0)]['geneSymbol'].values.tolist())), '\n') # genes with peak-annotated exons  

    ## STEP 4: get length of overlap
    exons['overlap_bp'] = exons.apply(lambda row: max(0, min(row['exon_end'], row['peak_end']) - max(row['exon_start'], row['peak_start'])) 
                            if row['peak_start'] != 0 and row['peak_end'] != 0 else 0, axis=1)
    
    #'overlap_bp_norm' based on the normalized overlap (by exon length)
    exons['overlap_bp_norm'] = exons.apply(lambda row: row['overlap_bp'] / (row['exon_end'] - row['exon_start']) if row['peak_start'] != 0 and row['peak_end'] != 0 else 0, axis=1)


    ## STEP 5: find if same histone peak is annotated to multiple exons
    duplicates = exons[exons.duplicated(subset=['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value'], keep=False) &
                    (exons['peak_start'] != 0) & (exons['peak_end'] != 0) & (exons['summit'] != 0) & (exons['M_value'] != 0)]
    
    ## STEP 5.1: Assign peaks to exons where TSS_exon is true (TSS-peak-leak)
    # Identifying duplicated rows
    duplicated_rows = exons.duplicated(subset=['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value'], keep=False)

    # Check if any of the duplicated rows have TSS_exon set to True
    mask = (duplicated_rows & duplicates['TSS_exon'])

    # Set values to zero where TSS_exon is False and any of the duplicated rows have TSS_exon set to True
    duplicates.loc[(duplicated_rows) & (~mask) & (~duplicates['TSS_exon']), ['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']] = 0

    exons.update(duplicates)

    print('Post TSS Peak-leak Filtering')
    print('All annotated exons:       ', len(exons[exons.chr_2 !=0][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates()))
    print('Annotated non-TSS exons:   ', len(exons[(exons.chr_2 !=0) & (exons.TSS_exon == False)][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates())) # final num non-TSS exons
    print('Annotated TSS exons:       ', len(exons[(exons.chr_2 !=0) & (exons.TSS_exon == True)][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates())) # final num non-TSS exons
    print('All annotated Genes:       ', len(set(exons[(exons.chr_2 !=0)]['geneSymbol'].values.tolist())), '\n') # genes with peak-annotated exons  

    ## STEP 5.2: Assign peaks to exons with at least 50% exon overlap
    exons.loc[exons['overlap_bp_norm'] < 0.5, ['chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'p_value']] = 0

    ## save non-TSS exons with peaks
    exons[(exons.TSS_exon == False)].to_csv(output2, sep='\t', header=False, index=False)

    print('Post Neighboring Exon Peak-leak Filtering')
    print('All annotated exons:       ', len(exons[exons.chr_2 !=0][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates()))
    print('Annotated non-TSS exons:   ', len(exons[(exons.chr_2 !=0) & (exons.TSS_exon == False)][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates())) # final num non-TSS exons
    print('Annotated TSS exons:       ', len(exons[(exons.chr_2 !=0) & (exons.TSS_exon == True)][['chr', 'exon_start', 'exon_end', 'strand']].drop_duplicates())) # final num non-TSS exons
    print('All annotated Genes:       ', len(set(exons[(exons.chr_2 !=0)]['geneSymbol'].values.tolist()))) # genes with peak-annotated exons  
    print('Candidate annotated Genes: ', len(set(exons[(exons.chr_2 !=0) & (exons.TSS_exon == False)]['geneSymbol'].values.tolist())), '\n')# final num of genes with annotated non-TSS exons  


    ## STEP 6: Assign peaks to flanks of exons

    ## STEP 6.1: make flanks of filtered exons
    
    # get non-TSS exons + peaks
    exons = exons[(exons.TSS_exon == False) & (exons.chr_2 != 0)]
    # save exons
    exons['flank_feature'] = 'flank'
    exons[['chr', 'exon_start', 'exon_end', "flank_feature", "score", "strand", "geneSymbol"]].drop_duplicates().to_csv(output_dir+f'{hm}_filtered_exons.bed', sep='\t', index=False, header=False)
    # save peaks
    exons['peak_feature'] = hm + '_peak'
    exons[['chr_2', 'peak_start', 'peak_end', "peak_feature", "M_value"]].drop_duplicates().to_csv(output_dir+f'{hm}_filtered_peaks.bed', sep='\t', index=False, header=False)

    # exon boundary external flanks
    os.system(f"bedtools flank -i {output_dir}/{hm}_filtered_exons.bed -g " + ref_genome + " -b 200 > 0_Files/flanks.bed" )

    # separate start,stop flank coords
    os.system("sed -n 'n;p' 0_Files/flanks.bed > 0_Files/stop.bed")
    os.system("sed -n 'p;n' 0_Files/flanks.bed > 0_Files/start.bed")

    # exon boundary internal flanks
    os.system("bedtools slop -i 0_Files/start.bed -g " + ref_genome + " -l 0 -r 200 > 0_Files/start_flanks.bed")
    os.system("bedtools slop -i 0_Files/stop.bed -g " + ref_genome +" -l 200 -r 0 > 0_Files/stop_flanks.bed")

    # combine start,stop flank coords
    os.system(f"paste -d'\n' 0_Files/start_flanks.bed 0_Files/stop_flanks.bed | sort -k1,1 -k2,2n > {output_dir}/{hm}_unannotated_flanks.bed")

    # remove intermediate files
    os.system("rm 0_Files/start*.bed")
    os.system("rm  0_Files/stop*.bed")
    os.system("rm 0_Files/flanks.bed")

    ## STEP 6.2: annotate flanks of filtered exons with peaks
    os.system(f'bedtools intersect -loj -a {output_dir}/{hm}_unannotated_flanks.bed -b {output_dir}/{hm}_filtered_peaks.bed | sort | uniq > {output_dir}/{hm}_annotated_flanks.bed') 

    ## STEP 6.3: Deal with peak-leak in annotated flanks

    flanks = pd.read_csv(f'{output_dir}/{hm}_annotated_flanks.bed', delimiter='\t', header=None)
    flanks.columns = ['chr', "flank_start", "flank_end", "feature", "score", "strand", "geneSymbol", 'chr_2', 'peak_start', 'peak_end', 'peak_feature', 'M_value']

    ## STEP 6.3.1: get length of overlap
    flanks['overlap_bp'] = flanks.apply(lambda row: max(0, min(row['flank_end'], row['peak_end']) - max(row['flank_start'], row['peak_start'])) 
                            if row['peak_start'] != 0 and row['peak_end'] != 0 else 0, axis=1)

    ## STEP 6.3.2: find if same histone peak is annotated to multiple exons
    duplicates = flanks[flanks.duplicated(subset=['chr', 'flank_start', 'flank_end', 'geneSymbol'], keep=False) &
                    (flanks['peak_start'] != -1) & (flanks['peak_end'] != -1) & (flanks['M_value'] != -1)]
    
    # Find the index of the row with the maximum 'overlap_bp' within each group
    ## if single flank has 1+ peak annotations, choose using max overlap
    max_overlap_idx = flanks.groupby(['chr', 'flank_start', 'flank_end', 'geneSymbol'])['overlap_bp'].idxmax()

    # Select the corresponding rows based on the index
    flanks = flanks.loc[max_overlap_idx]

    ## STEP 6.4: Deal with multiple peaks annotated to same flank (precaution)
    flanks.loc[flanks['chr_2'] == '.', ['peak_start', 'peak_end', 'peak_feature']] = '.'
    flanks.loc[flanks['chr_2'] == '.', ['M_value']] = 0
    flanks['M_value_abs'] = pd.to_numeric(flanks['M_value']).abs()

    # get all peaks that belong to each flank
    flank_peaks_group = flanks.groupby(['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'])['M_value_abs'] \
        .apply(lambda val: ','.join(str(v) for v in val)).reset_index()

    # # FILTER 1: if flank has 1+ peaks, keep peak with highest abs M-value
    flank_peaks_group['M_value_abs'] = flank_peaks_group['M_value_abs'].str.split(',')  # string -> list of strings
    flank_peaks_group['max_' + hm] = flank_peaks_group['M_value_abs'].apply(lambda x: max(map(float, x)))  # max MValue

    # # get the corresponding peak for each flank's max M-value
    flanks = pd.merge(flank_peaks_group[['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol', 'max_' + hm]], flanks, on=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'],
                                how='inner')
    # # flank_peaks_group = pd.merge(flank_peaks_group[['flanks', 'max_'+hm, '#peaks_'+hm]], peaks, on=['flanks'], how='inner')
    flanks = flanks[flanks['M_value_abs'] == flanks['max_' + hm]]

    # FILTER 2: If flank has 1+ peaks with same max |Mvalue|, keep one
    flanks.drop_duplicates(subset=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'], keep='first', inplace=True)
    flanks = flanks[['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol', 'peak_start', 'peak_end', 'peak_feature', f'max_{hm}']]

    ## save filtered flanks
    flanks.to_csv(f'{output_dir}/{hm}_annotated_flanks.bed', sep='\t', header=False, index=False)