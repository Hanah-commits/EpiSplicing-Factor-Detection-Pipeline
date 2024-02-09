import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os
from pathlib import Path


# STEP 0: Create directories to store RMATS files
output_dir = str(Path(os.getcwd())) + "/0_Files/RMATS/"
Path(output_dir).mkdir(parents=True, exist_ok=True)


# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
file = '/Users/hanah/EpiSplicing_RMATS/Output/neuro-H1/RMATS/SE.MATS.JC.txt'
rmats = pd.read_csv(file, delimiter='\t')
col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
rmats = rmats[col_list]

# use | dPSI | and only true values
rmats['IncLevelDifference'] = rmats['IncLevelDifference'].abs()
rmats = rmats[rmats['FDR'] <=0.05]


# FILTER 0: Filter out exons reported by RMATS that don't belong to transcripts with TSL 1-3

# write df into bed file
rmats['feature'] = "flank"
rmats['score'] = "."
rmats[['chr', 'exonStart_0base', 'exonEnd', 'feature', 'score', 'strand']].to_csv('0_Files/RMATS/rmats_query.bed', index=False, sep='\t', header=False )

# run bedtools
os.system('bedtools intersect -a 0_Files/RMATS/rmats_query.bed -b 0_Files/exon_coords.bed -wa | sort | uniq > 0_Files/RMATS/rmats_result.bed')

rmats_filtered = pd.read_csv('0_Files/RMATS/rmats_result.bed', delimiter='\t', header=None) # 6489
rmats_filtered.columns = ['chr', 'exonStart_0base', 'exonEnd', 'feature', 'score', 'strand']
rmats = pd.merge(rmats, rmats_filtered, on=['chr', 'exonStart_0base', 'exonEnd', 'feature', 'score', 'strand'], how='inner')
col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
rmats = rmats[col_list]

# housekeeping
os.system('rm 0_Files/RMATS/rmats_*.bed')

print(rmats[rmats.GeneID == 'ENSG00000237441.9'])


# FILTER 1: Get AS ( |dPSI| > 0.2, FDR < 0.05) and CS exons ( |dPSI| < 0.2, FDR < 0.05)
rmats_AS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() >= 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]
rmats_CS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() < 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]

# FILTER 2: If skipped exon is reported many times,  pick single dPSI score (can happen if down/upstream exons vary)

## get the largest dPSI value for AS exons (most differentially used score)
for i, df in enumerate([rmats_AS, rmats_CS]):
    # Create 'dPSI' and 'dPSI' columns
    df['dPSI'] = df.groupby('exonStart_0base')['IncLevelDifference'].transform(lambda x: ','.join(x.astype(str)))
    df['dPSI'] = df['dPSI'].str.split(',').apply(lambda x: max(map(float, x)) if x[0] else None)

    # Keep only rows where 'IncLevelDifference' is equal to 'dPSI'
    df = df[df['IncLevelDifference'] == df['dPSI']]

    # FILTER 3: Drop duplicate exon entries
    df = df.drop_duplicates(subset=["GeneID", "strand", "exonStart_0base", "exonEnd"], keep='first')

    # Assign the modified DataFrame back to the original variable
    if i == 0:
        rmats_AS = df
    else:
        rmats_CS = df

# FILTER 4: Get only the exons from rmats_CS that are unique to it (not in rmats_AS)
merged_df = pd.merge(rmats_CS, rmats_AS[["GeneID", "strand", "exonStart_0base", "exonEnd"]], on=["GeneID", "strand", "exonStart_0base", "exonEnd"], how='left', indicator=True)
rmats_CS = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])

# FILTER 5: Keep coords of single version of exon if A3SS/A5SS events exist (to prevent 2+ flanks per exon)

# #                   GeneID geneSymbol    chr strand  IncLevelDifference       FDR  exonStart_0base   exonEnd   dPSI
# # 4634  ENSG00000126456.15       IRF3  chr19      -               0.449  0.000202         49664442  49664673  0.449
# # 4636  ENSG00000126456.15       IRF3  chr19      -               0.584  0.000012         49664552  49664673  0.584
# # 4640  ENSG00000126456.15       IRF3  chr19      -               0.363  0.015466         49664586  49664673  0.363

def A3SS_A5SS_filter(group, subset_column):
    group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
    group.drop_duplicates(subset=['exonStart_0base'], keep='first', inplace=True)
    group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
    group.drop_duplicates(subset=['exonEnd'], keep='first', inplace=True)
    return group

for i, df in enumerate([rmats_AS, rmats_CS]):
    df = df.groupby('GeneID').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
    df = df.reset_index(drop=True)

    # Assign the modified DataFrame back to the original variable
    if i == 0:
        rmats_AS = df
    else:
        rmats_CS = df

# FILTER 6: Drop genes that only have exons with DEU scores < 0.2 (no alternate exons)
rmats = pd.concat([rmats_AS, rmats_CS],axis=0,sort=False).reset_index()
rmats = rmats.groupby('GeneID').filter(lambda x: (x['dPSI'] > 0.2).any())


## STEP 2: Prepare bedtools input

# temp output fiilee
df = rmats.copy()
df['feature'] = "Exon"
df['score'] = "."
df['exonStart_0base'] = pd.to_numeric(df['exonStart_0base']) + 1
df[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "GeneID", "dPSI"]].to_csv(f'0_Files/RMATS/SE_exons.tsv', index=False, sep='\t', header=True)
df_temp = df.copy()
del(df_temp['exonStart_0base'])
del(df['exonEnd'])

df.rename(columns={'exonStart_0base': 'exon_coord0'}, inplace = True)
df_temp.rename(columns={'exonEnd': 'exon_coord0'}, inplace=True)

df = pd.concat([df_temp, df]).sort_index(kind='merge')

keep_cols = ['chr', 'exon_coord0', 'strand']
df_bed = df[keep_cols]
df_bed = df_bed.drop_duplicates()
# to fit bedtools input requirements
df_bed['exon_coord1'] = pd.to_numeric(df_bed['exon_coord0']) + 1
df_bed['feature'] = "flank"
df_bed['score'] = "."


df_bed = df_bed[['chr', "exon_coord0", "exon_coord1", "feature", "score", "strand"]]
df_bed.to_csv(f'0_Files/RMATS/SE.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
df.to_csv(f'0_Files/RMATS/SE_exons.csv', index=False, sep='\t', header=True)

    
    










def view_dpsi(type):

    if type == "sub":
        data1 = rmats_AS['IncLevelDifference'].abs().values.tolist()
        data2 = rmats_CS['IncLevelDifference'].abs().values.tolist()

        kde1 = gaussian_kde(data1)
        kde2 = gaussian_kde(data2)

        # Generate points on the x-axis for the KDE plots
        x1 = np.linspace(min(data1), max(data1), 1000)
        x2 = np.linspace(min(data2), max(data2), 1000)

        # Calculate the KDE values for both data lists
        kde_values1 = kde1(x1)
        kde_values2 = kde2(x2)

        # Create two subplots side by side
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        # Plot the KDE for data1 on the first subplot
        axs[0].plot(x1, kde_values1)
        axs[0].set_xlabel('| dPSI | values of DJU events (Alternative)')
        axs[0].set_ylabel('Density')


        # Plot the KDE for data2 on the second subplot
        axs[1].plot(x2, kde_values2)
        axs[1].set_xlabel('| dPSI | values of non-DJU events (constitutive)')
        axs[1].set_ylabel('Density')


        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()

    elif type == "main":
        data = rmats['IncLevelDifference'].abs().values.tolist()
        kde = gaussian_kde(data)
       

        # Generate points on the x-axis for the KDE plots
        x1 = np.linspace(min(data), max(data), 1000)

        # Calculate the KDE values for both data lists
        kde_values1 = kde(x1)

        # Create two subplots side by side
        fig, axs = plt.subplots(1, figsize=(12, 5))

        # Plot the KDE for data1 on the first subplot
        axs.plot(x1, kde_values1)
        axs.set_xlabel('| dPSI | values of SE events')
        axs.set_ylabel('Density')

        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()
