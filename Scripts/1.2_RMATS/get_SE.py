import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os
from pathlib import Path
import sys


# STEP 0: Create directories to store RMATS files
output_dir = str(Path(os.getcwd())) + "/0_Files/RMATS/"
Path(output_dir).mkdir(parents=True, exist_ok=True)


# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
file = sys.argv[1] + 'RMATS/SE.MATS.JC.txt'
rmats = pd.read_csv(file, delimiter='\t')
col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
rmats = rmats[col_list]

print('Processing RMATS output: Skipped Exons \n')
print('# genes reported:                ', len(set(rmats.geneSymbol.values.tolist()))) # log

# use | dPSI | and only true values
rmats['IncLevelDifference'] = rmats['IncLevelDifference'].abs()
rmats = rmats[rmats['FDR'] <=0.05]

print('FDR-adj pvalue <= 0.05:          ', len(set(rmats.geneSymbol.values.tolist()))) # log

if len(rmats) == 0:
    print(' No skipped exons to process \n')
    sys.exit(0)

# FILTER 1: Get AS ( |dPSI| > 0.2, FDR < 0.05)
rmats_AS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() >= 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]

if len(rmats_AS) == 0:
    print(' No skipped exons to process \n')
    sys.exit(0)

# FILTER 2: If skipped exon is reported many times,  pick single dPSI score (can happen if down/upstream exons vary)

## get the largest dPSI value for AS exons (most differentially used score)
df = rmats_AS.copy()
# Create 'dPSI' and 'dPSI' columns
df['dPSI'] = df.groupby('exonStart_0base')['IncLevelDifference'].transform(lambda x: ','.join(x.astype(str)))
df['dPSI'] = df['dPSI'].str.split(',').apply(lambda x: max(map(float, x)) if x[0] else None)

# Keep only rows where 'IncLevelDifference' is equal to 'dPSI'
df = df[df['IncLevelDifference'] == df['dPSI']]

# FILTER 3: Drop duplicate exon entries
df = df.drop_duplicates(subset=["GeneID", "strand", "exonStart_0base", "exonEnd"], keep='first')

# Assign the modified DataFrame back to the original variable
rmats_AS = df

# FILTER 4: Keep coords of single version of exon if A3SS/A5SS events exist (to prevent 2+ flanks per exon)

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

rmats_AS = rmats_AS.groupby('GeneID').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
rmats_AS.reset_index(drop=True, inplace=True)

print('| IncLevelDifference | > 0.2:    ', len(set(rmats_AS.geneSymbol.values.tolist()))) # log

## STEP 2: Prepare bedtools input

# temp output fiilee
df = rmats_AS.copy()
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
