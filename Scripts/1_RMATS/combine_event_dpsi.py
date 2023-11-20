import pandas as pd
import os
import json

exon_types = ['AS']
event_types = ['A3SS', 'A5SS', 'SE', 'MXE']
flank_lens = [50, 100, 200]


# STEP 1: get combined dPSI scores of AS exons
A3SS = pd.read_csv(f'0_Files/Filtered_dPSI_A3SS_AS.csv', delimiter='\t', names=['chr', 'flanks', 'geneSymbol', 'A3SS_score'], skiprows=1)
A5SS = pd.read_csv(f'0_Files/Filtered_dPSI_A5SS_AS.csv', delimiter='\t',  names=['chr', 'flanks', 'geneSymbol', 'A5SS_score'], skiprows=1)
SE = pd.read_csv(f'0_Files/Filtered_dPSI_SE_AS.csv', delimiter='\t',  names=['chr', 'flanks', 'geneSymbol', 'SE_score'], skiprows=1)
MXE = pd.read_csv(f'0_Files/Filtered_dPSI_MXE_AS.csv', delimiter='\t',  names=['chr', 'flanks', 'geneSymbol', 'MXE_score'], skiprows=1)

merged_df_AS = pd.merge(A3SS, A5SS, on=['chr', 'flanks', 'geneSymbol'], how='outer')
merged_df_AS = pd.merge(merged_df_AS, SE, on=['chr', 'flanks', 'geneSymbol'], how='outer')
merged_df_AS = pd.merge(merged_df_AS, MXE, on=['chr', 'flanks', 'geneSymbol'], how='outer')
merged_df_AS.fillna(0, inplace=True)

# get | DPSI | scores
merged_df_AS[["A3SS_score", "A5SS_score", "SE_score", "MXE_score"]] = merged_df_AS[["A3SS_score", "A5SS_score", "SE_score", "MXE_score"]].abs()

# get single |DPSI| score per flank
merged_df_AS['dPSI'] = merged_df_AS[["A3SS_score", "A5SS_score", "SE_score", "MXE_score"]].max(axis=1)

# remove RNA genes (eg TMX2-CTNND1)
merged_df_AS[['start', 'stop']] = merged_df_AS['flanks'].str.split('-', n=1, expand=True)
duplicates_mask = merged_df_AS.duplicated(subset=['chr', 'start', 'stop'], keep=False)
merged_df_AS = merged_df_AS[~(duplicates_mask & (merged_df_AS['geneSymbol'].str.contains("-", na=False)))]

## output files
merged_df_AS[['flanks', 'geneSymbol', 'dPSI']].drop_duplicates().to_csv(f'0_Files/Filtered_dPSI_AS.csv', index=False, sep='\t')
merged_df_AS[['chr', 'start', 'stop']].drop_duplicates().to_csv(f'0_Files/filtered_flanks_AS.bed', index=False, sep='\t', header=False)


# STEP 2: get combined dPSI scores of CS exons

A3SS = pd.read_csv(f'0_Files/Filtered_dPSI_A3SS_CS.csv', delimiter='\t', names=['chr', 'flanks', 'geneSymbol', 'A3SS_score'], skiprows=1)
A5SS = pd.read_csv(f'0_Files/Filtered_dPSI_A5SS_CS.csv', delimiter='\t',  names=['chr', 'flanks', 'geneSymbol', 'A5SS_score'], skiprows=1)
SE = pd.read_csv(f'0_Files/Filtered_dPSI_SE_CS.csv', delimiter='\t',  names=['chr', 'flanks', 'geneSymbol', 'SE_score'], skiprows=1)
MXE = pd.read_csv(f'0_Files/Filtered_dPSI_MXE_CS.csv', delimiter='\t',  names=['chr', 'flanks', 'geneSymbol', 'MXE_score'], skiprows=1)

merged_df_CS = pd.merge(A3SS, A5SS, on=['chr', 'flanks', 'geneSymbol'], how='outer')
merged_df_CS = pd.merge(merged_df_CS, SE, on=['chr', 'flanks', 'geneSymbol'], how='outer')
merged_df_CS = pd.merge(merged_df_CS, MXE, on=['chr', 'flanks', 'geneSymbol'], how='outer')
merged_df_CS.fillna(0, inplace=True)

# get | DPSI | scores
merged_df_CS[["A3SS_score", "A5SS_score", "SE_score", "MXE_score"]] = merged_df_CS[["A3SS_score", "A5SS_score", "SE_score", "MXE_score"]].abs()

# get single |DPSI| score per flank
merged_df_CS['dPSI'] = merged_df_CS[["A3SS_score", "A5SS_score", "SE_score", "MXE_score"]].max(axis=1)

# check if column AS has larger dpsi scores across all flanks
# AS_lerger = (merged_df_AS['dPSI'] > merged_df_CS['dPSI']).all()

# remove CS flanks found in AS df (happens if flank is CS for one event)
merged_df = merged_df_CS.merge(merged_df_AS[['flanks', 'geneSymbol']], on=['flanks', 'geneSymbol'], how='left', indicator=True)
filtered_merged_df_CS = merged_df[merged_df['_merge'] == 'left_only'].drop('_merge', axis=1)

# remove RNA genes (eg TMX2-CTNND1)
filtered_merged_df_CS[['start', 'stop']] = filtered_merged_df_CS['flanks'].str.split('-', n=1, expand=True)
duplicates_mask = filtered_merged_df_CS.duplicated(subset=['chr', 'start', 'stop'], keep=False)
filtered_merged_df_CS = filtered_merged_df_CS[~(duplicates_mask & (filtered_merged_df_CS['geneSymbol'].str.contains("-", na=False)))]

# output files
filtered_merged_df_CS[['chr', 'start', 'stop']].drop_duplicates().to_csv(f'0_Files/filtered_flanks_CS.bed', index=False, sep='\t', header=False)
filtered_merged_df_CS[['flanks', 'geneSymbol', 'dPSI']].drop_duplicates().to_csv(f'0_Files/Filtered_dPSI_CS.csv', index=False, sep='\t')
