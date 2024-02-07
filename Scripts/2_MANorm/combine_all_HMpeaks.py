import pandas as pd
import json


with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]


concat = []
all_flanks = pd.read_csv('0_Files/flanks200.bed', delimiter='\t', header=None)
all_flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'geneSymbol']

peak_dfs = []
for hm in hms:
    flank_peaks = pd.read_csv(f'0_Files/MANorm/{hm}_annotated_flanks.bed', delimiter='\t', header=None)
    flank_peaks.columns = ['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol', 'peak_start', 'peak_end', 'peak_feature', hm]

    hm_flanks = all_flanks.merge(flank_peaks, on=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'], how='left')

    # Fill NaN values with 0 
    hm_flanks[['peak_start', 'peak_end', 'peak_feature',  hm]] = hm_flanks[['peak_start', 'peak_end', 'peak_feature',  hm]].fillna(0)

    hm_flanks = hm_flanks[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'geneSymbol', hm]]
    peak_dfs.append(hm_flanks)


# peak_dfs = [df.set_index('flanks') for df in peak_dfs]
peak_dfs = pd.concat(peak_dfs, axis=1)
peak_dfs = peak_dfs.loc[~(peak_dfs==0).all(axis=1)]

# Identify and keep only the first occurrence of each column name
unique_columns = ~peak_dfs.columns.duplicated(keep='first')
peak_dfs = peak_dfs.loc[:, unique_columns]
