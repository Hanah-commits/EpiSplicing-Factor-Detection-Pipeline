import pandas as pd
import numpy as np
import json
import sys


with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]


peaksfiles = [f'0_Files/MANorm/{hm}_annotated_flanks.bed' for hm in hms] 
peak_dfs = []

flanks = pd.read_csv(f'0_Files/DEXSEQ/dexseq_flanks200.bed', delimiter='\t', header=None)
flanks.columns = ['chr', "flank_start", "flank_stop", "feature", "score", "strand", "geneSymbol", 'dPSI']
flanks['flanks'] = flanks[['flank_start', 'flank_stop']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
flanks.drop_duplicates(inplace=True)

for file in peaksfiles:

    hm = file.split('_annotated')[0]
    peaks = pd.read_csv(file, delimiter='\t', header=None)
    peaks.drop([3, 4, 7, 10], axis=1, inplace=True)
    
    # assign 0 to flanks that have no peaks
    peaks.replace([-1, '.'], [0, 0], inplace=True)
    # peaks[3] = peaks[[1, 2]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    peaks.columns = ['chr',  "flank_start", "flank_stop", "strand", "geneSymbol", 'peak_start', 'peak_stop', 'M_value', 'overlap_bp']

    peaks['M_value_abs'] = pd.to_numeric(peaks['M_value']).abs()

    # get all peaks that belong to each flank
    flank_peaks_group = peaks.groupby(['chr', 'flank_start', 'flank_stop'])['M_value_abs'] \
        .apply(lambda val: ','.join(str(v) for v in val)).reset_index()

    # # FILTER 1: if flank has 1+ peaks, keep peak with highest abs M-value
    flank_peaks_group['M_value_abs'] = flank_peaks_group['M_value_abs'].str.split(',')  # string -> list of strings
    flank_peaks_group['max_' + hm] = flank_peaks_group['M_value_abs'].apply(lambda x: max(map(float, x)))  # max MValue
    # # flank_peaks_group['#peaks_'+hm] = flank_peaks_group['M_value'].apply(lambda x: len(x))  # no of peaks/ flank

    # # get the corresponding peak for each flank's max M-value
    flank_peaks_group = pd.merge(flank_peaks_group[['flank_start', 'flank_stop', 'max_' + hm]], peaks, on=['flank_start', 'flank_stop'],
                                how='inner')
    # # flank_peaks_group = pd.merge(flank_peaks_group[['flanks', 'max_'+hm, '#peaks_'+hm]], peaks, on=['flanks'], how='inner')
    flank_peaks_group = flank_peaks_group[flank_peaks_group['M_value_abs'] == flank_peaks_group['max_' + hm]]

    # FILTER 2: If flank has 1+ peaks with same max |Mvalue|, keep one
    flank_peaks_group.drop_duplicates(subset=['flank_start', 'flank_stop'], keep='first', inplace=True)

    flank_peaks_group.rename(columns={'M_value': hm}, inplace=True)
    peak_dfs.append(flank_peaks_group[['chr', 'flanks', 'geneSymbol', 'strand', 'dPSI', hm]]) 

peak_dfs = [df.set_index('flanks') for df in peak_dfs]
peak_dfs = pd.concat(peak_dfs, axis=1)

# Identify and keep only the first occurrence of each column name
unique_columns = ~peak_dfs.columns.duplicated(keep='first')
peak_dfs = peak_dfs.loc[:, unique_columns]

# no peak for hm A in flank
peak_dfs.fillna(0, inplace=True)

peak_dfs.to_csv(f'0_Files/DEXSEQ/Filtered_MValues_dexseq.csv', sep='\t')
