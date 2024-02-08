import pandas as pd
import json

with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]

deu_flanks = pd.read_csv('0_Files/DEXSEQ/dexseq_flanks200.bed', delimiter='\t', header=None)
deu_flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'geneSymbol', 'dPSI']
deu_flanks.drop(columns=['feature', 'score'], inplace=True)
deu_flanks.drop_duplicates(inplace=True)

dhm_flanks = pd.read_csv('0_Files/MANorm/DHM_peaks_annotation.tsv', delimiter='\t')

## STEP 0: deal with flank have 1+ dpsi score (flanks of A3SS/A5SS exons -> duplicate flanks with dPSI scores of diff exons)

deu_flanks['dPSI_abs'] = pd.to_numeric(deu_flanks['dPSI']).abs()

# get all peaks that belong to each flank
flanks_dpsi_group = deu_flanks.groupby(['chr', 'flank_start', 'flank_end', 'geneSymbol'])['dPSI_abs'] \
    .apply(lambda val: ','.join(str(v) for v in val)).reset_index()

## FILTER 1: if flank has 1+ dPSI, keep deu with highest abs dPSI
flanks_dpsi_group['dPSI_abs'] = flanks_dpsi_group['dPSI_abs'].str.split(',')  # string -> list of strings
flanks_dpsi_group['max_dPSI'] = flanks_dpsi_group['dPSI_abs'].apply(lambda x: max(map(float, x)))  # max dPSI

# get the corresponding DEU for each flank's max dPSI
flanks_dpsi_group = pd.merge(flanks_dpsi_group[['flank_start', 'flank_end', 'max_dPSI']], deu_flanks, on=['flank_start', 'flank_end'],
                            how='inner')
flanks_dpsi_group = flanks_dpsi_group[flanks_dpsi_group['dPSI_abs'] == flanks_dpsi_group['max_dPSI']]

## FILTER 2: If flank has 1+ peaks with same max |dPSI|, keep one
flanks_dpsi_group.drop_duplicates(subset=['chr', 'flank_start', 'flank_end', 'geneSymbol'], keep='first', inplace=True)


## STEP 1: combine DEU and DHM scores
dhm_flanks = dhm_flanks.merge(flanks_dpsi_group, on=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'], how='left')

# Fill NaN values with 0 
dhm_flanks['dPSI'] = dhm_flanks['dPSI'].fillna(0)
dhm_flanks.drop_duplicates(inplace=True)


## STEP 2: keep only genes where DEU scores are available
deu_genes = list(set(deu_flanks.geneSymbol.values.tolist()))
dhm_flanks = dhm_flanks[dhm_flanks.geneSymbol.isin(deu_genes)]

## save
dhm_flanks.to_csv('0_Files/DEXSEQ/DEU_DHM_dexseq_flanks.tsv', sep='\t', index=False)