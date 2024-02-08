import pandas as pd
import json

with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]

deu_flanks = pd.read_csv('0_Files/RMATS/rmats_flanks200.bed', delimiter='\t', header=None)
deu_flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'geneSymbol', 'dPSI']
deu_flanks.drop(columns=['feature', 'score'], inplace=True)
deu_flanks.drop_duplicates(inplace=True)

dhm_flanks = pd.read_csv('0_Files/MANorm/DHM_peaks_annotation.tsv', delimiter='\t')

# combine DEU and DHM scores
dhm_flanks = dhm_flanks.merge(deu_flanks, on=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'], how='left')

# Fill NaN values with 0 
dhm_flanks['dPSI'] = dhm_flanks['dPSI'].fillna(0)

# keep only genes where DEU scores are available
deu_genes = list(set(deu_flanks.geneSymbol.values.tolist()))
dhm_flanks = dhm_flanks[dhm_flanks.geneSymbol.isin(deu_genes)]

# save
dhm_flanks.to_csv('0_Files/RMATS/DEU_DHM_rmats_flanks.tsv', sep='\t', index=False)
