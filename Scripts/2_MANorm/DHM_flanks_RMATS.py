import pandas as pd
import json
import sys

with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]

try:
    deu_flanks = pd.read_csv('0_Files/RMATS/rmats_flanks200.bed', delimiter='\t', header=None)
except:
    print('No RMATS exons to annotate \n \n')
    sys.exit()

deu_flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'geneSymbol', 'dPSI']
deu_flanks.drop(columns=['feature', 'score'], inplace=True)
deu_flanks.drop_duplicates(inplace=True)

print ('Annotating RMATS exons with HM peaks \n')
print('TSS Filtering:                   ', len(set(deu_flanks.geneSymbol.values.tolist()))) # log

dhm_flanks = pd.read_csv('0_Files/MANorm/DHM_peaks_annotation.tsv', delimiter='\t')

## STEP 1: combine DEU and DHM scores
dhm_flanks = dhm_flanks.merge(deu_flanks, on=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'], how='left')

# Fill NaN values with 0 
dhm_flanks['dPSI'] = dhm_flanks['dPSI'].fillna(0)

## STEP 2:  keep only genes where DEU scores are available
deu_genes = list(set(deu_flanks.geneSymbol.values.tolist()))
dhm_flanks = dhm_flanks[dhm_flanks.geneSymbol.isin(deu_genes)]

print('TSL Filtering:                   ', len(set(dhm_flanks.geneSymbol.values.tolist())), '\n \n') # log

# # save
if len(dhm_flanks) > 0:
    dhm_flanks.to_csv('0_Files/RMATS/DEU_DHM_rmats_flanks.tsv', sep='\t', index=False)
