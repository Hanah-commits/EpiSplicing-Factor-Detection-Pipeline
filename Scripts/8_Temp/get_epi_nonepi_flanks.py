import pandas as pd
import json
import glob
from collections import Counter
from pathlib import Path
import os

# STEP 0: Create directories to store MAJIQ files
output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

# STEP 1: Get the HMs available for current comparison
with open('paths.json') as f:
        d = json.load(f)

hms = d["Histone modifications"]

# STEP 2: Get the epigenes detected by the DEU-TOOLS for all the HMs
hm_epigenes = {}
hm_nonepigenes = {}
for hm in hms:
    print('\n', hm)

    # get epigenes
    file_path = glob.glob(f'0_Files/RMATS/{hm}/{hm}_truepos_epigenes.txt')
    try: ## ../0_Files/RMATS/H3K27ac/H3K27ac_truepos_epigenes.txt
        with open(file_path[0], 'r') as file:
            epi= [gene.strip() for gene in file]               
            hm_epigenes[hm] = list(set(epi))
    except: ## no epigenes detected by current tool for current HM
        continue

    # get nonepigenes
    file_path = glob.glob(f'0_Files/RMATS/rmats_nonepigenes.txt')
    try: ## ../0_Files/DEXSEQ/dexseq_nonepigenes.txt
        with open(file_path[0], 'r') as file:
            nonepi= [gene.strip() for gene in file]                
            hm_nonepigenes[hm] = list(set(nonepi))
    except: ## no nonepigenes detected by current tool for current HM
        continue

## STEP 4a: Get DHM-DEU values for all epigenes
df_epi = pd.read_csv('0_Files/RMATS/dPSI_Mval_epi_rmats.csv', delimiter='\t')
df_nonepi = pd.read_csv('0_Files/RMATS/dPSI_Mval_nonepi_rmats.csv', delimiter='\t')

epi_dfs = []
nonepi_dfs = []

## STEP 5: Get AS exon flanks of hm-specific epigenes and non-epigenes
for hm in hms:
    df = df_epi[df_epi.gene_name.isin(hm_epigenes[hm])]
    df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
    df = df[df.dPSI != 0]
    # df.to_csv(f'{output_dir}{hm}_epigenes.tsv', sep='\t', index=False)
    
    df['type'] = hm
    epi_dfs.append(df)

    df = df_nonepi[df_nonepi.gene_name.isin(hm_nonepigenes[hm])]
    df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
    df = df[df.dPSI != 0]
    # df.to_csv(f'{output_dir}{hm}_nonepigenes.tsv', sep='\t', index=False)

    df['type'] = hm
    nonepi_dfs.append(df)

## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

# epi and nonepi AS flanks separately into bed files
i = 0
for df in [epi_dfs, nonepi_dfs]:
    df = pd.concat(df,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    
    if i == 0:
        df.to_csv(f'{output_dir}epi_flanks.bed', sep='\t', index=False, header=False)
    else:
        df.to_csv(f'{output_dir}nonepi_flanks.bed', sep='\t', index=False, header = False)

    i += 1