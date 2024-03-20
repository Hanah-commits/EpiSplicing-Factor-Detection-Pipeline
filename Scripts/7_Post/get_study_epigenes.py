import pandas as pd
import json
import glob
from collections import Counter
from pathlib import Path
import os


# STEP 0: Create directories to store MAJIQ files
output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

with open('paths_multi.json') as f:
        data = json.load(f)

# STEP 1: Get the list of histone modifications available in the study
epigenes = {}
nonepigenes = {}

all_hms = set()
for process in data['list_of_processes']:
    all_hms.update(data[process]['Histone modifications'])

# STEP 2: Get the epigenes detected by the DEU-TOOLS for all the HMs
for hm in list(all_hms):
    print('\n', hm)
    # Get the list of output directories for the current hm
    processes = [process for process in data['list_of_processes'] if hm in data[process]['Histone modifications']]
    output_directories = [data[process]['Output directory'] for process in processes]

    hm_epigenes = {
        'DEXSEQ': [],
        'MAJIQ': [],
        'RMATS': []
    }

    hm_nonepigenes = {
        'DEXSEQ': [],
        'MAJIQ': [],
        'RMATS': []        
    }

    ## 
    for dir in output_directories:
        for tool in ['DEXSEQ', 'MAJIQ', 'RMATS']:

            # get epigenes
            file_path = glob.glob(f'{dir}*0_Files/{tool}/{hm}/{hm}_truepos_epigenes.txt')
            try: ## ../0_Files/DEXSEQ/H3K27ac/H3K27ac_truepos_epigenes.txt
                with open(file_path[0], 'r') as file:
                    epi= [gene.strip() for gene in file]                
                    hm_epigenes[tool].extend(epi)
                    hm_epigenes[tool] = list(set(hm_epigenes[tool]))
            except: ## no epigenes detected by current tool for current HM
                continue

            # get nonepigenes
            file = tool.lower()
            file_path = glob.glob(f'{dir}*0_Files/{tool}/{file}_nonepigenes.txt')
            try: ## ../0_Files/DEXSEQ/dexseq_nonepigenes.txt
                with open(file_path[0], 'r') as file:
                    nonepi= [gene.strip() for gene in file]                
                    hm_nonepigenes[tool].extend(nonepi)
                    hm_nonepigenes[tool] = list(set(hm_nonepigenes[tool]))
            except: ## no nonepigenes detected by current tool for current HM
                continue

    
    ## STEP 3a: FINAL EPIGENES
    # Count occurrences of genes across all three tools
    epigene_counts = Counter(value for sublist in hm_epigenes.values() for value in sublist)

    # Filter genes that appear in all three tools
    overlap_epigenes = [gene for gene, count in epigene_counts.items() if count >= 2]
    

    print('Epigenes:\n', overlap_epigenes)
    epigenes[hm] = overlap_epigenes

    ## STEP 3b: FINAL NON-EPIGENES
    # Count occurrences of genes across all three tools
    nonepigene_counts = Counter(value for sublist in hm_nonepigenes.values() for value in sublist)

    # Filter genes that appear in all three tools
    overlap_nonepigenes = [gene for gene, count in nonepigene_counts.items() if count > 2]
    

    print('Non-epigenes:\n', len(overlap_nonepigenes))
    nonepigenes[hm] = overlap_nonepigenes

## STEP 4a: Get DHM-DEU values for all epigenes:
opdirs = []
for process in data['list_of_processes']:
    opdirs.append(data[process]['Output directory'])  

DHM_vals_epi = []
for dir in list(opdirs):

    for tool in ['DEXSEQ', 'MAJIQ', 'RMATS']:
        try:
            file = tool.lower()
            file_path = glob.glob(f'{dir}*0_Files/{tool}/dPSI_Mval_epi_{file}.csv')
            DHM_vals_epi.append(pd.read_csv(file_path[0],delimiter='\t'))
        except:
            continue

## STEP 4b: Get DHM-DEU values for all non-epigenes:
DHM_vals_nonepi = []
for dir in list(opdirs):

    for tool in ['DEXSEQ', 'MAJIQ', 'RMATS']:
        try:
            file = tool.lower()
            file_path = glob.glob(f'{dir}*0_Files/{tool}/dPSI_Mval_nonepi_{file}.csv')
            DHM_vals_nonepi.append(pd.read_csv(file_path[0],delimiter='\t'))
        except:
            continue

df_epi = pd.concat(DHM_vals_epi,axis=0,sort=False)
df_nonepi = pd.concat(DHM_vals_nonepi,axis=0,sort=False)

epi_dfs = []
nonepi_dfs = []

## STEP 5: Get AS exon flanks of hm-specific epigenes and non-epigenes
for hm in list(all_hms):
    df = df_epi[df_epi.gene_name.isin(epigenes[hm])]
    df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
    df = df[df.dPSI != 0]
    df.to_csv(f'{output_dir}{hm}_epigenes.tsv', sep='\t', index=False)
    
    df['type'] = hm
    epi_dfs.append(df)

    df = df_nonepi[df_nonepi.gene_name.isin(nonepigenes[hm])]
    df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
    df = df[df.dPSI != 0]
    df.to_csv(f'{output_dir}{hm}_nonepigenes.tsv', sep='\t', index=False)

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