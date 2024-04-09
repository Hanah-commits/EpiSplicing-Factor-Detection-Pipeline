import pandas as pd
import json
import glob
from pathlib import Path
import os


# STEP 0: Create directories to store MAJIQ files
output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

with open('paths_multi.json') as f:
        data = json.load(f)

# STEP 1: Get the list of histone modifications available in the study
epigenes = {}

all_hms = set()
for process in data['list_of_processes']:
    all_hms.update(data[process]['Histone modifications'])

# STEP 2: Get the epigenes detected by the DEU-TOOLS for all the HMs
for hm in list(all_hms):
    print('\n', hm)
    # Get the list of output directories for the current hm
    processes = [process for process in data['list_of_processes'] if hm in data[process]['Histone modifications']]
    output_directories = [data[process]['Output directory'] for process in processes]

    hm_epigenes = []

    ## 
    for dir in output_directories:
        # get epigenes
        file_path = glob.glob(f'{dir}*0_Files/RMATS/{hm}/{hm}_truepos_epigenes.txt')
        try: ## ../0_Files/RMATS/H3K27ac/H3K27ac_truepos_epigenes.txt
            with open(file_path[0], 'r') as file:
                epi= [gene.strip() for gene in file]                
                hm_epigenes.extend(epi)
                hm_epigenes= list(set(hm_epigenes))
        except: ## no epigenes detected by RMATS for current HM
            continue

    
    ## STEP 3: TODO: Plot frequency of occurrence of genes across all ten conditions
    
    
    print('Epigenes:\n', hm_epigenes)
    epigenes[hm] = hm_epigenes


## STEP 4: Get DHM-DEU values for all epigenes:
opdirs = []
for process in data['list_of_processes']:
    opdirs.append(data[process]['Output directory'])  

DHM_vals_epi = []
for dir in list(opdirs):
    try:
        file_path = glob.glob(f'{dir}*0_Files/RMATS/dPSI_Mval_epi_rmats.csv')
        DHM_vals_epi.append(pd.read_csv(file_path[0],delimiter='\t'))
    except:
        continue

df_epi = pd.concat(DHM_vals_epi,axis=0,sort=False)

epi_dfs = []


## STEP 5: Get AS exon flanks of hm-specific epigenes and non-epigenes
for hm in list(all_hms):
    df = df_epi[df_epi.gene_name.isin(epigenes[hm])]
    df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
    df = df[df.dPSI != 0]
    df.to_csv(f'{output_dir}{hm}_epigenes.tsv', sep='\t', index=False)
    
    df['type'] = hm
    epi_dfs.append(df)

## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

# epi AS flanks into bed files

df = pd.concat(epi_dfs,axis=0,sort=False)
df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
df.to_csv(f'{output_dir}epi_flanks.bed', sep='\t', index=False, header=False)
