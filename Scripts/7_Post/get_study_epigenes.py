import pandas as pd
import json
import glob
from collections import Counter


with open('paths_multi.json') as f:
        data = json.load(f)

# STEP 1: Get the list of histone modifications available in the study
all_epigenes = []
all_hms = set()
for process in data['list_of_processes']:
    all_hms.update(data[process]['Histone modifications'])

# STEP 2: Get the epigenes detected by the DEU-TOOLS for all the HMs
for hm in list(all_hms):
    print(hm, '\n\n')
    # Get the list of output directories for the current hm
    processes = [process for process in data['list_of_processes'] if hm in data[process]['Histone modifications']]
    output_directories = [data[process]['Output directory'] for process in processes]

    hm_epigenes = {
        'DEXSEQ': [],
        'MAJIQ': [],
        'RMATS': []
    }

    ## 
    for dir in output_directories:
        for tool in ['DEXSEQ', 'MAJIQ', 'RMATS']:
            file_path = glob.glob(f'{dir}*0_Files/{tool}/{hm}/{hm}_truepos_epigenes.txt')
            try: ## ../0_Files/DEXSEQ/H3K27ac/H3K27ac_truepos_epigenes.txt
                with open(file_path[0], 'r') as file:
                    epigenes= [gene.strip() for gene in file]                
                    hm_epigenes[tool].extend(epigenes)
                    hm_epigenes[tool] = list(set(hm_epigenes[tool]))
            except: ## no epigenes detected by current tool for current HM
                continue
    
    ## STEP 3: FINAL EPIGENES
    # Count occurrences of genes across all three tools
    epigene_counts = Counter(value for sublist in hm_epigenes.values() for value in sublist)

    # Filter genes that appear in all three tools
    overlap_epigenes = [gene for gene, count in epigene_counts.items() if count >= 2]
    

    print(overlap_epigenes)
    all_epigenes.append(overlap_epigenes)

## STEP 4: Get DHM-DEU values for all epigenes:
opdirs = []
for process in data['list_of_processes']:
    opdirs.append(data[process]['Output directory'])  

DHM_vals = []
for dir in list(opdirs):

    for tool in ['DEXSEQ', 'MAJIQ', 'RMATS']:
        try:
            file = tool.lower()
            file_path = glob.glob(f'{dir}*0_Files/{tool}/dPSI_Mval_epi_{file}.csv')
            DHM_vals.append(pd.read_csv(file_path[0],delimiter='\t'))
        except:
            continue

all_epigenes = list(set([item for sublist in all_epigenes for item in sublist]))

df = pd.concat(DHM_vals,axis=0,sort=False)
df = df[df['gene_name'].isin(all_epigenes)].reset_index(drop=True)
df = df.drop_duplicates(subset=['idx', 'dPSI'], keep='first').reset_index(drop=True)
df.fillna(0, inplace=True)
df.to_csv('final_Epigenes.tsv', sep='\t', index=False)


