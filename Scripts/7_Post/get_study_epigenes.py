import pandas as pd
import json
import glob
from collections import Counter


with open('paths_multi.json') as f:
        data = json.load(f)

# STEP 1: Get the list of histone modifications available in the study
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

        
            