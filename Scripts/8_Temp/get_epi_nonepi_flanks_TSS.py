import pandas as pd
import json
import glob
from collections import Counter
from pathlib import Path
import os
import math
import sys


def get_epigenes_study(tool):

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
            file_path = glob.glob(f'{dir}*0_Files/{tool}/{hm}/{hm}_truepos_epigenes.txt')
            try: ## ../0_Files/{tool}/H3K27ac/H3K27ac_truepos_epigenes.txt
                with open(file_path[0], 'r') as file:
                    epi= [gene.strip() for gene in file]                
                    hm_epigenes.extend(epi)
                    hm_epigenes= list(set(hm_epigenes))
            except: ## no epigenes detected by tool for current HM
                continue

        
        ## STEP 3: TODO: Plot frequency of occurrence of genes across all ten conditions
        
        
        print('Epigenes:\n', hm_epigenes)
        epigenes[hm] = hm_epigenes


    ## STEP 4: Get TSS_exons of all epigenes:

    df_epi = pd.read_csv('0_Files/Post-processing/TSS_flanks200.bed', delimiter='\t', header=None)
    df_epi.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_id', 'gene_name']

    epi_dfs = []


    ## STEP 5: Get TSS exons of hm-specific epigenes and non-epigenes
    for hm in list(all_hms):
        df = df_epi[df_epi.gene_name.isin(epigenes[hm])]
        
        df.loc[:, 'type'] = hm
        epi_dfs.append(df)

    ## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

    # epi AS flanks into bed files

    df = pd.concat(epi_dfs,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    df.to_csv(f'{output_dir}epi_flanks.bed', sep='\t', index=False, header=False)

def get_nonepigenes(tool):

    prefix = tool.lower()

    output_dir = str(Path(os.getcwd())) + "/Post-processing/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    with open('paths_multi.json') as f:
        data = json.load(f)

    # STEP 1: Get the HMs available for current comparison
    hms = set()
    for process in data['list_of_processes']:
        hms.update(data[process]['Histone modifications'])

    # STEP 2: Get the nonepigenes detected by the tool for all the HMs
    nonepigenes = {}
    for hm in list(hms):
        print('\n', hm)
        # Get the list of output directories for the current hm
        processes = [process for process in data['list_of_processes'] if hm in data[process]['Histone modifications']]
        output_directories = [data[process]['Output directory'] for process in processes]

        hm_nonepigenes = []

        for dir in output_directories:
            file_path = glob.glob(f'{dir}*0_Files/{tool}/{prefix}_nonepigenes.txt')
            try: ## ../0_Files/{tool}/dexseq_nonepigenes.txt
                with open(file_path[0], 'r') as file:
                    nonepi= list(set([gene.strip() for gene in file]))                
                    hm_nonepigenes.extend(nonepi)
            except: ## no nonepigenes detected by current tool for current HM
                continue

        ## STEP 3: FINAL NONEPIGENES
        # Count occurrences of nonepigenes 
        nonepigene_counts = Counter(gene for gene in hm_nonepigenes)

        # Filter genes detected in at least 80% of the analyses that have chipseq data for this hm
        overlap_nonepigenes = [gene for gene, count in nonepigene_counts.items() if count >= math.floor(len(output_directories) * 0.8)]
        
        print('Non-epigenes:\n', len(overlap_nonepigenes))
        nonepigenes[hm] = overlap_nonepigenes



    ## STEP 4: Get TSS_exons of all non-epigenes:

    df_nonepi = pd.read_csv('0_Files/Post-processing/TSS_flanks200.bed', delimiter='\t')
    df_nonepi.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_id', 'gene_name']

    nonepi_dfs = []


    ## STEP 5: Get TSS exon flanks of hm-specific epigenes and non-epigenes
    for hm in list(hms):
        df = df_nonepi[df_nonepi.gene_name.isin(nonepigenes[hm])]
        
        df.loc[:, 'type'] = hm
        nonepi_dfs.append(df)

    ## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

    # epi and nonepi AS flanks separately into bed files
    df = pd.concat(nonepi_dfs,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    df.to_csv(f'{output_dir}nonepi_flanks.bed', sep='\t', index=False, header = False)


if __name__ == "__main__":
    tool = sys.argv[1]
    type = sys.argv[2]

    if type == 'epi':
        get_epigenes_study(tool)
    else:
        get_nonepigenes(tool)