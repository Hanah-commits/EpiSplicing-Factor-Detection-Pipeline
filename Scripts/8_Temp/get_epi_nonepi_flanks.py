import pandas as pd
import json
import glob
from collections import Counter
from pathlib import Path
import os
import math


def get_epigenes():

    # STEP 0: Create directories to store MAJIQ files
    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # STEP 1: Get the HMs available for current comparison
    with open('paths.json') as f:
            d = json.load(f)

    hms = d["Histone modifications"]

    # STEP 2: Get the epigenes detected by the DEU-TOOLS for all the HMs
    hm_epigenes = {}
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

    ## STEP 4a: Get DHM-DEU values for all epigenes
    df_epi = pd.read_csv('0_Files/RMATS/dPSI_Mval_epi_rmats.csv', delimiter='\t')

    epi_dfs = []

    ## STEP 5: Get AS exon flanks of hm-specific epigenes and non-epigenes
    for hm in hms:
        df = df_epi[df_epi.gene_name.isin(hm_epigenes[hm])]
        df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
        df = df[df.dPSI != 0]
        # df.to_csv(f'{output_dir}{hm}_epigenes.tsv', sep='\t', index=False)
        
        df['type'] = hm
        epi_dfs.append(df)

    ## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes
    # epi AS flanks into bed files
    df = pd.concat(epi_dfs,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    df.to_csv(f'{output_dir}epi_flanks.bed', sep='\t', index=False, header=False)


def get_nonepigenes():

    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"

    with open('paths_multi.json') as f:
        data = json.load(f)

    # STEP 1: Get the HMs available for current comparison
    with open('paths.json') as f:
            d = json.load(f)

    hms = d["Histone modifications"]

    # STEP 2: Get the nonepigenes detected by the RMATS for all the HMs
    nonepigenes = {}
    for hm in hms:
        print('\n', hm)
        # Get the list of output directories for the current hm
        processes = [process for process in data['list_of_processes'] if hm in data[process]['Histone modifications']]
        output_directories = [data[process]['Output directory'] for process in processes]

        hm_nonepigenes = []
        DHM_vals_nonepi = []

        for dir in output_directories:
            file_path = glob.glob(f'{dir}*0_Files/RMATS/rmats_nonepigenes.txt')
            try: ## ../0_Files/RMATS/dexseq_nonepigenes.txt
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

    ## STEP 4: Get DHM-DEU values for all non-epigenes:
        
    opdirs = []
    DHM_vals_nonepi = []
    nonepi_dfs = []

    for process in data['list_of_processes']:
        opdirs.append(data[process]['Output directory'])  

    for dir in list(opdirs):
        try:
            file_path = glob.glob(f'{dir}*0_Files/RMATS/dPSI_Mval_nonepi_rmats.csv')
            DHM_vals_nonepi.append(pd.read_csv(file_path[0],delimiter='\t'))
        except:
            continue

    df_nonepi = pd.concat(DHM_vals_nonepi,axis=0,sort=False)

    ## STEP 5: Get AS exon flanks of hm-specific epigenes and non-epigenes
    for hm in hms:

        df = df_nonepi[df_nonepi.gene_name.isin(nonepigenes[hm])]
        df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
        df = df[df.dPSI != 0]
        # df.to_csv(f'{output_dir}{hm}_nonepigenes.tsv', sep='\t', index=False)

        df['type'] = hm
        nonepi_dfs.append(df)

    ## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

    # epi and nonepi AS flanks separately into bed files
    df = pd.concat(nonepi_dfs,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    df.to_csv(f'{output_dir}nonepi_flanks.bed', sep='\t', index=False, header = False)


get_epigenes()
get_nonepigenes()