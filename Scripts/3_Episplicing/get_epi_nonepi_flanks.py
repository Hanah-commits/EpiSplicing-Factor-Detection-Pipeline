import pandas as pd
import json
import glob
from collections import Counter
from pathlib import Path
import os
import math
import sys

def get_epigenes_study(tool):

    prefix = tool.lower()

    # STEP 0: Create directories to store MAJIQ files
    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    with open('paths.json') as f:
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


    ## STEP 4: Get DHM-DEU values for all epigenes:
    opdirs = []
    for process in data['list_of_processes']:
        opdirs.append(data[process]['Output directory'])  

    DHM_vals_epi = []
    for dir in list(opdirs):
        try:
            file_path = glob.glob(f'{dir}*0_Files/{tool}/dPSI_Mval_epi_{prefix}.csv')
            DHM_vals_epi.append(pd.read_csv(file_path[0],delimiter='\t'))
        except:
            continue

    df_epi = pd.concat(DHM_vals_epi,axis=0,sort=False)

    epi_dfs = []


    ## STEP 5: Get AS exon flanks of hm-specific epigenes and non-epigenes
    for hm in list(all_hms):
        df = df_epi[df_epi.gene_name.isin(epigenes[hm])]
        df = df.drop_duplicates(subset=['idx'], keep='first').reset_index(drop=True)
        df = df[(df.dPSI != 0) & (df[hm] !=0)]
        df.to_csv(f'{output_dir}{hm}_epigenes.tsv', sep='\t', index=False)
        
        df[ 'type'] = hm
        epi_dfs.append(df)

    ## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

    # epi AS flanks into bed files

    df = pd.concat(epi_dfs,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    df.to_csv(f'{output_dir}epi_flanks.bed', sep='\t', index=False, header=False)

def get_nonepigenes(tool):

    prefix = tool.lower()

    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    with open('paths.json') as f:
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
        DHM_vals_nonepi = []

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

    ## STEP 4: Get DHM-DEU values for all non-epigenes:
        
    opdirs = []
    DHM_vals_nonepi = []
    nonepi_dfs = []

    for process in data['list_of_processes']:
        opdirs.append(data[process]['Output directory'])  

    for dir in list(opdirs):
        try:
            file_path = glob.glob(f'{dir}*0_Files/{tool}/dPSI_Mval_nonepi_{prefix}.csv')
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

        df[ 'type'] = hm
        nonepi_dfs.append(df)

    ## STEP 6: Prepare RBPmap input : Get AS exons of hm-specific epigenes and non-epigenes

    # epi and nonepi AS flanks separately into bed files
    df = pd.concat(nonepi_dfs,axis=0,sort=False)
    df = df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']]
    df = df.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    df.to_csv(f'{output_dir}nonepi_flanks.bed', sep='\t', index=False, header = False)

def get_epi_nonspliced_genes(tool):

    file = tool.lower()

    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    with open('paths.json') as f:
            data = json.load(f)

    # STEP 1: Get all non-DEU exons with DHM annotation (from constitutively-spliced genes)from pairwise analyses:
    hms = set()
    for process in data['list_of_processes']:
        hms.update(data[process]['Histone modifications'])

    all_nonDEU_dfs = []
    for process in data['list_of_processes']:
        
        odir = data[process]['Output directory']
        pr_hms = data[process]['Histone modifications']

        tissues = os.path.basename(os.path.dirname(odir)).split('_')[1:3]
        tissues_substr = ['endo', 'ecto', 'meso', 'neuro', 'H1']
        tissues = [next((sub for sub in tissues_substr if sub in tissue), tissue) for tissue in tissues]
        tissues = '-'.join(tissues).title()
        print("\n", tissues)

        flanks = pd.read_csv(f'{odir}0_Files/{tool}/DEU_DHM_{file}_flanks.tsv', delimiter='\t')

        print('Exon flanks (DEU !=0 & DHM !=0) before filtering out differentially spliced genes \n')
        for hm in pr_hms:
            print(hm, len(flanks[(flanks.dPSI == 0) & (flanks[hm]!=0)]))

        ## FILTER: keep only genes with non-DEU events
        flanks = flanks.groupby("gene_name").filter(lambda g:  (g["dPSI"] == 0).all())
        flanks = flanks[flanks[pr_hms].ne(0).any(axis=1)]

        print('Exon flanks  (DEU !=0 & DHM !=0) after filtering out differentially spliced genes \n')
        for hm in pr_hms:
            print(hm, len(flanks[flanks[hm]!=0]))

        all_nonDEU_dfs.append(flanks)

    epi_nonDEU_flanks = pd.concat(all_nonDEU_dfs).drop_duplicates().fillna(0)
    if len(epi_nonDEU_flanks) > 0:
        epi_nonDEU_flanks['type'] = epi_nonDEU_flanks[list(hms)].apply(lambda row: ",".join(row.index[row.ne(0)]), axis=1)
        epi_nonDEU_flanks['type'] = epi_nonDEU_flanks['type'].str.split(",")
        epi_nonDEU_flanks = epi_nonDEU_flanks.explode('type', ignore_index=True)

    else:
        print('\n No Exon flanks  (DEU !=0 & DHM !=0)  combined across all 10 analyses with DHM!=0 \n')

    # STEP 2: Save flanks
    epi_nonDEU_flanks = epi_nonDEU_flanks[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']].drop_duplicates()
    
    print('Exon flanks  (DEU !=0 & DHM !=0)  combined across all 10 analyses:\n')
    print(epi_nonDEU_flanks['type'].value_counts())

    epi_nonDEU_flanks = epi_nonDEU_flanks.groupby(['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name'])['type'].apply(','.join).reset_index()
    op_dir = '0_Files/Post-processing/'
    epi_nonDEU_flanks.to_csv(f'{output_dir}epi_nonspliced_flanks.bed', sep='\t', index=False, header=False)

if __name__ == "__main__":
    tool = "RMATS"

    get_epigenes_study(tool) # DEU !=0 & DHM !=0
    get_nonepigenes(tool) # DEU != 0 & DHM == 0
    get_epi_nonspliced_genes(tool) # DEU ==0 & DHM != 0
