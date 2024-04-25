import pandas as pd
import json
import glob
from collections import Counter
from pathlib import Path
import os
import math
import sys

def get_TSS_exons(tool):

    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/"

    with open('paths_multi.json') as f:
            data = json.load(f)


    output_directories = [data[process]['Output directory'] for process in data['list_of_processes']]
    fasta = list(set([data[process]['Reference fasta'] for process in data['list_of_processes']]))[0]
    ref_genome= fasta+".fai"

    i = 0
    for dir in output_directories:

        if tool == 'RMATS':

            # exon boundary external flanks
            os.system(f"bedtools flank -i {dir}0_Files/RMATS/rmats_exons_coords.bed -g  {ref_genome} -b 200 > {output_dir}flanks.bed" )

            # separate start,stop flank coords
            os.system(f"sed -n 'n;p' {output_dir}flanks.bed > {output_dir}stop.bed")
            os.system(f"sed -n 'p;n' {output_dir}flanks.bed > {output_dir}start.bed")

            # exon boundary internal flanks
            os.system(f"bedtools slop -i {output_dir}start.bed -g  {ref_genome} -l 0 -r 200 > {output_dir}start_flanks.bed")
            os.system(f"bedtools slop -i {output_dir}stop.bed -g  {ref_genome} -l 200 -r 0 > {output_dir}stop_flanks.bed")

            # combine start,stop flank coords
            os.system(f"paste -d'\n' {output_dir}start_flanks.bed {output_dir}stop_flanks.bed | sort -k1,1 -k2,2n > {output_dir}AS_flanks200_pr{i}.bed")

            # remove intermediate files
            os.system(f"rm {output_dir}start*.bed")
            os.system(f"rm  {output_dir}stop*.bed")
            os.system(f"rm {output_dir}flanks.bed")

        elif tool == 'MAJIQ':

            AS_flanks = pd.read_csv(f'{dir}0_Files/MAJIQ/Filtered_dPSI.csv', delimiter='\t')

            ## STEP 1: fetch CS exons
            # add label
            AS_flanks['feature'] = 'AS'
            AS_flanks['score'] = '.'

            AS_flanks = AS_flanks.drop_duplicates(subset=['seqid', 'start', 'stop', 'gene_id'])

            #@ STEP 2: Get exon flanks

            AS_flanks[['seqid', "start", "stop", "feature", "score", "strand", "gene_id", "mean_dpsi_per_lsv_junction"]].to_csv(f'{output_dir}AS_flanks200_pr{i}.bed', index=False, sep='\t', header=False)

        elif tool == 'DEXSEQ':

            # exon boundary external flanks
            os.system(f"bedtools flank -i {dir}0_Files/DEXSEQ/dexseq_exons_coords.bed -g {ref_genome} -b 200 >{output_dir}flanks.bed" )

            # separate start,stop flank coords
            os.system(f"sed -n 'n;p' {output_dir}flanks.bed > {output_dir}stop.bed")
            os.system(f"sed -n 'p;n' {output_dir}flanks.bed > {output_dir}start.bed")

            # exon boundary internal flanks
            os.system(f"bedtools slop -i {output_dir}start.bed -g {ref_genome} -l 0 -r 200 > {output_dir}start_flanks.bed")
            os.system(f"bedtools slop -i {output_dir}stop.bed -g {ref_genome} -l 200 -r 0 > {output_dir}stop_flanks.bed")

            # combine start,stop flank coords
            os.system(f"paste -d'\n' {output_dir}start_flanks.bed {output_dir}stop_flanks.bed | sort -k1,1 -k2,2n > {output_dir}AS_flanks200_pr{i}.bed")

            # remove intermediate files
            os.system(f"rm {output_dir}start*.bed")
            os.system(f"rm {output_dir}stop*.bed")
            os.system(f"rm {output_dir}flanks.bed")

        i += 1


    # combine AS flanks from all pairwise analyses into one
    os.system(f'cat {output_dir}AS_flanks200_pr*.bed | sort | uniq > {output_dir}AS_flanks200.bed')

    # keep only flanks overlapping with TSS exon flanks
    os.system(f'bedtools intersect -wa -a {output_dir}AS_flanks200.bed -b {output_dir}TSS_flanks200.bed -s | sort | uniq > {output_dir}{tool}TSS_flanks200.bed ')

    os.system(f'rm {output_dir}AS_flanks200*.bed')


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

    ## STEP 3: Get Alternative TSS_exons of all epigenes:

    df_epi = pd.read_csv(f'0_Files/Post-processing/{tool}TSS_flanks200.bed', delimiter='\t', header=None)
    df_epi.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_id', 'gene_name']

    epi_dfs = []


    ## STEP 5: Get TSS exons of hm-specific epigenes and non-epigenes
    for hm in list(all_hms):
        df = df_epi[df_epi.gene_name.isin(epigenes[hm])]
        
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


    ## STEP 4: Get alternative TSS_exons of all non-epigenes:

    df_nonepi = pd.read_csv(f'0_Files/Post-processing/{tool}TSS_flanks200.bed', delimiter='\t')
    df_nonepi.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_id', 'gene_name']

    nonepi_dfs = []


    ## STEP 5: Get TSS exon flanks of hm-specific epigenes and non-epigenes
    for hm in list(hms):
        df = df_nonepi[df_nonepi.gene_name.isin(nonepigenes[hm])]
        
        df[ 'type'] = hm
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

    # get all TSS exons with DEU reported by the current tool
    get_TSS_exons(tool)

    if type == 'epi':
        get_epigenes_study(tool)
    else:
        get_nonepigenes(tool)