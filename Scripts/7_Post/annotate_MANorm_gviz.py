import os
import sys
import json
import pandas as pd
from pathlib import Path
import numpy as np
from argparse import ArgumentParser
from statsmodels.stats.multitest import multipletests


def manorm_to_bigwig(): 

    # Get the process name, use it in the output directory
    def get_argument_parser():
        p = ArgumentParser()
        p.add_argument("output_dir")
        p.add_argument("--process", "-p",
            help="The name of the process")
        return p


    p = get_argument_parser()
    args = p.parse_args()
    proc = args.process


    with open('paths.json') as f:
            data = json.load(f)
    d = data[proc]

    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]
    hms = d["Histone modifications"]
    fasta = d['Reference fasta']
    ref_genome= fasta+".fai"

    prefix = args.output_dir + 'MANorm/'

    print('Annotating candidate exons with HM peaks from MANorm \n') #log

    for hm in ['H3K36me3']:
        print(hm)

        # filtered_peaks = f'{prefix}{hm}_{tissue1}_vs_{tissue2}_DHM.bed'
        all_peaks = f'{prefix}{hm}_{tissue1}_peak_vs_{hm}_{tissue2}_peak_all_MAvalues.xls'
        prefix = '0_Files/gviz'

        # Filter unique peaks using adjusted pvalue
        unique_peaks = pd.read_csv(all_peaks, delimiter='\t')
        unique_peaks['adj_pval'] = multipletests(unique_peaks['P_value'], method='fdr_bh')[1]
        unique_peaks = unique_peaks[unique_peaks['adj_pval'] < 0.05]
        output_file1 = f'{hm}_{tissue1}.bedgraph'
        output_file2 = f'{hm}_{tissue2}.bedgraph'
        unique_peaks[['chr', 'start', 'end', f'normalized_read_density_in_{hm}_{tissue1}_peak']].to_csv(f'{prefix}/{output_file1}', sep='\t',header=False, index=False)
        unique_peaks[['chr', 'start', 'end', f'normalized_read_density_in_{hm}_{tissue2}_peak']].to_csv(f'{prefix}/{output_file2}', sep='\t',header=False, index=False)
        
        # Sort peaks
        os.system(f'sort -k1,1 -k2,2n {prefix}/{output_file1} > {prefix}/tmp.bedgraph && mv {prefix}/tmp.bedgraph {prefix}/{output_file1}')
        os.system(f'sort -k1,1 -k2,2n {prefix}/{output_file2} > {prefix}/tmp.bedgraph && mv {prefix}/tmp.bedgraph {prefix}/{output_file2}')

        # Convert to bigwig
        os.system(f'bedGraphToBigWig {prefix}/{output_file1} HelperFunctions/hg38.chrom.sizes {prefix}/{hm}_{tissue1}.bw')
        os.system(f'bedGraphToBigWig {prefix}/{output_file2} HelperFunctions/hg38.chrom.sizes {prefix}/{hm}_{tissue2}.bw')

        os.system(f'rm {prefix}/*.bedgraph')


def eclip_peaks_to_bigwig(): 


    tissue1 = "HepG2"
    tissue2 = "K562"

    for rbp in ['PTBP1', 'U2AF1']:
        print(rbp)
        for tissue in [tissue1, tissue2]:
            peaks_file = f'~/data/validation/eclip/{rbp}_{tissue}.bed'
            prefix = '0_Files/gviz'

            output_file = f'{rbp}_{tissue}.bedgraph'
    
            # convert to bedgraph
            os.system(f'cut -f 1,2,3,7 {peaks_file} > {prefix}/{output_file}')

            # sort peaks
            os.system(f'sort -k1,1 -k2,2n {prefix}/{output_file} > {prefix}/tmp.bedgraph && mv {prefix}/tmp.bedgraph {prefix}/{output_file}')

            # Convert to bigwig
            os.system(f'bedGraphToBigWig {prefix}/{output_file} HelperFunctions/hg38.chrom.sizes {prefix}/{rbp}_{tissue}.bw')


def manorm_summit_to_bigwig(): 

    # Get the process name, use it in the output directory
    def get_argument_parser():
        p = ArgumentParser()
        p.add_argument("output_dir")
        p.add_argument("--process", "-p",
            help="The name of the process")
        return p


    p = get_argument_parser()
    args = p.parse_args()
    proc = args.process

    with open('paths.json') as f:
            data = json.load(f)
    d = data[proc]

    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]
    hms = d["Histone modifications"]
    fasta = d['Reference fasta']
    ref_genome= fasta+".fai"

    prefix = args.output_dir + 'MANorm/'

    print('Annotating candidate exons with HM peaks from MANorm \n') #log

    for hm in ['H3K36me3']:
        print(hm)

        all_peaks = f'{prefix}{hm}_{tissue1}_peak_vs_{hm}_{tissue2}_peak_all_MAvalues.xls'

        prefix = '0_Files/gviz'
        Path(prefix).mkdir(parents=True, exist_ok=True)

        # Filter unique peaks using adjusted pvalue
        unique_peaks = pd.read_csv(all_peaks, delimiter='\t')
        unique_peaks['summit+1'] = unique_peaks['summit']+1
        print(unique_peaks[unique_peaks['Peak_Group'] != 'H3K4me3_K562_peak_unique'])
        unique_peaks['adj_pval'] = multipletests(unique_peaks['P_value'], method='fdr_bh')[1]
        unique_peaks = unique_peaks[unique_peaks['adj_pval'] < 0.05]

        print(unique_peaks[unique_peaks['Peak_Group'] != 'H3K4me3_K562_peak_unique'])
        output_file1 = f'{hm}_{tissue1}.bedgraph'
        output_file2 = f'{hm}_{tissue2}.bedgraph'
        # unique_peaks[['chr', 'summit', 'summit+1', f'normalized_read_density_in_{hm}_{tissue1}_peak']].to_csv(f'{prefix}/{output_file1}', sep='\t',header=False, index=False)
        # unique_peaks[['chr', 'summit', 'summit+1', f'normalized_read_density_in_{hm}_{tissue2}_peak']].to_csv(f'{prefix}/{output_file2}', sep='\t',header=False, index=False)
        
        # # Sort peaks
        # os.system(f'sort -k1,1 -k2,2n {prefix}/{output_file1} > {prefix}/tmp.bedgraph && mv {prefix}/tmp.bedgraph {prefix}/{output_file1}')
        # os.system(f'sort -k1,1 -k2,2n {prefix}/{output_file2} > {prefix}/tmp.bedgraph && mv {prefix}/tmp.bedgraph {prefix}/{output_file2}')

        # # Convert to bigwig
        # os.system(f'bedGraphToBigWig {prefix}/{output_file1} HelperFunctions/hg38.chrom.sizes {prefix}/{hm}_{tissue1}.bw')
        # os.system(f'bedGraphToBigWig {prefix}/{output_file2} HelperFunctions/hg38.chrom.sizes {prefix}/{hm}_{tissue2}.bw')

        # os.system(f'rm {prefix}/*.bedgraph')


if __name__ == "__main__":

    # Use process_id from config file as argument
    # Eg. if pr11 has K562-HepG2 pair: python 7_Post/annotate_MANorm_gviz.py pr11
    manorm_peaks_to_bigwig()
    eclip_peaks_to_bigwig()
    # manorm_summit_to_bigwig()