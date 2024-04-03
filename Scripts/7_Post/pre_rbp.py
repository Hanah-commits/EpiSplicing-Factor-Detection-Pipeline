import pandas as pd
import os
import json
from argparse import ArgumentParser

def flanked_exons():      

    with open('paths.json') as f:
        d = json.load(f)

    fasta = d['Reference fasta']

    op_dir = "0_Files/Post-processing"
    input_files = []

    for i in range(2):
        name = 'epi' if i == 0 else 'nonepi'

        # extend exon boundary 
        os.system(f"bedtools slop -i {op_dir}/{name}_exons.bed -g  {fasta}.fai  -b 200 | sort | uniq > {op_dir}/{name}_flanked_exons.bed" )

        flanked_exons = pd.read_csv(f'{op_dir}/{name}_flanked_exons.bed', delimiter='\t', header=None)
        flanked_exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']
        flanked_exons['flanks'] = flanked_exons[['exon_start', 'exon_end']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)

        
        # Get input sequences for RBPmap
        n = 0
        i = 1
        while n < len(flanked_exons):
                file_name = f'../RBPmap/rbp_input_{name}{i}.csv'
                input_files.append(file_name)
                if n+5000 <= len(flanked_exons):
                        flanked_exons[['chr', 'flanks', 'strand']].iloc[n:n+5000].to_csv(file_name, index=False, sep=':', header=False)
                        n += 5000
                else:
                        flanked_exons[['chr', 'flanks', 'strand']].iloc[n:len(flanked_exons)+1].to_csv(file_name, index=False, sep=':', header=False)
                        break
                i += 1


        # # get flanks for feature matrix preparation step
        flanked_exons[['gene_name', 'flanks']].to_csv(f'{op_dir}/query_flanks_{name}.csv', sep='\t', index=False)

    input_files = [os.getcwd() + '/' + file for file in input_files]
    with open('../RBPmap/input.txt', 'w') as f:
        for item in input_files:
            f.write("%s\n" % item)

def exon_body():  

    op_dir = "0_Files/Post-processing"
    input_files = []
    
    for i in range(2):
        name = 'epi' if i == 0 else 'nonepi'

        exons = pd.read_csv(f'{op_dir}/{name}_exons.bed', delimiter='\t', header=None)
        exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']
        exons['flanks'] = exons[['exon_start', 'exon_end']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
  
        # Get input sequences for RBPmap
        n = 0
        i = 1
        while n < len(exons):
                file_name = f'../RBPmap/rbp_input_{name}{i}.csv'
                input_files.append(file_name)
                if n+5000 <= len(exons):
                        exons[['chr', 'flanks', 'strand']].iloc[n:n+5000].to_csv(file_name, index=False, sep=':', header=False)
                        n += 5000
                else:
                        exons[['chr', 'flanks', 'strand']].iloc[n:len(exons)+1].to_csv(file_name, index=False, sep=':', header=False)
                        break
                i += 1


        # # get flanks for feature matrix preparation step
        exons[['gene_name', 'flanks']].to_csv(f'{op_dir}/query_flanks_{name}.csv', sep='\t', index=False)

    input_files = [os.getcwd() + '/' + file for file in input_files]
    with open('../RBPmap/input.txt', 'w') as f:
        for item in input_files:
            f.write("%s\n" % item)


def exon_flanks():

    op_dir = "0_Files/Post-processing"
    input_files = []

    for i in range(2):
        
        name = 'epi' if i == 0 else 'nonepi'

        flanks = pd.read_csv(f'{op_dir}/{name}_flanks.bed', delimiter='\t', header=None)
        flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']
        flanks['flanks'] = flanks[['flank_start', 'flank_end']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        
        # Get input sequences for RBPmap
        n = 0
        i = 1
        while n < len(flanks):
                file_name = f'../RBPmap/rbp_input_{name}{i}.csv'
                input_files.append(file_name)
                if n+5000 <= len(flanks):
                        flanks[['chr', 'flanks', 'strand']].iloc[n:n+5000].to_csv(file_name, index=False, sep=':', header=False)
                        n += 5000
                else:
                        flanks[['chr', 'flanks', 'strand']].iloc[n:len(flanks)+1].to_csv(file_name, index=False, sep=':', header=False)
                        break
                i += 1


        # # get flanks for feature matrix preparation step
        flanks[['gene_name', 'flanks']].to_csv(f'{op_dir}/query_flanks_{name}.csv', sep='\t', index=False)

    input_files = [os.getcwd() + '/' + file for file in input_files]
    with open('../RBPmap/input.txt', 'w') as f:
        for item in input_files:
            f.write("%s\n" % item)


if __name__ == "__main__":

    p = ArgumentParser()
    p.add_argument("mode", help='flanks, exon, flanked')
    args = p.parse_args()
    mode = args.mode
     
    if mode == 'flanks':
        exon_flanks()
    elif mode == 'exon':
        exon_body()
    elif mode == 'flanked':
        flanked_exons()
    else:
        raise ValueError('Wrong value for paramter "mode".')