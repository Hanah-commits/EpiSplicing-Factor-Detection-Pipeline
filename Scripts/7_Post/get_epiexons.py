import pandas as pd
import json
import os

def bedtools_input(hm, all_flanks):

    file = '0_Files/dPSI_Mval_epi_' + hm + '.csv'
    hm_flanks = pd.read_csv(file, delimiter='\t')
    epi_flanks = all_flanks[all_flanks['flanks'].isin(list(set(hm_flanks['flanks'].values)))]

    #bedtools input
    epi_flanks[['seqid', 'start', 'stop', 'strand']].to_csv('0_Files/flanks_'+ hm + '.bed', sep='\t', index=False, header=False)
            

def run_bedtools(hm):

    file = '0_Files/flanks_'+ hm + '.bed'
    os.system('bedtools intersect -loj -s -a 0_Files/exon_coords.bed -b ' + file + '| sort | uniq > 0_Files/epiexons_' + hm + '.bed')


def post_bedtools(hm):
    """
    #TODO
    """


if __name__ == "__main__":

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    
    flanks = pd.read_csv('0_Files/all_flanks.csv', delimiter='\t')

    for hm in hms[:1]:
        bedtools_input(hm, flanks)
        run_bedtools(hm)