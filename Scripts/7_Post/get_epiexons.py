import pandas as pd
import json
import os

def bedtools_input(hm, all_flanks):

    file = '0_Files/dPSI_Mval_epi_' + hm + '.csv'
    hm_flanks = pd.read_csv(file, delimiter='\t')
    epi_flanks = all_flanks[all_flanks['flanks'].isin(list(set(hm_flanks['flanks'].values)))]
    epi_flanks['score'] = '.'
    epi_flanks['feature'] = 'flanks'

    #bedtools input
    epi_flanks[['seqid', 'start', 'stop', 'feature', 'score', 'strand']].to_csv('0_Files/flanks_'+ hm + '.bed', sep='\t', index=False, header=False)
            

def run_bedtools(hm):

    file = '0_Files/flanks_'+ hm + '.bed'
    os.system('bedtools intersect -loj -s -a 0_Files/exon_coords.bed -b ' + file + '| sort | uniq > 0_Files/epiexons_' + hm + '.bed')


def post_bedtools(hm):

    exons = pd.read_csv('0_Files/epiexons_' + hm + '.bed', delimiter='\t', header=None)
    # # drop feature, strand, score etc.
    # exons.drop([3, 4, 5, 9, 10, 11], axis=1, inplace=True)
    # # assign 0 to exons that have no flanks
    # exons.replace([-1, '.'], [0, 0], inplace=True)
    # exons = exons.set_axis(['seqid', 'exon_start', 'exon_stop', 'chr', 'flank_start', 'flank_stop'],
    #                axis=1) #=True)
    
    # # keep exons that have flank annotation
    # exons = exons[exons['flank_start'] !=0]

    # #keep seq id and co-ordinates of exons
    # exons = exons[['seqid', 'exon_start', 'exon_stop']]
    # exons.drop_duplicates(inplace=True)

    exons.to_csv('0_Files/epiexons_' + hm + '.bed', sep='\t', header=False, index=None)


if __name__ == "__main__":

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    
    flanks = pd.read_csv('0_Files/all_flanks.csv', delimiter='\t')

    for hm in hms[:1]:
        # bedtools_input(hm, flanks)
        # run_bedtools(hm)
        post_bedtools(hm)