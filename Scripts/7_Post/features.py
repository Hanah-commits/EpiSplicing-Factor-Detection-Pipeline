import pandas as pd
from argparse import ArgumentParser

def feature_matrix(filename1, filename2, mode):

    exons = pd.read_csv(filename1, delimiter='\t', header=None)

    if mode == 'flanks':
        exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name', 'type']
    elif mode == 'exon' or mode == 'flanked':
        exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']

    rbp = pd.read_csv(filename2, delimiter=',')

    features = pd.concat([exons, rbp], axis=1)

    name = ''
    if i == 0:
        name = 'epi'
    else:
        name = 'nonepi'

    # # FILTER 1: If RBP motif has 1+ PSSMs, keep only one
    features = features.loc[:, ~features.columns.duplicated()]
    features.to_csv('0_Files/Post-processing/features_' + name + '.csv', sep='\t', index=False)


if __name__ == "__main__":

    p = ArgumentParser()
    p.add_argument("mode", help='flanks, exon, flanked')
    args = p.parse_args()
    mode = args.mode
     
    if mode == 'flanks':
        exons_files = ['0_Files/Post-processing/epi_flanks.bed', '0_Files/Post-processing/nonepi_flanks.bed']
    elif mode == 'exon' or mode == 'flanked':
        exons_files = ['0_Files/Post-processing/epi_exons.bed', '0_Files/Post-processing/nonepi_exons.bed']
    else:
        raise ValueError('Wrong value for paramter "mode".')

    Zscore_files = ['0_Files/Post-processing/FilteredZscores_epi.csv', '0_Files/Post-processing/FilteredZscores_nonepi.csv']

    for i in range(2):
        feature_matrix(exons_files[i], Zscore_files[i], mode)

