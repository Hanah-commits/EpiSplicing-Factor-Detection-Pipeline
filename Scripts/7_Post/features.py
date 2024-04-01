import pandas as pd

def feature_matrix(filename1, filename2):

    exons = pd.read_csv(filename1, delimiter='\t', header=None)
    exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'type', 'strand', 'gene_name']
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

    exons_files = ['0_Files/Post-processing/epi_exons.bed', '0_Files/Post-processing/nonepi_exons.bed']
    Zscore_files = ['0_Files/Post-processing/FilteredZscores_epi.csv', '0_Files/Post-processing/FilteredZscores_nonepi.csv']
    query_files = ['0_Files/Post-processing/query_flanks_epi.csv', '0_Files/Post-processing/query_flanks_nonepi.csv']

    for i in range(len(query_files)):
        feature_matrix(exons_files[i], Zscore_files[i])
