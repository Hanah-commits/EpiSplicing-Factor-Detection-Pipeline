import pandas as pd

prefix = '0_Files/Post-processing'
unscaled = ['features_epi.csv', 'features_nonepi.csv', 'all_features.csv']


epi_features = pd.read_csv('0_Files/Post-processing/features_epi.csv', delimiter='\t')
nonepi_features = pd.read_csv('0_Files/Post-processing/features_nonepi.csv', delimiter='\t')

epi_features['label'] = 'epigene'
nonepi_features['label'] = 'non-epigene'

all_features = pd.concat([epi_features, nonepi_features], axis=0)
all_features.drop(['chr', 'exon_start', 'exon_end', 'feature', 'type', 'strand', 'gene_name'], axis=1, inplace=True)
col = all_features.pop("label")
all_features.insert(0, col.name, col)

all_features.to_csv('0_Files/Post-processing/features_all.csv', sep='\t', index=False)