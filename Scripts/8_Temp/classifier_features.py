import pandas as pd
import json

with open('paths.json') as f:
    data = json.load(f)

hms = data["Histone modifications"]

prefix = '0_Files/Post-processing'
unscaled = ['features_epi.csv', 'features_nonepi.csv', 'all_features.csv']


epi_features = pd.read_csv('0_Files/Post-processing/features_epi.csv', delimiter='\t')
nonepi_features = pd.read_csv('0_Files/Post-processing/features_nonepi.csv', delimiter='\t')

# keep nonepi flanks for hms available in current study
nonepi_features = nonepi_features[nonepi_features['type'].apply(lambda x: any(item in hms for item in x.split(',')))]
nonepi_features.to_csv('0_Files/Post-processing/features_nonepi.csv', sep='\t', index=False)

epi_features['label'] = 'epigene'
nonepi_features['label'] = 'non-epigene'

all_features = pd.concat([epi_features, nonepi_features], axis=0)

# remove genes with both labels
common_genes = list(set(epi_features.gene_name.values.tolist()) & set(nonepi_features.gene_name.values.tolist()))
all_features = all_features[~((all_features.gene_name.isin(common_genes)) & (all_features.label == 'non-epigene'))]
all_features = all_features[~((all_features.gene_name.isin(common_genes)) & (all_features.label == 'epigene'))]

all_features.drop(['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name'], axis=1, inplace=True)


col = all_features.pop("label")
all_features.insert(0, col.name, col)

all_features.to_csv('0_Files/Post-processing/features_all.csv', sep='\t', index=False)

