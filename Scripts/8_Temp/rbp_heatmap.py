import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

op_dir = '0_Files/Post-processing'

epi_features = pd.read_csv(f'{op_dir}/features_epi.csv', delimiter='\t')
nonepi_features = pd.read_csv(f'{op_dir}/features_nonepi.csv', delimiter='\t')

epi_features['label'] = 'epigene'
nonepi_features['label'] = 'non-epigene'

all_features = pd.concat([epi_features, nonepi_features], axis=0)

# remove genes with both labels
common_genes = list(set(epi_features.gene_name.values.tolist()) & set(nonepi_features.gene_name.values.tolist()))
all_features = all_features[~((all_features.gene_name.isin(common_genes)) & (all_features.label == 'non-epigene'))]
all_features = all_features[~((all_features.gene_name.isin(common_genes)) & (all_features.label == 'epigene'))]

hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K4me1", "H3K36me3"]
color_dict = dict(zip(hms,["purple", "red", "magenta", "orange", "green", "blue"]))

for hm in hms:
    print(hm)
    features = all_features[all_features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.set_index('label')
    

    ## EPI-ENRICHED RBPS
    with open(f'{op_dir}/enriched_epi_{hm}.txt', 'r') as file:
        sfs = [line.strip() for line in file.readlines()]

#     # sfs = [ val for val in list(features.columns) if val not in ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name', 'type']]

    # features.loc[:, sfs] = features.loc[:, sfs].applymap(lambda val: 0 if val < 3 else val)

    plt.close('all')
    plt.figure(figsize=(10, 8))

    label_colors = dict(zip(['epigene', 'non-epigene'], [color_dict[hm], "grey"]))
    row_colors = features.index.map(label_colors)
    cluster = sns.clustermap(features[sfs], annot=False, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='viridis',  cbar_pos=(0.9, .2, .03, .4))
    heatmap = cluster.ax_heatmap
    cbar = heatmap.collections[0].colorbar # custom y tick colorbar

    # Remove y-axis ticks and tick labels
    heatmap.set_yticks([])
    heatmap.set_yticklabels([])

    # Set x-axis tick labels
    heatmap.set_xticks(range(len(sfs)))
    heatmap.set_xticklabels(sfs, size=6, rotation=90)

    # add title and axes labels
    heatmap.set_title(f'{hm} : Epi-enriched RBPs')
    heatmap.set_ylabel('')

    # Position the legend next to the plot
    plt.savefig(f'{op_dir}/{hm}_epi.png')


    ## NONEPI ENRICHED RBPS
    with open(f'{op_dir}/enriched_nonepi_{hm}.txt', 'r') as file:
        sfs = [line.strip() for line in file.readlines()]

    plt.close('all')
    plt.figure(figsize=(10, 8))

    label_colors = dict(zip(['epigene', 'non-epigene'], [color_dict[hm], "grey"]))
    row_colors = features.index.map(label_colors)
    cluster = sns.clustermap(features[sfs], annot=False, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='viridis',  cbar_pos=(0.9, .2, .03, .4))
    heatmap = cluster.ax_heatmap
    cbar = heatmap.collections[0].colorbar # custom y tick colorbar

    # Remove y-axis ticks and tick labels
    heatmap.set_yticks([])
    heatmap.set_yticklabels([])

    # Set x-axis tick labels
    heatmap.set_xticks(range(len(sfs)))
    heatmap.set_xticklabels(sfs, size=6, rotation=90)

    # add title and axes labels
    heatmap.set_title(f'{hm} : Nonepi-enriched RBPs')
    heatmap.set_ylabel('')

    # Position the legend next to the plot
    plt.savefig(f'{op_dir}/{hm}_nonepi.png')