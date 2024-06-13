import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def binding_sites():
        
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
        if len(sfs) > 0:
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

        if len(sfs) > 0:
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



def expression_levels():
        
    op_dir = '0_Files/RBPs'

    # TODO: log2 transform counts
    counts = pd.read_csv(f'{op_dir}/normalized_RBP_counts.tsv', delimiter='\t')


    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K4me1", "H3K36me3"]
    
    for hm in hms:
        print(hm)

        for type in ['epi', 'nonepi']:
        
            # ## Get Epi/NonEpi proteins
            sfs = []
            with open(f'0_Files/Post-processing/enriched_{type}_{hm}.txt', 'r') as file:
                sfs.extend([line.strip() for line in file.readlines()])

            sfs = list(set(sfs))
            counts_hm = counts[counts.gene.isin(sfs)].sort_values('gene') 
            # # #TODO Not all RBPs in norm counts
            # # print(len(sfs))
            # # print(len(counts_hm))
            sfs_index = counts_hm['gene'].values.tolist()
            counts_hm = counts_hm.set_index('gene')
            

            if len(sfs_index) > 0:
                plt.close('all')

                # Calculate the range of the data to ensure the colormap is centered on 0
                vmin = counts_hm.loc[sfs_index].values.min()
                vmax = counts_hm.loc[sfs_index].values.max()
                max_abs = max(abs(vmin), abs(vmax))

                # Create the heatmap
                fig, heatmap = plt.subplots(figsize=(10, 10))
                sns.heatmap(
                    counts_hm.loc[sfs_index],
                    annot=False,
                    linewidths=0.3,
                    cmap='coolwarm',
                    center=0,
                    vmin=-max_abs,
                    vmax=max_abs,
                    ax=heatmap,
                    cbar_kws={'shrink': 0.5}  # Adjusts the colorbar size
                )

                # Customize the colorbar
                cbar = heatmap.collections[0].colorbar
                cbar.set_ticks(np.linspace(-max_abs, max_abs, num=5))  # Customize the number of ticks
                cbar.set_ticklabels([f'{i:.2f}' for i in np.linspace(-max_abs, max_abs, num=5)])  # Customize tick labels format
                cbar.ax.yaxis.set_ticks_position('right')  # Position the colorbar ticks to the left

                # Add a title to the colorbar
                cbar.set_label('log2(RPKM+1)')

                # Set y-axis tick labels
                heatmap.set_yticks(range(len(sfs_index)))
                heatmap.set_yticklabels(sfs_index, size=6)
                heatmap.yaxis.tick_left()

                # Set x-axis tick labels
                heatmap.set_xticks(range(len(counts_hm.columns)))
                heatmap.set_xticklabels([col.split('bam')[0] for col in counts_hm.columns.tolist()], size=7, rotation=90)

                # add title and axes labels
                substr = 'Epi' if type == 'epi' else 'Non-epi'
                heatmap.set_title(f'{hm} : {substr}splicing RBPs')
                heatmap.set_ylabel('')

                # adjust layout for better appearance
                plt.tight_layout()

                plt.savefig(f'0_Files/Post-processing/{hm}/{type}/log2_RPKM.png')
        

binding_sites()
expression_levels()