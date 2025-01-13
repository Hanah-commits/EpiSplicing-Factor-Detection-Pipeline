import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from sklearn.decomposition import PCA
import seaborn as sns


def plot_epigenes():
    epigenes = pd.read_csv('0_Files/Post-processing/Analyses/epigenes/Epigene_numbers.tsv', delimiter='\t').iloc[:, : 11]
    epigenes.set_index('HM', inplace=True)
    epigenes = epigenes.T
    epigenes = epigenes.replace('X', -1)
    epigenes = epigenes.astype(float)

    # COLOR MAP
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))


    print(epigenes)
    ax = epigenes.plot.bar(stacked=True, color=color_dict)

    # set title, axes etc
    hfont = {'fontname':'Calibri'}
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(10)
    plt.title('Number of Epigenes Across All Embryonic Cell Line Pairs', fontsize=12)
    plt.xticks(rotation=45)
    plt.ylabel('Number of Epigenes', fontsize=10)
    plt.legend(loc='upper right')
    plt.savefig('0_Files/Post-processing/Analyses/epigenes/Epigenes.png',bbox_inches='tight', dpi=300)
    # plt.close()


def plot_epiflanks():
    epiflanks = pd.read_csv('0_Files/Post-processing/Analyses/epigenes/Epiflanks_numbers.tsv', delimiter='\t')
    epiflanks.set_index('HM', inplace=True)
    # epiflanks = epiflanks.T
    epiflanks = epiflanks.replace('X', -1)
    epiflanks = epiflanks.astype(float)

    # COLOR MAP
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    # Create a list of colors for each bar segment
    colors_tp = [color_dict[idx] for idx in epiflanks.index]  # Colors for TP
    colors_fp = ["#e3e1d3"] * len(epiflanks)  # Grey for FP

    print(epiflanks)

    # Plot the stacked bar plot
    fig, ax = plt.subplots(figsize=(8, 5))

    epiflanks['Non-epispliced Exon Flanks'].plot(kind='bar', stacked=True, color=colors_fp, ax=ax, width=0.8, label='# Nonepispliced Exon Flanks') # Plot Nonepi flanks
    epiflanks['Epispliced Exon Flanks'].plot(kind='bar', stacked=True, color=colors_tp, ax=ax, width=0.8,  bottom=epiflanks['Non-epispliced Exon Flanks']) ## Plot Epiflanks on top

    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy() 
        ax.text(x+width/2, 
                y+height/2, 
                '{:.0f}'.format(height), 
                horizontalalignment='center', 
                verticalalignment='center')

    # Add legend and labels
    # set title, axes etc
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(10)
    plt.title('Number of Flanks of Epispliced and Non-epispliced Exons Available for All Histone Marks', fontsize=12)
    plt.xticks(rotation=45)
    plt.ylabel('Number of Exon Flanks', fontsize=10)
    plt.xlabel('')
    plt.savefig('0_Files/Post-processing/Analyses/epigenes/Epiflanks.png',bbox_inches='tight', dpi=300)


def PCA_plot(hm):

    os.makedirs("0_Files/Post-processing/Analyses/PCA/", exist_ok=True)  # Create the directory if it doesn't exist

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info
    # encode labels
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)

    sf = [val for val in features.columns if val != 'label']

    # keep only strong binding events
    sf_data = features[sf]
    sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val) 
    X, y = sf_data.values, features['label'].values

    # get custom color map
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    ## Plot # PCs
    threshold = 0.9
    pca = PCA(n_components=threshold)
    pca.fit(X)
    
    # calculate cumulative variance explained
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
    n_components = len(cumulative_variance)


    plt.rcParams["figure.figsize"] = (14, 8)
    fig, ax = plt.subplots()

    # x-axis: # components
    xi = np.arange(1, n_components + 1, step=1)

    # plot cumulative variance
    ax.plot(xi, cumulative_variance, marker='o', linestyle='-', color=color_dict[hm], lw=2, label='Cumulative Variance Explained')

    # plot threshold
    ax.axhline(y=threshold, color='red', linestyle='--', lw=1.5, label=f'{str(threshold)}% Variance Threshold')

    # edit plot
    ax.set_ylim(0.0, 1.1)
    ax.set_xlim(0, n_components + 1)
    ax.set_xlabel('Number of Components', fontsize=8)
    ax.set_ylabel('Cumulative Variance Explained', fontsize=8)
    ax.set_title(f'Number of Principal Components Needed to Explain Variance - {hm}', fontsize=10)
    ax.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.legend(loc='lower right', fontsize=12)

    # edit x-axis ticks
    ax.set_xticks(np.arange(1, n_components + 1, step=1))
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='both', labelsize=8)

    plt.savefig(f"0_Files/Post-processing/Analyses/PCA/PCA_{hm}_explained_variance.png", bbox_inches='tight', dpi=300)
    plt.close()
    

def heatmap_allRBPs(hm):

    os.makedirs("0_Files/Post-processing/Analyses/Heatmaps/Corr_Features/", exist_ok=True)  # Create the directory if it doesn't exist
    os.makedirs("0_Files/Post-processing/Analyses/Heatmaps/Binding_Scores/", exist_ok=True)

    print(f'\n\n{hm}')

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)    

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info
    

    sf = [val for val in features.columns if val != 'label']

    # keep only strong binding events
    features[sf] = features[sf]
    features[sf] = features[sf].applymap(lambda val: 0 if val < 2 else val)


    modes = ['corr', 'binding']

    for mode in modes:

        if mode == 'corr':
            # Plot corr coeffs
            corr = features[sf].corr().abs()

            plt.figure(figsize=(8, 6))
            ax = sns.heatmap(corr, annot=False, xticklabels=corr.columns, yticklabels=corr.index, cmap='coolwarm', linewidths=0.1)
            plt.title(f'Correlation of Predicted RBP Binding Scores - {hm}', fontsize=8)

            plt.xticks(np.arange(len(sf))+ 0.5, labels=corr.columns, fontsize=2, rotation=90)
            plt.yticks(np.arange(len(sf))+ 0.5, labels=corr.columns, fontsize=2, rotation=0)

            cbar = ax.collections[0].colorbar
            cbar.ax.tick_params(labelsize=5) 
            cbar.set_label("|Pearson's R|", fontsize=8)

            plt.savefig(f"0_Files/Post-processing/Analyses/Heatmaps/Corr_Features/corr_{hm}.png",dpi=300) #bbox_inches='tight')
            plt.close()

        else:

            hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
            color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

            label_colors = dict(zip(['epigene', 'non-epigene'], [color_dict[hm], "grey"]))
            features.set_index('label', inplace=True)    
            row_colors = features.index.map(label_colors)

            plt.figure(figsize=(8, 6))
            cluster = sns.clustermap(features[sf], annot=False, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='winter', cbar_pos=(0.9, .2, .03, .4))
            heatmap = cluster.ax_heatmap

            # Remove y-axis ticks and tick labels
            heatmap.set_yticks([])
            heatmap.set_yticklabels([])

            # Set x-axis tick labels
            heatmap.set_xticks(range(len(sf)))
            heatmap.set_xticklabels(sf, size=2, rotation=90)

            # add title and axes labels
            heatmap.set_title(f'Predicted Binding Scores of RBPs in Exon Flanks - {hm}', fontsize=10)
            heatmap.set_ylabel('Flanks of Epispliced and Non-epispliced Exons', fontsize ='8')
            cax = cluster.figure.axes[-1]
            cax.tick_params(labelsize=8)
            
            # plt.ylabel()
            plt.savefig(f"0_Files/Post-processing/Analyses/Heatmaps/Binding_Scores/bindingscores_{hm}.png", dpi=300) #bbox_inches='tight', 
            plt.close()


def ridgeplot_splice_site_scores(df, hm, ss_type):

    op_dir = '0_Files/Post-processing/Analyses'
    op_dir = os.path.join(op_dir, 'MAXENTSCAN')
    os.makedirs(op_dir, exist_ok=True)  # Create the directory if it doesn't exist

    print(df)

    # set color palette
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))
    custom_palette = { 'Non-epispliced Exon': '#e3e1d3'}
    custom_palette['Epispliced Exon'] = color_dict[hm]

    # change df to long format
    df_melted = df.melt(id_vars='Exon Class', value_vars=[f"{ss_type}'ss"], 
                        var_name='Variable', value_name='Value')

    # Create a FacetGrid for the ridge plot
    g = sns.FacetGrid(df_melted, row='Variable', hue='Exon Class', aspect=4, height=3, palette=custom_palette)

    # map kdeplot to each subplot
    g.map(sns.kdeplot, 'Value', fill=True, alpha=0.6, bw_adjust=0.6, levels = 10, common_norm = False)

    #  customize plot
    g.set_titles('')
    g.set_xlabels(f"{ss_type}' Splice Site Strength Scores")
    g.set_ylabels('Density')
    g.add_legend()
    g.despine(left=True)
    plt.subplots_adjust(hspace=0.5)
    g.fig.suptitle(f" Epispliced vs Non-epispliced Exons - {hm}", fontsize=10)
    g.fig.subplots_adjust(top=0.9)

    plt.savefig(f'{op_dir}/{hm}_ridgeplot_{ss_type}prime_splicesite.png', bbox_inches='tight', dpi=300)


def splice_site_strength_epi_nonepi():
    op_dir = '0_Files/Post-processing/Analyses/MAXENTSCAN/scores'

    epi_ss = pd.read_csv(f'{op_dir}/epi_exons_splicesite_scores.bed', delimiter='\t', header=None)
    nonepi_ss = pd.read_csv(f'{op_dir}/nonepi_exons_splicesite_scores.bed', delimiter='\t', header=None)

    epi_ss.columns = ['chr', 'start', 'stop', 'gene_name', 'type', 'strand', "5'ss", "3'ss"]
    nonepi_ss.columns = ['chr', 'start', 'stop', 'gene_name', 'type', 'strand', "5'ss", "3'ss"]

    epi_ss['Exon Class'] = 'Epispliced Exon'
    nonepi_ss['Exon Class'] = 'Non-epispliced Exon'

    all_ss = pd.concat([epi_ss, nonepi_ss], axis=0)

    # remove exons with both labels
    common_genes = list(set(epi_ss.gene_name.values.tolist()) & set(nonepi_ss.gene_name.values.tolist()))
    all_ss = all_ss[~((all_ss.gene_name.isin(common_genes)) & (all_ss['Exon Class'] == 'Non-epispliced Exon'))]
    all_ss = all_ss[~((all_ss.gene_name.isin(common_genes)) & (all_ss['Exon Class'] == 'Epispliced Exon'))]

    hms = [  "H3K27ac","H3K4me3","H3K9me3", "H3K36me3", "H3K27me3"]

    for hm in hms:
        print(hm) 
        features = all_ss[all_ss['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]


        for score in ["5", "3"]:
            ridgeplot_splice_site_scores(features, hm, score)


def heatmap_epiRBPs_correlation(hm):

    os.makedirs("0_Files/Post-processing/epiRBPS/corr_heatmaps/", exist_ok=True)  # Create the directory if it doesn't exist

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)    

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info

    rbps_file = open(f"0_Files/Post-processing/epiRBPS/epiRBPs_{hm}.txt", "r")
    epi_rbps = [rbp for rbp in rbps_file.read().split('\n') if rbp]
    features = features[epi_rbps]
    features = features.applymap(lambda val: 0 if val < 2 else val) # keep only strong binding events

    # Plot corr coeffs
    corr = features.corr()
    ax = sns.heatmap(corr, annot=False, cmap='coolwarm')
    plt.title(f'Correlation Coefficients of Episplicing RBPs Associated with {hm}', fontsize=8)
    plt.xticks(fontsize=5, rotation=90)
    plt.yticks(fontsize=5, rotation=0)
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=6)
    plt.savefig(f"0_Files/Post-processing/epiRBPS/corr_heatmaps/corr_{hm}.png",dpi=300) #bbox_inches='tight')
    plt.close()


def heatmap_imptRBPs_binding(hm):

    os.makedirs("0_Files/Post-processing/imptRBPS/bindinscores_heatmaps/", exist_ok=True)  # Create the directory if it doesn't exist

    #get epi and nonepiRBPs to plot
    features_to_plot = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPS/epiRBPs_{hm}.txt", "r")
    features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPS/nonepiRBPs_{hm}.txt", "r")
    features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)
    features.set_index('label', inplace=True)    


    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info
    features = features.applymap(lambda val: 0 if val < 2 else val) # keep only strong binding events

    # plot binding scores

    # custom colormap
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    label_colors = dict(zip(['epigene', 'non-epigene'], [color_dict[hm], "grey"]))
    row_colors = features.index.map(label_colors)

    cluster = sns.clustermap(features[features_to_plot], annot=False, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='winter',  cbar_pos=(0.9, .2, .03, .4))
    heatmap = cluster.ax_heatmap

    # remove y-axis ticks and tick labels
    heatmap.set_yticks([])
    heatmap.set_yticklabels([])

    # set x-axis tick labels
    heatmap.set_xticks(range(len(features_to_plot)))
    heatmap.set_xticklabels(features_to_plot, size=8, rotation=45)

    # add title and axes labels
    heatmap.set_title(f'{hm} : Predicted Binding Scores of Episplicing and Non-episplicing RBPs', fontsize=10)
    heatmap.set_ylabel(f'Flanks of Epispliced and Non-epispliced Exons', fontsize =10, labelpad=-505)
    cax = cluster.figure.axes[-1]
    cax.tick_params(labelsize=8)

    colorbar = cluster.cax
    colorbar.set_ylabel("Predicted Binding Scores", rotation=90, labelpad=10)
    
    plt.savefig(f"0_Files/Post-processing/imptRBPS/bindinscores_heatmaps/bindingscores_{hm}.png", bbox_inches='tight', dpi=300)
    plt.close()


def prep_log2_norm_counts():
    # read normalised counts tsv
    counts = pd.read_csv('0_Files/Post-processing/Analyses/expression/counts/rpkm_values_rbps.tsv', delimiter='\t')

    # remove unwanted prefix from column names
    counts.columns = counts.columns.str.replace(r'validation.', '')

    #remove unwanted suffx
    celllines = set(col.split('_')[0] for col in counts.columns if col != 'gene')

    # Rename columns: H1_ENCHAAQAH.bam -> H1_rep{i}
    rename_mapping = {}
    for prefix in celllines:
        replicate_counter = 1
        for col in counts.columns:
            if col.startswith(prefix):
                clean_prefix = prefix.replace("dermalcell", "") #ectodermalcell -> ecto
                clean_prefix = prefix.replace("cell", "") #neuronalcell -> neuro
                rename_mapping[col] = f"{clean_prefix}_rep{replicate_counter}"
                replicate_counter += 1

    counts.rename(columns=rename_mapping, inplace=True)


    # replace rbp names with ids used in this study
    replacement_dict = {
        'CELF4':'BRUNOL4',
        'CELF5':'BRUNOL5',
        'CELF6':'BRUNOL6',
        'ELAVL1': 'HuR',
        'HNRNPLL':'HNRPLL'
    }

    counts['gene'] = counts['gene'].replace(replacement_dict)


    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3"]    
    for hm in hms:
        print(hm)
        sfs = []
        with open(f'0_Files/Post-processing/epiRBPS/epiRBPs_{hm}.txt', 'r') as file:
            sfs.extend([line.strip() for line in file.readlines()])      


if __name__ == "__main__":

    plot_epigenes() # Figure 1
    plot_epiflanks() # Figure 2
    splice_site_strength_epi_nonepi() # Figure 14 (Suppl)
    prep_log2_norm_counts() # Fig 19 (Suppl)

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']

    for hm in hms:
        heatmap_allRBPs(hm) # Figure 11 (Supp)
        PCA_plot(hm) # Figure 12 (Supp)
        heatmap_epiRBPs_correlation(hm) # Fig 16 (Supp)
        heatmap_imptRBPs_binding(hm) # Fig 17 (Supp)


    