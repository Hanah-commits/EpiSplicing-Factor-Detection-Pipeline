import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

def violinplot(data, hm, type):

    op_dir = '0_Files/Post-processing'

    # Create a directory specific to the hm value
    hm_dir = os.path.join(op_dir, hm, type)
    os.makedirs(hm_dir, exist_ok=True)  # Create the directory if it doesn't exist

    type = type.upper()
    df = pd.DataFrame(data)
    cols = [col for col in df.columns if col != 'label']
    df['label'] = df['label'].replace('epigene', 'Epi Exon Flanks')
    df['label'] = df['label'].replace('non-epigene', 'NonEpi Exon Flanks')
    df = df.rename(columns={'label': 'Sequence Class'})


    # set color palette
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["darkorchid", "red", "magenta", "darkorange", "green", "blue"]))
    custom_palette = { 'NonEpi Exon Flanks': 'lightgray'}
    custom_palette['Epi Exon Flanks'] = color_dict[hm]

    # Melt the DataFrame to long format
    df_melted = df.melt(id_vars='Sequence Class', value_vars=cols, 
                        var_name='Variable', value_name='Value')

    # Create the violin plot
    plt.figure(figsize=(12, 6))

    sns.violinplot(x='Variable', y='Value', hue='Sequence Class', data=df_melted, scale='area', split=False, palette=custom_palette) #, cut=0)


    plt.title(f'Binding Affinities of {type} RBPs - {hm}')
    plt.ylabel('Binding Scores')
    plt.xlabel('')
    plt.tick_params(axis='x', labelsize=7, labelrotation=90)
    plt.savefig(f'{hm_dir}/{hm}_violinplot_{type}.png')


def rideplot(data, hm, type):

    op_dir = '0_Files/Post-processing'

    type = type.upper()
    df = pd.DataFrame(data)
    cols = [col for col in df.columns if col != 'label']
    df['label'] = df['label'].replace('epigene', 'Epi Exon Flanks')
    df['label'] = df['label'].replace('non-epigene', 'NonEpi Exon Flanks')
    df = df.rename(columns={'label': 'Sequence Class'})


    # set color palette
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))
    custom_palette = { 'NonEpi Exon Flanks': 'lightgray'}
    custom_palette['Epi Exon Flanks'] = color_dict[hm]

    # Melt the DataFrame to long format
    df_melted = df.melt(id_vars='Sequence Class', value_vars=cols, 
                        var_name='Variable', value_name='Value')

    # Create a FacetGrid for the ridge plot
    g = sns.FacetGrid(df_melted, row='Variable', hue='Sequence Class', aspect=4, height=1.5, palette=custom_palette)

    # Map the kdeplot to each subplot
    g.map(sns.kdeplot, 'Value', fill=True, alpha=0.6, bw_adjust=0.6)

    # Add titles and adjust layout
    g.set_titles('')
    g.set_xlabels('')
    g.set_ylabels('')
    g.despine(left=True)
    plt.subplots_adjust(hspace=0.5)
    g.fig.suptitle(f'Binding Affinities of {type} RBPs - {hm}', fontsize=16)
    g.fig.subplots_adjust(top=0.9)  # Adjust title position

    plt.savefig(f'{op_dir}/{hm}_ridgeplot_{type}.png')


def rideplot_indiv(data, hm, type):

    op_dir = '0_Files/Post-processing'

    # Create a directory specific to the hm value
    hm_dir = os.path.join(op_dir, hm, type)
    os.makedirs(hm_dir, exist_ok=True)  # Create the directory if it doesn't exist

    type = type.upper()
    df = pd.DataFrame(data)
    cols = [col for col in df.columns if col != 'label']
    df['label'] = df['label'].replace('epigene', 'Epi Exon Flanks')
    df['label'] = df['label'].replace('non-epigene', 'NonEpi Exon Flanks')
    df = df.rename(columns={'label': 'Sequence Class'})


    # set color palette
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))
    custom_palette = { 'Epi Exon Flanks': color_dict[hm]}
    custom_palette['NonEpi Exon Flanks'] = 'gray'

    # Melt the DataFrame to long format
    df_melted = df.melt(id_vars='Sequence Class', value_vars=cols, 
                        var_name='Variable', value_name='Value')

    # Create a ridge plot for each variable
    for variable in cols:
            plt.figure(figsize=(8, 4))
            sns.kdeplot(data=df_melted[df_melted['Variable'] == variable], x='Value', hue='Sequence Class',
                        fill=True, alpha=0.6, palette=custom_palette, 
                        hue_order=['NonEpi Exon Flanks', 'Epi Exon Flanks'], 
                        common_norm=False,
                        bw_adjust=0.4)
 

            plt.title(f'Ridge Plot for {variable}')
            plt.xlabel('Binding Affinity (ZScore)')
            plt.ylabel('Density')

            plt.tight_layout()
            plt.savefig(f'{hm_dir}/{hm}_ridgeplot_{variable}_{type}.png')
            plt.close()


def enriched_epi_nonepi():
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

    hms = [  "H3K27ac","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1", "H3K27me3"]

    for hm in hms:
        print(hm) 
        features = all_features[all_features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
        # features = features.set_index('label')
        

        ## EPI-ENRICHED RBPS
        with open(f'{op_dir}/enriched_epi_{hm}.txt', 'r') as file:
            sfs = [line.strip() for line in file.readlines()]

        # features.loc[:, sfs] = features.loc[:, sfs].applymap(lambda val: 0 if val < 2 else val)
        cols = sfs + ['label']
        df = features[cols]
        if len(sfs) > 0:
            df.loc[:, sfs] = df.loc[:, sfs].applymap(lambda val: 0 if val < 2 else val)
            rideplot_indiv(df, hm, 'epi')


        ## NONEPI ENRICHED RBPS
        with open(f'{op_dir}/enriched_nonepi_{hm}.txt', 'r') as file:
            sfs = [line.strip() for line in file.readlines()]

        cols = sfs + ['label']
        df = features[cols]
        if len(sfs) >0:
            df.loc[:, sfs] = df.loc[:, sfs].applymap(lambda val: 0 if val < 2 else val)
            rideplot_indiv(df, hm, 'nonepi')

enriched_epi_nonepi()