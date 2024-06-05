import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def violinplot(data, hm, type):

    op_dir = '0_Files/Post-processing'

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
    plt.savefig(f'{op_dir}/{hm}_violinplot_{type}.png')

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
            violinplot(df, hm, 'epi')


        ## NONEPI ENRICHED RBPS
        with open(f'{op_dir}/enriched_nonepi_{hm}.txt', 'r') as file:
            sfs = [line.strip() for line in file.readlines()]

        cols = sfs + ['label']
        df = features[cols]
        if len(sfs) >0:
            violinplot(df, hm, 'nonepi')
  
enriched_epi_nonepi()