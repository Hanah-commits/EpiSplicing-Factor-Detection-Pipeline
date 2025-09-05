from pathlib import Path
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def ridgeplot_splice_site_scores(df, hm, ss_type, op_dir):

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
    g.set_xlabels(f"{ss_type}' Splice Site Strength Scores", fontsize=12)
    g.set_ylabels('Density', fontsize=12)
    g.add_legend()
    g.despine(left=True)
    plt.subplots_adjust(hspace=0.5)
    # g.fig.suptitle(f" Epispliced vs Non-epispliced Exons - {hm}", fontsize=14)
    g.fig.subplots_adjust(top=0.9)

    plt.savefig(f'{op_dir}/{hm}_ridgeplot_{ss_type}prime_splicesite.png', bbox_inches='tight', dpi=300)


def splice_site_strength_epi_nonepi(base_op_dir):
    """
    Input: Move output of in HelperFunctions/MAXENTSCAN/master.sh to base_op_dir/MAXENTSCAN/scores
    """

    op_dir = f"{base_op_dir}/MAXENTSCAN"
    Path(op_dir).mkdir(parents=True, exist_ok=True)

    epi_ss = pd.read_csv(f'{op_dir}/scores/epi_exons_splicesite_scores.bed', delimiter='\t', header=None)
    nonepi_ss = pd.read_csv(f'{op_dir}/scores/nonepi_exons_splicesite_scores.bed', delimiter='\t', header=None)

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
            ridgeplot_splice_site_scores(features, hm, score, op_dir)


def prep_log2_norm_counts(base_op_dir):
    """
    # Input: Move output of ./HelperFunctions/RBP_Expression/main.sh to base_op_dir/expression/counts
    """
    output_dir = f"{base_op_dir}/expression"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # read normalised counts tsv
    counts = pd.read_csv(f'{output_dir}/counts/tpm_values_rbps.tsv', delimiter='\t')

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
                clean_prefix = clean_prefix.replace("dermal", "") #ectodermal -> ecto
                clean_prefix = clean_prefix.replace("cell", "") #neuronalcell -> neuro
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
        with open(f'0_Files/Post-processing/epiRBPs/epiRBPs_{hm}.txt', 'r') as file:
            sfs.extend([line.strip() for line in file.readlines()])      

        sfs.sort()
        sf_counts = counts[counts['gene'].isin(sfs)]
        del sf_counts['gene']

        # Plot heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(
            sf_counts,
            cmap="seismic",
            annot=False,
            linewidths=1,
            cbar_kws={"label": "log2(TPM+1)"},
            xticklabels=True,
            yticklabels=sfs
        )
        # plt.title(f"Expression of Episplicing RBPs - {hm}", fontsize=10)
        plt.xlabel("")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{hm}.png', bbox_inches='tight', dpi=300)


def ridgeplot_exon_lengths(base_op_dir):
    """
    Input: Move output of in HelperFunctions/MAXENTSCAN/master.sh to base_op_dir/MAXENTSCAN/scores
    """
    op_dir = '0_Files/Post-processing/Analyses/MAXENTSCAN/scores'
    
    try:
        epi_ss = pd.read_csv(f'{op_dir}/epi_exons_splicesite_scores.bed', delimiter='\t', header=None)
        nonepi_ss = pd.read_csv(f'{op_dir}/nonepi_exons_splicesite_scores.bed', delimiter='\t', header=None)
    except Exception as e:
        print('Move output of in HelperFunctions/MAXENTSCAN/master.sh to op_dir/MAXENTSCAN/scores')
        print(e)

    epi_ss.columns = ['chr', 'start', 'stop', 'gene_name', 'type', 'strand', "5'ss", "3'ss"]
    nonepi_ss.columns = ['chr', 'start', 'stop', 'gene_name', 'type', 'strand', "5'ss", "3'ss"]

    epi_ss['Exon Class'] = 'Epispliced Exon'
    nonepi_ss['Exon Class'] = 'Non-epispliced Exon'

    all_ss = pd.concat([epi_ss, nonepi_ss], axis=0)

    # remove exons with both labels
    common_genes = list(set(epi_ss.gene_name.values.tolist()) & set(nonepi_ss.gene_name.values.tolist()))
    all_ss = all_ss[~((all_ss.gene_name.isin(common_genes)) & (all_ss['Exon Class'] == 'Non-epispliced Exon'))]
    all_ss = all_ss[~((all_ss.gene_name.isin(common_genes)) & (all_ss['Exon Class'] == 'Epispliced Exon'))]

    plt_op_dir = f"{base_op_dir}/epigenes"
    Path(plt_op_dir).mkdir(parents=True, exist_ok=True)

    hms = [  "H3K27ac","H3K4me3","H3K9me3", "H3K36me3", "H3K27me3"]

    combined_data = pd.DataFrame()
    for hm in hms:
        print(hm) 
        features = all_ss[all_ss['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))].copy()
        features.loc[:,'len'] = abs(features.loc[:, 'start'] - features.loc[:, 'stop'])

        features = features[features['Exon Class'] == 'Epispliced Exon']
        features['hm'] = hm
        combined_data = pd.concat([combined_data, features], ignore_index=True)

    
    # set color palette
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))
    custom_palette = { 'Non-epispliced Exon': '#e3e1d3'}
    custom_palette['Epispliced Exon'] = color_dict[hm]

    # change df to long format
    df_melted = combined_data.melt(id_vars=['Exon Class', 'hm'], value_vars=["len"], 
                        var_name='Variable', value_name='Value')

    plt.figure(figsize=(12, 6))
    for hm in hms:
        subset = df_melted[df_melted['hm'] == hm]

        # Density curve
        sns.kdeplot(data=subset, x='Value', label=hm, color=color_dict[hm], fill=True, alpha=0.4, bw_adjust=0.8, linewidth=1.5)

    # Add labels and legend
    # plt.title('Lengths of Epispliced Exons')
    plt.xlabel('Lengths (bp)', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.legend(title='Histone Mark', fontsize=10)
    plt.savefig(f'{plt_op_dir}/Exon_lengths.png', bbox_inches='tight', dpi=300)


def microexons_table(base_op_dir):
    """
    Input: Move output of in HelperFunctions/MAXENTSCAN/master.sh to base_op_dir/MAXENTSCAN/scores
    """
    
    op_dir = '0_Files/Post-processing/Analyses'

    try:
        epi_ss = pd.read_csv(f'{op_dir}/MAXENTSCAN/scores/epi_exons_splicesite_scores.bed', delimiter='\t', header=None)
        nonepi_ss = pd.read_csv(f'{op_dir}/MAXENTSCAN/scores/nonepi_exons_splicesite_scores.bed', delimiter='\t', header=None)
    except Exception as e:
        print(f'Move output of in HelperFunctions/MAXENTSCAN/master.sh to {op_dir}/MAXENTSCAN/scores')
        print(e)

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

    microgenes = []

    combined_data = pd.DataFrame()
    for hm in hms:
        print(hm) 
        features = all_ss[all_ss['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))].copy()
        features.loc[:,'len'] = abs(features.loc[:, 'start'] - features.loc[:, 'stop'])

        # ######## microexons ##############
        features = features[(features['len'] <= 30) & (features['Exon Class'] == 'Epispliced Exon')]
        microgenes.extend(features.gene_name.values.tolist())
        # print(features)
        
        combined_data = pd.concat([combined_data, features], ignore_index=True)

    combined_data['exon coordinates'] = combined_data['chr'].astype(str) + ':' + combined_data['start'].astype(str) + '-' + combined_data['stop'].astype(str)
    combined_data.drop_duplicates().to_csv(f'{op_dir}/epigenes/microexons.tsv', sep='\t', index=False)


def last_exon_epi_overlap(base_op_dir):
    """
    # epi_flanks_TSS.bed: /0_Files/Post-processing/epi_flanks.bed with overlap (bedtools intersect) with (alternative) last exons -> 200bp window of transcrp. termination site
    # epi_flanks_non_TSS.bed: /0_Files/Post-processing/epi_flanks.bed with overlap (bedtools intersect) with (alternative) last exons -> 200bp window of transcrp. termination site
    """
    
    op_dir = f"{base_op_dir}/epigenes"
    Path(op_dir).mkdir(parents=True, exist_ok=True)
    cols = ['chr', 'start', 'stop', 'feature', 'score', 'strand', 'gene_name', 'type']

    # epi flanks with and without overlap with (alternative) last exons -> 200bp window of transcrp. termination site
    epi_flanks_non_TSS = pd.read_csv(f'{op_dir}/exons/epi_flanks_non_TSS.bed', delimiter='\t', names=cols)
    epi_flanks_TSS = pd.read_csv(f'{op_dir}/exons/epi_flanks_TSS.bed', delimiter='\t', names=cols)

    # original flanks
    epi_flanks = pd.read_csv(f'0_Files/Post-processing/epi_flanks.bed', delimiter='\t', names=cols)
    nonepi_flanks = pd.read_csv(f'0_Files/Post-processing/nonepi_flanks.bed', delimiter='\t', names=cols)

    ##FILTER 1: Remove epi and nonepiflank overlaps
    common_genes = list(set(epi_flanks.gene_name.values.tolist()) & set(nonepi_flanks.gene_name.values.tolist()))
    epi_flanks_non_TSS = epi_flanks_non_TSS[~(epi_flanks_non_TSS.gene_name.isin(common_genes))]
    epi_flanks_TSS = epi_flanks_TSS[~(epi_flanks_TSS.gene_name.isin(common_genes))]

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    # figure and axes
    fig, axes = plt.subplots(1, len(hms), figsize=(3 * len(hms), 3.5))  
    # fig.suptitle('Percentage of Internal Epispliced Exons', fontsize=16)

    # shared legend
    handles = []
    labels = []

    for i, hm in enumerate(hms):
        epi_flanks_non_TSS_hm = epi_flanks_non_TSS[epi_flanks_non_TSS['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))].drop_duplicates()
        epi_flanks_TSS_hm = epi_flanks_TSS[epi_flanks_TSS['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))].drop_duplicates()

        size_of_groups = [len(epi_flanks_non_TSS_hm), len(epi_flanks_TSS_hm)]
    
        # set up subplot
        ax = axes[i]

        # custom colors
        colors = [color_dict[hm], 'oldlace']

        # pie chart
        wedges, text = ax.pie(
            size_of_groups,
            colors=colors,
            wedgeprops={'linewidth': 1, 'edgecolor': None}
        )

        # circle in the center to make it a donut chart
        center_circle = plt.Circle((0, 0), 0.8, color='white', fc='white')
        ax.add_artist(center_circle)

        largest_percentage = max(size_of_groups) / sum(size_of_groups) * 100
        ax.text(
            0, 0, f"{largest_percentage:.1f}%", 
            ha='center', va='center', fontsize=25, color=color_dict[hm]
        )

        # add the color and label to the legend list
        handles.append(wedges[0])  
        labels.append(hm)

    # shared legend
    fig.legend(
    handles, labels,
    title="", loc="upper center", bbox_to_anchor=(0.5, 0.15), ncol=len(hms), fontsize=15
    )
    # spacing between subplots and adjust layout
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.savefig(f'{op_dir}/Last_exon_overlap.png')
    plt.close()


def validation_epigenes(base_op_dir):

    hms = ['H3K27ac','H3K36me3', 'H3K9me3', 'H3K4me3']
    color_map = dict(zip(hms,["#9A71F8", "#B0D212", "#FF9900", "#ED588A"]))

    data = dict(zip(hms,[6,29,1,3])) # number of epigenes/hm in K562-HepG2 cell line from pr_K562_HepG2_timestamp/output.log
    labels = list(data.keys())
    sizes = list(data.values())
    
    colors = [color_map[label] for label in labels]

    # pie chart
    plt.figure(figsize=(6, 6))
    wedges, texts = plt.pie(
        sizes,
        labels=labels,
        colors=colors,
        # autopct='%1.1f%%',
        startangle=140
    )

    # value labels
    for i, text in enumerate(texts):
        text.set_text(f"{labels[i]}: {sizes[i]}")

    # Add a title
    # plt.title("Number of Epispliced Genes: HepG2-K562")
    plt.savefig('0_Files/Post-processing/Analyses/epigenes/Validation-Epigenes.png',bbox_inches='tight', dpi=300)


if __name__ == "__main__":

    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/Analyses"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    splice_site_strength_epi_nonepi(output_dir)
    prep_log2_norm_counts(output_dir)
    ridgeplot_exon_lengths(output_dir)
    microexons_table(output_dir)
    last_exon_epi_overlap(output_dir)
    validation_epigenes(output_dir)

    