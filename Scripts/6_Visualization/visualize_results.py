from pathlib import Path
import pandas as pd
import numpy as np
import os
import json
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.utils import shuffle
from sklearn.decomposition import PCA
from upsetplot import UpSet, from_indicators
from scipy import stats
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests
import matplotlib.patches as mpatches
from itertools import combinations
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def plot_epigenes(base_op_dir):
    """
    Input: Epigene_numbers.tsv 
    Create table of num epigenes/hm in each condition using values from tail -10 ../Output/<tissue1>_<tissue2>_timestamp/output.log
    """

    output_dir = f"{base_op_dir}/epigenes"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    epigenes = pd.read_csv(f'{output_dir}/Epigene_numbers.tsv', delimiter='\t').iloc[:, : 11]
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
    ax.set_yticks(np.arange(0, 180, 25))
    # plt.title('Number of Epigenes Across All Embryonic Cell Line Pairs', fontsize=10)
    plt.xticks(rotation=45)
    plt.ylabel('Number of Epigenes', fontsize=12)
    plt.legend(loc='upper right')
    plt.savefig(f'{output_dir}/1_Epigenes.png',bbox_inches='tight', dpi=300)
    # plt.close()


def epigene_overlap(base_op_dir):

    output_dir = f"{base_op_dir}/epigenes"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    epigenes_df = pd.read_csv('0_Files/Post-processing/epi_flanks.bed', delimiter='\t', names=['chr', 'start', 'stop', 'feature', 'score', 'strand', 'gene', 'HM'])
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_map = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    data = {}
    for hm in hms:
        data[hm] = sorted(list(set(epigenes_df[epigenes_df['HM'].str.contains(hm)]['gene'].values.tolist())))

    # convert to dataframe
    all_elements = sorted(list(set().union(*data.values())))  # unique elements
    df = pd.DataFrame({key: [item in values for item in all_elements] for key, values in data.items()}, index=all_elements)

    # save as table
    df.to_csv(f'{output_dir}/epigenes_overlap.tsv', sep='\t')

    # convert to upset format
    upset_data = from_indicators(df.columns, df)

    # plot
    fig = plt.figure(figsize=(8, 6))
    upset = UpSet(upset_data, show_percentages=False, show_counts=True, sort_by="cardinality")
    plot = upset.plot(fig=fig)
    plot['intersections'].set_ylabel('Intersection size', fontsize=12)
    plot["totals"].set_xlabel("Num. Epigenes",  fontsize=12)
    # plt.suptitle("Overlap of Epispliced Genes",  fontsize=10, x=0.5, y=0.98, ha='center')
    

    #  y-axis labels (set names) and color them
    axes = fig.axes  #  all subplot axes
    set_labels = axes[1].get_yticklabels()  #  y-axis tick labels (set names)

    # Set tick label and dot colors to match based on the color_map
    for label in set_labels:
        text = label.get_text()
        if text in color_map:
            label.set_color(color_map[text])

    plt.savefig(f'{output_dir}/Epigenes_overlap.png',bbox_inches='tight', dpi=300)
    plt.close()


def heatmap_imptRBPs_binding(hm, base_op_dir):

    output_dir = f"{base_op_dir}/heatmaps"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    #get epi and nonepiRBPs to plot
    features_to_plot = []
    epi_RBPs = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPs/epiRBPs_{hm}.txt", "r")
    features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))
    epi_RBPs = features_to_plot.copy()
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPs/nonepiRBPs_{hm}.txt", "r")
    features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features = features.sort_values(by="label", key=lambda col: col != "epigene")
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
    tick_label_color_dict = dict(zip(hms,["#6034C6", "#19A4E0", "#56B50E", "#DA5E0C", "#E11156"]))

    label_colors = dict(zip(['epigene', 'non-epigene'], [color_dict[hm], "grey"]))
    row_colors = features.index.map(label_colors)

    cluster = sns.clustermap(features[features_to_plot], annot=None, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='bwr',  cbar_pos=(0.9, .2, .03, .4))
    heatmap = cluster.ax_heatmap

    # remove y-axis ticks and tick labels
    heatmap.set_yticks([])
    heatmap.set_yticklabels([])

    # set x-axis tick labels
    heatmap.set_xticks(range(len(features_to_plot)))
    heatmap.set_xticklabels(features_to_plot, size=12, rotation=90)

    #  customize x-axis tick label color
    for label in heatmap.get_xticklabels():
        text = label.get_text()
        label.set_color(tick_label_color_dict[hm] if text in epi_RBPs else 'dimgrey')

    # add title and axes labels
    # heatmap.set_title(f'{hm} : Predicted Binding Scores of Episplicing and Non-episplicing RBPs', fontsize=10)
    heatmap.set_ylabel(f'Flanks of Epispliced and Non-epispliced Exons', fontsize =16, labelpad=-520)
    cax = cluster.cax
    cax.tick_params(labelsize=12)
    cax.set_ylabel("Predicted Binding Scores", rotation=90, labelpad=10, fontsize=16)
    plt.savefig(f"{output_dir}/bindingscores_{hm}.png", bbox_inches='tight', dpi=300, pad_inches=0.2)
    plt.close()


def process_dataframe(df, hms):
    # Internal column filtering
    df = df.drop(['Unnamed: 1', 'dPSI'], axis=1)

    # Dropping p-values/coeffs of hm-hm correlations
    df = df.groupby('gene_name').last()

    # Reset the index to ensure it starts from 0
    df['gene_name'] = df.index
    df.reset_index(drop=True, inplace=True)

    # Drop genes where no dPSI-HM correlations exist
    df.dropna(subset=hms, how='all', inplace=True)
    
    # Cleanup
    df.reset_index(drop=True, inplace=True)

    # rearrange columns
    df = df[['gene_name'] + hms]

    return df

    
def plot_correlation_manhattan(base_op_dir, tool = 'RMATS'):

    prefix = tool.lower()

    op_dir = f"{base_op_dir}/epigenes"
    Path(op_dir).mkdir(parents=True, exist_ok=True)

    with open('paths.json') as f:
            data = json.load(f)

    # STEP 1: Get the list of histone modifications available in the study
    hms = set()
    for process in data['list_of_processes']:
        hms.update(data[process]['Histone modifications'])

    # STEP 2: Get the corr coeffs and pvalues
    for hm in list(hms):
        print('\n', hm)
        # Get the list of output directories for the current hm
        processes = [process for process in data['list_of_processes'] if hm in data[process]['Histone modifications']]
        output_directories = [data[process]['Output directory'] for process in processes]

        ## STEP 2A: Get corr coeffs and pvalues in all HM analyses
        hm_epigenes = {}
        corr_dfs_hm = {}
        for dir in output_directories:

            # key
            tissues = os.path.basename(os.path.dirname(dir)).split('_')[1:3]
            tissues_substr = ['endo', 'ecto', 'meso', 'neuro', 'H1']
            tissues = [next((sub for sub in tissues_substr if sub in tissue), tissue) for tissue in tissues]
            tissues = '-'.join(tissues).title()

            # get epigenes
            file_path = glob.glob(f'{dir}*0_Files/{tool}/{hm}/{hm}_truepos_epigenes.txt') # only DEU != 0 , DHM != 0
            try: ## ../0_Files/{tool}/H3K27ac/H3K27ac_truepos_epigenes.txt
                with open(file_path[0], 'r') as file:
                    true_epigenes= [gene.strip() for gene in file]              
                    hm_epigenes[tissues] = true_epigenes
            except: 
                continue

            try:
                file_path = glob.glob(f'{dir}*0_Files/{tool}')
                coeff = process_dataframe(pd.read_csv(f'{file_path[0]}/coeff.csv',delimiter='\t'), [hm])
                coeff.columns = ['gene_name', 'R']
                pval =process_dataframe(pd.read_csv(f'{file_path[0]}/pvals.csv', delimiter='\t'), [hm])
                pval.columns = ['gene_name', 'pval']
                pval["pval"] = multipletests(pval["pval"].values, method="fdr_bh")[1]
                # df of R and pval
                corr_df = pd.merge(coeff,pval, on='gene_name')
                corr_df.loc[corr_df["R"] == 0, "pval"] = 1.0 #double-check
                corr_dfs_hm[tissues] = corr_df
            except Exception as e:
                print(e)
                continue
        
        ## STEP 2B: Plot
        plot_manhattan(corr_dfs_hm, hm_epigenes, op_dir, hm)


def plot_manhattan(corr_dfs_hm, hm_epigenes, op_dir,  hm):

        # reorder dict
        col_order = ['Ecto-H1', 'Ecto-Meso', 'Endo-Ecto', 'Ecto-Neuro', 'Endo-Meso', 'Endo-Neuro', 'Endo-H1', 'Meso-Neuro', 'Meso-H1', 'Neuro-H1']
        corr_dfs_hm = {k: corr_dfs_hm[k] for k in col_order if k in corr_dfs_hm.keys()}

        # Thresholds
        r_thresh = 0.5
        p_thresh = 0.05

        def categorize(row):
            
            if row["R"] >= 0.5:
                if row["pval"] <= 0.05:
                    if row["gene_name"] in hm_epigenes[row['source']]:
                        return f"R >= 0.5 & ¬(¬DEU & DHM:{hm})"
                    else:
                        return "R >= 0.5"
                else:
                    return "R >= 0.5"
            else:
                return "R < 0.5"

        combined = pd.concat(corr_dfs_hm, names=["source"]).reset_index(level=0)
        combined["category"] = combined.apply(categorize, axis=1)

        # dataset mapping for x-axis
        datasets = combined["source"].unique()
        x_map = {name: i for i, name in enumerate(datasets)}

        # define markers/colors for categories
        marker_map = {
            "R < 0.5" : "x",
            "R >= 0.5": ".",
            f"R >= 0.5 & ¬(¬DEU & DHM:{hm})": "D"  # diamond
        }

        size_map = {
            "R < 0.5" : 15,
            "R >= 0.5": 10,
            f"R >= 0.5 & ¬(¬DEU & DHM:{hm})": 15
        }

        alpha_map = {
            "R < 0.5" : 0.5,
            "R >= 0.5": 0.5,
            f"R >= 0.5 & ¬(¬DEU & DHM:{hm})": 0.5
        }

        # plot with jitter
        plt.figure(figsize=(10,6))
        for cat in combined["category"].unique():
            df_sub = combined[combined["category"] == cat]
            x_vals = df_sub["source"].map(x_map) + np.random.uniform(-0.25, 0.25, size=len(df_sub))
            sc = plt.scatter(
                x_vals, 
                df_sub["R"], 
                marker=marker_map[cat],
                c=df_sub['pval'],
                cmap="coolwarm_r", 
                alpha=alpha_map[cat],
                s=size_map[cat],
                label=cat,
                vmin=0,
                vmax=1.0
            )

        # threshold lines
        plt.axhline(0.5, color="#BB4430", linestyle="--")

        # colorbar
        cbar = plt.colorbar(sc)   # use the same scatter handle
        cbar.set_label("Adjusted p-value", fontsize=16)

        # formatting
        # Set y-axis ticks in intervals of 0.10
        # plt.yticks(np.arange(-1, 1.0, 0.10))
        plt.xticks(range(len(datasets)), datasets, rotation=45, fontsize=13, ha="right", color='black')
        plt.ylim(-1.05, 1.05)
        plt.ylabel("Pearson Correlation Coefficient R", fontsize=16)
        # plt.legend()

        # Legend below plot
        handles = [plt.Line2D([], [], marker=marker_map[cat], color="black", linestyle="", label=cat) for cat in ["R < 0.5", "R >= 0.5", f"R >= 0.5 & ¬(¬DEU & DHM:{hm})"]]
        plt.legend(
            fontsize=12,
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.1),   # below the x-axis
            ncol=3,                        # put categories in one row
            frameon=True
        )

        plt.tight_layout()
        plt.savefig(f'{op_dir}/{hm}_manhattan.png', bbox_inches='tight', dpi=500)


def RBP_binding_comparison_epispliced_vs_nonepispliced(base_op_dir):

    op_dir = f"{base_op_dir}/binding_scores"
    Path(op_dir).mkdir(parents=True, exist_ok=True)

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t').set_index('label')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t').set_index('label')

    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)
    
    # compute row means
    def row_means(df, cols, ax=0):
        return df[cols].mean(axis=ax)
    
    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))
    edgecolor_dict = dict(zip(hms,["#6034C6", "#19A4E0", "#56B50E", "#DA5E0C", "#E11156"]))

    for hm in hms:
        #get epi and nonepiRBPs to plot
        features_to_plot = []
        rbps_file = open(f"0_Files/Post-processing/epiRBPs/epiRBPs_{hm}.txt", "r")
        features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))
        epi_RBPs = features_to_plot.copy()
        rbps_file = open(f"0_Files/Post-processing/nonepiRBPs/nonepiRBPs_{hm}.txt", "r")
        features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))
        nonepi_RBPs = [rbp for rbp in features_to_plot if rbp not in epi_RBPs]
            
        # extract features for hm
        features_hm = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
        features_hm = features_hm[features_to_plot]
        features_hm = features_hm.applymap(lambda val: 0 if val < 2 else val) # keep only strong binding events

        results = []

        # STEP 1: Compile comparisons
        pairs = {
            "Epi exon flanks vs Non-epi exon flanks (Epi RBPs)": (row_means(features_hm.loc["epigene"], epi_RBPs), row_means(features_hm.loc["non-epigene"], epi_RBPs)),
            "Epi RBPs vs Non-epi RBPs (Epi exon flanks)": (row_means(features_hm.loc["epigene"], epi_RBPs), row_means(features_hm.loc["epigene"], nonepi_RBPs)),
            "Epi exon flanks vs Non-epi exon flanks (Non-epi RBPs)": (row_means(features_hm.loc["epigene"], nonepi_RBPs), row_means(features_hm.loc["non-epigene"], nonepi_RBPs)),
            "Epi RBPs vs Non-epi RBPs (Non-epi exon flanks)": (row_means(features_hm.loc["non-epigene"], epi_RBPs), row_means(features_hm.loc["non-epigene"], nonepi_RBPs))
        }

        # STEP 2: Welsh's Test
        for label, (x, y) in pairs.items():
            stat, p = stats.ttest_ind(x, y, equal_var=False)
            results.append({"Comparison": label, "p-value": p, "group1": x, "group2": y})
        results_df = pd.DataFrame(results)
        results_df["p-value"] = multipletests(results_df["p-value"], method="fdr_bh")[1]


        # STEP 3: Plot all comparisons
        plot_data = []
        for r in results:
            for val in r["group1"]:
                plot_data.append([r["Comparison"], "Group1", val])
            for val in r["group2"]:
                plot_data.append([r["Comparison"], "Group2", val])

        plot_df = pd.DataFrame(plot_data, columns=["Comparison","Group","Value"])

        g = sns.catplot(
            data=plot_df,
            x="Group", y="Value",
            col="Comparison",
            kind="box",
            col_wrap=2, sharey=False,
            showcaps=True, boxprops={'facecolor':'none'},
            whiskerprops={'linewidth':1.5}, medianprops={'color':'grey'},
            whis=2.0
        )
     
        for ax, comp in zip(g.axes.flat, results_df["Comparison"]):
            # get p-value
            p_val = results_df.loc[results_df["Comparison"] == comp, "p-value"].values[0]

            # Create pairs for pvalue annotation
            pairs = [("Group1", "Group2")]

            annotator = Annotator(ax, pairs, data=plot_df[plot_df["Comparison"]==comp], x="Group", y="Value")
            annotator.configure(text_format="simple", loc='inside', fontsize=12)
            annotator.set_pvalues([p_val])
            annotator.annotate()

            # set ylimit
            ax.set_ylim(0, 4)
            ax.set_yticks(np.arange(0, 4.5, 0.5))

            # include datapoints    
            sns.swarmplot(
                data=plot_df.loc[plot_df["Comparison"] == ax.get_title().split(" = ")[-1]],
                x="Group", y="Value",
                ax=ax,
                color="black",
                size=3,
                dodge=True
            )


        # custom coloring rules per comparison
        color1 = color_dict[hm]
        color2 = '#e3e1d3'
        color3 = "#B6B39A"
        color4 = edgecolor_dict[hm]

        settings_dict = {}
        for comp in set(plot_df["Comparison"].values.tolist()):
            
            if "Epi exon flanks vs Non-epi exon flanks (Epi RBPs)" in comp:
                facecolors = [color1, color2]
                edgecolors = [color4, color4]
            elif "Epi exon flanks vs Non-epi exon flanks (Non-epi RBPs)" in comp:
                facecolors = [color1, color2]
                edgecolors = [color3, color3]
            elif "Epi RBPs vs Non-epi RBPs (Epi exon flanks)" in comp:
                facecolors = [color1, color1]
                edgecolors = [color4, color3]
            elif "Epi RBPs vs Non-epi RBPs (Non-epi exon flanks)" in comp:
                facecolors = [color2, color2]
                edgecolors = [color4, color3]
            elif "Epi exon flanks (Epi RBPs) vs Non-epi exon flanks (Non-epi RBPs)" in comp:
                facecolors = [color1, color2]
                edgecolors = [color4, color3]
            elif "Epi exon flanks (Non-epi RBPs) vs Non-epi exon flanks (Epi RBPs)" in comp:
                facecolors = [color1, color2]
                edgecolors = [color3, color4]

            settings_dict[comp] = {
                'fc' : facecolors,
                'ec' : edgecolors
            }


        for comp, colors in settings_dict.items():        
            ax = g.axes_dict[comp]  # get the exact axis for current comparison
            facecolors = colors['fc']
            edgecolors = colors['ec']

            boxes = [p for p in ax.patches if isinstance(p, mpatches.PathPatch)]
    
            if len(boxes) != len(facecolors) or len(boxes) != len(edgecolors):
                raise ValueError(f"Mismatch in number of patches and colors for {comp}")
            
            # adjust colors of boxes, whiskers and caps
            for i, (patch, fc, ec) in enumerate(zip(boxes, facecolors, edgecolors)):

                patch.set_facecolor(fc)
                patch.set_edgecolor(ec)
                patch.set_linewidth(3)
                patch.set_linestyle("--")

                line_offset = i * 6
                whiskers = ax.lines[line_offset:line_offset + 2]
                for whisker in whiskers:
                    whisker.set_color(ec)
                    whisker.set_linestyle("--")
                    whisker.set_linewidth(3)

                caps = ax.lines[line_offset + 2 : line_offset + 4]
                for cap in caps:
                    cap.set_color(ec)
                    cap.set_linestyle("--")
                    cap.set_linewidth(3)

        # legend
        legend_elements = [
            mpatches.Patch(facecolor=color1, edgecolor=color1, alpha=0.7, label=f"{hm}: Epispliced exon flanks"),
            mpatches.Patch(facecolor=color2, edgecolor=color2, label="Non-epispliced exon flanks"),
            mpatches.Patch(facecolor="white", edgecolor=color4, linewidth=1, linestyle="--", label=f"{hm}: Episplicing RBPs"),
            mpatches.Patch(facecolor="white", edgecolor=color3, linewidth=1,  linestyle="--", label="Non-episplicing RBPs"),
        ]
        g.fig.legend(handles=legend_elements, loc="upper center", ncol=4, frameon=True, bbox_to_anchor=(0.5, 1.05))
        g.fig.subplots_adjust(top=0.85)

        # Round limits to nearest multiple of 0.5
        ymin, ymax = plot_df["Value"].min(), plot_df["Value"].max()
        ymin = 0.5 * np.floor(ymin / 0.5)
        ymax = 0.5 * np.ceil(ymax / 0.5)

        # remove titles, labels and ticks, adjust ylim
        for ax in g.axes.flat:
            ax.set_title("")
            ax.set_ylabel("")
            ax.set_xticks([])
            ax.tick_params(axis="y", labelsize=12)
            ax.set_ylim(ymin, ymax)

        g.set_axis_labels("", " ")

        # Adjust x-coordinate to add padding
        g.fig.text(
            x=0.01,
            y=0.5,
            s="Average Predicted Binding Scores",
            va="center",
            ha="center",
            rotation="vertical",
            fontsize=18
        )

        plt.tight_layout()

        plt.savefig(f'{op_dir}/{hm}_epispliced_vs_nonepispliced.png',bbox_inches='tight', dpi=300)


def imptRBP_binding_comparison(RBP_type, base_op_dir):

    op_dir = f"{base_op_dir}/binding_scores"
    Path(op_dir).mkdir(parents=True, exist_ok=True)
    
    features1 = pd.read_csv('0_Files/Post-processing/features_all_exon_classes_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_exon_classes_47.csv', delimiter='\t')

    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)
    features['label'] = np.where(features.label == 'epigene','DEU & DHM', 
                        np.where(features.label == 'non-epigene', 'DEU & ¬DHM', '¬DEU & DHM'))
    features.set_index(['label'], inplace=True)
   

    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))
    if RBP_type == 'epi':
        edgecolor_dict = dict(zip(hms,["#6034C6", "#19A4E0", "#56B50E", "#DA5E0C", "#E11156"]))
    else:
        edgecolor_dict = dict(zip(hms,["#B6B39A"]*len(hms)))

    for hm in hms:

        print(hm)

        #get RBPs to plot
        features_to_plot = []
        epi_RBPs = []
        rbps_file = open(f"0_Files/Post-processing/{RBP_type}RBPs/{RBP_type}RBPs_{hm}.txt", "r")
        features_to_plot.extend(sorted([rbp for rbp in rbps_file.read().split('\n') if rbp]))
        class_var = f'{hm} : Exon_class'
            
        # extract features for hm
        features_hm = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
        features_hm = features_hm.drop_duplicates()
        features_hm = features_hm[features_to_plot]
        features_hm= features_hm.map(lambda val: 0 if val < 2 else val) # keep only strong binding events
        
        ## Compare distributions
        exon_classes = features_hm.index.unique()
        pairs = list(combinations(exon_classes, 2))
        pvals = []
        for g1, g2 in pairs:
            vals1= features_hm[features_hm.index == g1].mean(axis=0)
            vals2= features_hm[features_hm.index == g2].mean(axis=0)
            stat, p = stats.ttest_ind(vals1, vals2, equal_var=False)
            pvals.append(p)
        pvals_adj= multipletests(pvals, method="fdr_bh")[1]

        ## Plot boxplot
        custom_palette = dict(zip(exon_classes.tolist(), [color_dict[hm], '#e3e1d3', color_dict[hm]]))  #custom cmap
        # datapoints to plot
        plot_values = []
        for exon_class in exon_classes:
            plot_values.append(features_hm[features_hm.index == exon_class].drop_duplicates().mean(axis=0))

        flanks_mean = pd.DataFrame(plot_values).T
        flanks_mean.columns = exon_classes
        flanks_mean = flanks_mean.reset_index().melt(id_vars="index", value_vars=exon_classes,
                                var_name=class_var, value_name="Binding_mean")

        plt.figure(figsize=(6,5))
        ax = sns.boxplot(data=flanks_mean, x=class_var, y="Binding_mean",
        palette=custom_palette, saturation=1, showfliers=False, width=0.6, whis=2.0,
        boxprops=dict(edgecolor=edgecolor_dict[hm], linewidth=2, linestyle="--"),
        medianprops=dict(color="grey"),
        whiskerprops=dict(color=edgecolor_dict[hm], linewidth=2, linestyle="--"),
        capprops=dict(color=edgecolor_dict[hm], linewidth=2, linestyle="--"))
        sns.swarmplot(data=flanks_mean, x=class_var, y="Binding_mean", 
                    color="black",
                    size=3,
                    dodge=True,
                    ax = ax)

        ## annotate with pvalues
        annotator = Annotator(ax, pairs, data=flanks_mean, x=class_var, y="Binding_mean")
        annotator.configure(text_format="simple", loc='inside', fontsize=14)
        annotator.set_pvalues(pvals_adj.tolist())
        annotator.annotate()

        # remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # adjust ylimits
        ymin = 0.0
        ymax = 4.0
        ymin = np.floor(flanks_mean['Binding_mean'].min() * 4) / 4 
        ax.set_ylim(ymin, ymax)
        # ax.set_yticks(np.arange(ymin, ymax + 0.25, 0.25))
        ax.tick_params(axis="both", labelsize=12)
       
        plt.xlabel(f"Exon Classes - {hm}", fontsize=14, labelpad=15)
        plt.ylabel(f"Average Predicted Binding Scores", fontsize=14, labelpad=15)
        plt.tight_layout()
        plt.savefig(f'{op_dir}/{hm}_{RBP_type}RBP_binding.png', bbox_inches='tight', dpi=300)


def plot_epiflanks(base_op_dir):
    """
    Input: Epigene_numbers.tsv 
    Create table of num epiflanks/hm in each condition using values from std output of post-rbp.py
    """
    output_dir = f"{base_op_dir}/epigenes"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    epiflanks = pd.read_csv(f'{output_dir}/Epiflanks_numbers.tsv', delimiter='\t')
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
    # plt.title('Number of Flanks of Epispliced and Non-epispliced Exons Available for All Histone Marks', fontsize=10)
    plt.xticks(rotation=45)
    plt.ylabel('Number of Exon Flanks', fontsize=12)
    plt.xlabel('')
    plt.savefig(f'{output_dir}/EpiFlanks.png',bbox_inches='tight', dpi=300)


def PCA_plot(hm, base_op_dir):

    op_dir = f"{base_op_dir}/PCA"
    Path(op_dir).mkdir(parents=True, exist_ok=True)

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
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
    ax.set_xlabel('Number of Components', fontsize=14)
    ax.set_ylabel('Cumulative Variance Explained', fontsize=14)
    # ax.set_title(f'Number of Principal Components Needed to Explain Variance - {hm}', fontsize=10)
    ax.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.legend(loc='lower right', fontsize=12)

    # edit x-axis ticks
    ax.set_xticks(np.arange(1, n_components + 1, step=1))
    ax.tick_params(axis='x', rotation=90)
    ax.tick_params(axis='both', labelsize=10)

    plt.savefig(f"{op_dir}/PCA_{hm}_explained_variance.png", bbox_inches='tight', dpi=300)
    plt.close()
    

def heatmap_epiRBPs_correlation(hm, base_op_dir):

    output_dir = f"{base_op_dir}/heatmaps"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)    

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info

    rbps_file = open(f"0_Files/Post-processing/epiRBPs/epiRBPs_{hm}.txt", "r")
    epi_rbps = [rbp for rbp in rbps_file.read().split('\n') if rbp]
    features = features[epi_rbps]
    features = features.applymap(lambda val: 0 if val < 2 else val) # keep only strong binding events

    # Plot corr coeffs
    corr = features.corr()
    ax = sns.heatmap(corr, annot=False, cmap='coolwarm')
    # plt.title(f'Correlation Coefficients of Episplicing RBPs Associated with {hm}', fontsize=8)
    plt.xticks(fontsize=5.5, rotation=90)
    plt.yticks(fontsize=5.5, rotation=0)
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=6)
    cbar = ax.collections[0].colorbar
    cbar.set_label("Pearson's R", fontsize=10)
    
    plt.savefig(f"{output_dir}/corr_{hm}.png",dpi=300) #bbox_inches='tight')
    plt.close()


def heatmap_allRBPs(hm, base_op_dir):

    output_dir = f"{base_op_dir}/heatmaps/corr_all_features"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print(f'\n\n{hm}')

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
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

    # Plot corr coeffs
    corr = features[sf].corr().abs()

    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(corr, annot=False, xticklabels=corr.columns, yticklabels=corr.index, cmap='coolwarm', linewidths=0.1)
    # plt.title(f'Correlation of Predicted RBP Binding Scores - {hm}', fontsize=8)

    plt.xticks(np.arange(len(sf))+ 0.5, labels=corr.columns, fontsize=2, rotation=90)
    plt.yticks(np.arange(len(sf))+ 0.5, labels=corr.columns, fontsize=2, rotation=0)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=5) 
    cbar.set_label("Pearson's R", fontsize=12)

    plt.savefig(f"{output_dir}/corr_{hm}.png",dpi=300) #bbox_inches='tight')
    plt.close()


if __name__ == "__main__":

    output_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/Analyses"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # plot_epigenes(output_dir) # number of epigenes across all ten analyses
    # plot_epiflanks(output_dir) 

    epigene_overlap(output_dir)
    plot_correlation_manhattan(output_dir)
    RBP_binding_comparison_epispliced_vs_nonepispliced(output_dir)
    for RBP_type in ['epi', 'nonepi']:
        imptRBP_binding_comparison(RBP_type, output_dir) 

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:
        heatmap_imptRBPs_binding(hm, output_dir)
        PCA_plot(hm, output_dir)
        heatmap_epiRBPs_correlation(hm, output_dir)
        heatmap_allRBPs(hm, output_dir)