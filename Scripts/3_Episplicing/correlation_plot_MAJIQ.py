import pandas as pd
import json
import os
import seaborn as sns
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pathlib import Path
from venn import venn



import numpy as np

def p_adjust_bh(p_values):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p_values = np.asfarray(p_values)
    num_tests = len(p_values)
    
    # Sort p-values in ascending order and keep track of the original indices
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    
    # Calculate the adjusted p-values using the Benjamini-Hochberg procedure
    adjusted_p_values = np.empty(num_tests)
    for i in range(num_tests):
        adjusted_p_values[i] = min(sorted_p_values[i] * (num_tests / (i + 1)), 1.0)
    
    # Restore the original order of p-values
    adjusted_p_values = adjusted_p_values[np.argsort(sorted_indices)]

    return adjusted_p_values



def strong_corr(df):
    coeff_cols = df.columns.tolist()[1:]  # skipping gene-id column
    strong_corr_genes = []
    for col in coeff_cols:
        df [col] = df[col].abs()
        strong_corr_genes.extend((df.loc[df[col] >= 0.5, 'gene']).values.tolist())

    strong_corr_genes = sorted(list(set(strong_corr_genes)))
    return strong_corr_genes


def find_final_epigenes(df):

    corr_genes = (df.loc[df['coeff'].abs() >= 0.5, 'gene']).values.tolist()
    significant_genes = (df.loc[df['adj_pvals'].abs() <= 0.05, 'gene']).values.tolist()

    epi_genes = sorted(list(set(significant_genes) & set(corr_genes)))
    return epi_genes


def remove_false_positives(both_flanks):

    ## some junctions are reported as non-dju in one LSV but as dju in another LSV: 
    ## we keep them as dju events (Eg, MLLT1, Neuor-H1, Exons 7,8,9)

    # Identify rows with duplicated 'flanks' values
    duplicates = both_flanks['flanks'].duplicated(keep=False)

    # Filter the df to keep only the duplicated rows where 'type' is 'dju' or where 'flanks' has unique values.
    filtered_both_flanks = both_flanks[(duplicates & (both_flanks['type'] == 'dju')) | ~duplicates]
    
    return filtered_both_flanks



def make_df(hm, control_flanks):

    print(hm)

    epi_file = '0_Files/dPSI_Mval_epi_' + hm + '_majiq.csv'
    both_hm_flanks = pd.read_csv(epi_file, delimiter='\t')
    dju_genes = list(set(both_hm_flanks['gene'].values.tolist()))
    both_hm_flanks["type"] = both_hm_flanks.apply(lambda row: 'dju' if row['dPSI'] != 0 else 'non-dju', axis=1)

    # impute missing data points
    both_hm_flanks.fillna(0,inplace=True)

    # get absolute DJU, DHM values
    both_hm_flanks['dPSI'] = both_hm_flanks['dPSI'].abs()
    both_hm_flanks[hm] = both_hm_flanks[hm].abs()
    
    # round off
    both_hm_flanks['dPSI'] = both_hm_flanks['dPSI'].round(2)
    both_hm_flanks[hm] = both_hm_flanks[hm].round(2)

    # create directory to save files
    Path(f'0_Files/{hm}/').mkdir(parents=True, exist_ok=True)

    p_vals = []
    r_coeff = []
    filter_out = []
    # # get all genes
    for gene in dju_genes:
        gene_df = both_hm_flanks[both_hm_flanks['gene'] == gene]
        r, p = sp.stats.pearsonr(x=gene_df['dPSI'], y=gene_df[hm])
        p_vals.append(p)
        r_coeff.append(r)
        # plot_3(gene_df, gene, hm, path=f'0_Files/{hm}/')

        ## FILTER 1: remve genes whee CS exons also have DHM peakss
        if ((gene_df[hm] != 0) & (gene_df['type'] == 'non-dju')).any():
            filter_out.append(gene)

    ## FILTER 2: Remove genes with weak correlation between DEU-DHM
    adj_pvals = p_adjust_bh(p_vals)
    # make temp df
    df = pd.DataFrame({'gene': dju_genes, 'coeff': r_coeff, 'adj_pvals': adj_pvals})
    true_genes = find_final_epigenes(df)
    true_genes = [g for g in dju_genes if g not in filter_out]

    with open(f'0_Files/{hm}_truepos_epigenes.txt', 'w') as f:
        for line in true_genes:
            f.write("%s\n" % line)

    both_hm_flanks = both_hm_flanks[both_hm_flanks['gene'].isin(true_genes)]
    both_hm_flanks.to_csv('0_Files/dPSI_Mval_epi_' + hm + '_majiq.csv', sep='\t', index=False)

    # get correlation plot of true epigenes
    for gene in true_genes:
        gene_df = both_hm_flanks[both_hm_flanks['gene'] == gene]
        plot_3(gene_df, gene, hm, path=f'0_Files/{hm}/')
    
 

def plot_3(gene_df, gene, hm, path):
    
    # plot deu vs dhm with regression line
    sns.regplot(x='dPSI',y=hm,data=gene_df,fit_reg=True, scatter_kws={'alpha':0.3})

    # call the scipy function for pearson correlation
    r, p = sp.stats.pearsonr(x=gene_df['dPSI'], y=gene_df[hm])
    # annotate the pearson correlation coefficient text to 2 decimal places
    ax = plt.gca()
    plt.text(0.45, 0.9, 'R = {:.2f}\np = {:.4f}'.format(r,p), transform=ax.transAxes)

    # num of flanks
    num = len(gene_df)
    # define titles and axes labels
    plt.title(f'{gene}  - {num} flanks')
    plt.xlabel('DJU')
    plt.ylabel(f'DHM - {hm}')

    plt.savefig(path+f'{gene}.png', bbox_inches='tight', dpi=300)
    plt.close()



def venn_plot(values, labels, info):

    # # get epigenes common between h3k27ac and h3k27me3
    # print(list(set(values[0]).intersection(values[1])))

    unique_epi =  set(i for j in values for i in j)
    title = info + ' : ' + str(len(unique_epi)) + ' Epigenes'

    labels = [elem.split('_adj')[0]+ ' : ' + str(len(values[ind])) for ind,elem in enumerate(labels)]
    values = [set(i) for i in values]

    data = dict(zip(labels, values))
    venn(data, cmap="plasma")

    plt.title(title)
    plt.savefig('0_Files/hm_overlap_' + info + '.png')


def plot_venn():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    tissue1 = d['tissue1'].capitalize()
    tissue2 = d['tissue2'].capitalize()

    hm_epigenes = []
    
    for hm in hms:
        with open(f'0_Files/{hm}_truepos_epigenes.txt') as file:
            epigenes = [line.rstrip() for line in file]
            hm_epigenes.append(epigenes)
    
    # visualise overlap
    venn_plot(hm_epigenes, hms, info = tissue1 + ' - ' + tissue2)
    


if __name__ == "__main__":

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    # read dPSI and M-values
    flanks = pd.read_csv('0_Files/Filtered_MValues_majiq.csv', delimiter='\t')

    for hm in hms:
        # make correlation    plot of true epoigenes
        make_df(hm, flanks)

    # venn diagram of true epigenes
    plot_venn()