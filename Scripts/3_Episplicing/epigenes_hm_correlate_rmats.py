import os
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import rankdata
import matplotlib.pyplot as plt
from venn import venn
import seaborn as sns
import numpy as np
import scipy as sp
from pathlib import Path


def pearsonr_pval(x, y):
    return pearsonr(x, y)[1]

def pearsonr_coeff(x, y):
    return pearsonr(x, y)[0]


def adjust_pvalue(df):
    pval_cols = df.columns.tolist()[1:]  # skipping gene-id column
    new_cols = []
    col_names = []
    for col in pval_cols:

        # get indices of null values
        na_idx = df[df[col].isnull()].index.tolist()

        # adjust non-null p values
        pvals = df[col].values.tolist()
        pvals = [x for x in pvals if str(x) != 'nan']
        adj_pval = p_adjust_bh(pvals).tolist()

        # insert null at original indices
        for idx in na_idx:
            adj_pval.insert(idx, None)

        new_cols.append(adj_pval)
        col_names.append(col + '_adj')

    # adjusted p values as new df
    df1 = pd.DataFrame(columns=col_names)
    df1['geneSymbol'] = df['geneSymbol'].values.tolist()
    for i in range(len(new_cols)):
        df1[col_names[i]] = new_cols[i]

    return df1


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def strong_corr(df):
    coeff_cols = df.columns.tolist()[1:]  # skipping gene-id column
    strong_corr_genes = []
    for col in coeff_cols:
        df [col] = df[col].abs()
        strong_corr_genes.extend((df.loc[df[col] >= 0.5, 'geneSymbol']).values.tolist())

    strong_corr_genes = sorted(list(set(strong_corr_genes)))
    return strong_corr_genes

    

def find_epigenes(df, corr_genes):
    pval_cols = df.columns.tolist()[:-1]  # skipping gene-id column
    epi_genes = []
    for col in pval_cols:
        epi_genes.extend((df.loc[df[col] <= 0.05, 'geneSymbol']).values.tolist())

    epi_genes = sorted(list(set(epi_genes) & set(corr_genes)))
    return epi_genes


def find_final_epigenes(df):

    corr_genes = (df.loc[df['coeff'].abs() >= 0.5, 'gene']).values.tolist()
    significant_genes = (df.loc[df['adj_pvals'].abs() <= 0.05, 'gene']).values.tolist()

    epi_genes = sorted(list(set(significant_genes) & set(corr_genes)))
    return epi_genes


def make_hm_plots(hm, both_hm_flanks):

    both_hm_flanks = both_hm_flanks.copy() # Make a copy to avoid the SettingWithCopyWarning
    dju_genes = list(set(both_hm_flanks['geneSymbol'].values.tolist()))
    both_hm_flanks["type"] = both_hm_flanks.apply(lambda row: 'dju' if row['dPSI'] != 0 else 'non-dju', axis=1)

    # impute missing data points
    both_hm_flanks.fillna(0,inplace=True)

    # get absolute DJU, DHM values
    both_hm_flanks.loc[:, 'dPSI'] = both_hm_flanks['dPSI'].abs()
    both_hm_flanks.loc[:, hm] = both_hm_flanks[hm].abs()
    
    # round off
    both_hm_flanks.loc[:, 'dPSI'] = both_hm_flanks['dPSI'].round(2)
    both_hm_flanks.loc[:, hm] = both_hm_flanks[hm].round(2)

    # create directory to save files
    Path(f'0_Files/RMATS/{hm}/').mkdir(parents=True, exist_ok=True)

    p_vals = []
    r_coeff = []
    filter_out = []
    # # get all genes
    for gene in dju_genes:
        gene_df = both_hm_flanks[both_hm_flanks['geneSymbol'] == gene]
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

    with open(f'0_Files/RMATS/{hm}_truepos_epigenes.txt', 'w') as f:
        for line in true_genes:
            f.write("%s\n" % line)

    both_hm_flanks = both_hm_flanks[both_hm_flanks['geneSymbol'].isin(true_genes)]
    both_hm_flanks.to_csv('0_Files/dPSI_Mval_epi_' + hm + '_rmats.csv', sep='\t', index=False)

    print(hm, '  ', len(true_genes))


    # get correlation plot of true epigenes
    for gene in true_genes:
        gene_df = both_hm_flanks[both_hm_flanks['geneSymbol'] == gene]
        corr_plot(gene_df, gene, hm, path=f'0_Files/RMATS/{hm}/')
    
    return true_genes
 

def corr_plot(gene_df, gene, hm, path):
    
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


def indiv_hms():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    tissue1 = d['tissue1'].capitalize()
    tissue2 = d['tissue2'].capitalize()

        # read dPSI and M-values
    flanks = pd.read_csv('0_Files/Filtered_MValues_rmats.csv', delimiter='\t')
    flanks_meta = flanks.copy()
    flanks.drop(['chr', 'strand'], axis=1, inplace=True)

    # FILTER 1: drop genes with less than 3 flanks
    flanks = flanks[flanks.groupby('geneSymbol').geneSymbol.transform(len) > 2]

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    AS_flanks = flanks[flanks.dPSI != 0]
    AS_flanks = AS_flanks.copy() # Make a copy to avoid the SettingWithCopyWarning
    AS_flanks.replace(0, None, inplace=True) # to make comparison easier in next step

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    cols = ['geneSymbol'] + hms
    grouped = AS_flanks[cols].groupby('geneSymbol')
    non_epi = []
    for gene, group in grouped:
        if group[hms].isnull().all().all():
            non_epi.append(gene)

    filtered_flanks = flanks[~flanks['geneSymbol'].isin(non_epi)].copy()
    filtered_flanks.set_index('flanks', inplace=True)

    # if one or more flanks of a gene have less than four peaks
    filtered_flanks.fillna(0, inplace=True)

    # Step 0: Find genes with strong DEU-DHM correlations
    coeff = filtered_flanks.groupby('geneSymbol').corr(method=pearsonr_coeff)

    # internal column filtering
    coeff.to_csv('0_Files/coeff.csv', sep='\t')
    coeff = pd.read_csv('0_Files/coeff.csv', delimiter='\t')
    coeff.drop(['Unnamed: 1', 'dPSI'], axis=1, inplace=True)
    
    # dropping p-values of hm-hm correlations
    coeff = coeff.iloc[::5, :]

    # drop genes where no dPSI-HM correlations exist
    coeff.dropna(subset=hms, how='all', inplace=True)

    # cleanup
    coeff.reset_index(inplace=True)
    del coeff['index']

    ## corr indiv-hms
    coeff_cols = coeff.columns.tolist()[1:]  # skipping gene-id column)
    coeff_genes = []
    for col in coeff_cols:
        tmp_df = coeff[['geneSymbol', col]]
        strong_corr_genes = strong_corr(coeff)
        coeff_genes.append(strong_corr_genes)
        

    # Step 1: Obtain p values
    pval = filtered_flanks.groupby('geneSymbol').corr(method=pearsonr_pval)

    # Step 2: Keep only relevant correlations

    # internal column filtering
    pval.to_csv('0_Files/pvals.csv', sep='\t')
    pval = pd.read_csv('0_Files/pvals.csv', delimiter='\t')
    pval.drop(['Unnamed: 1', 'dPSI'], axis=1, inplace=True)

    # dropping p-values of hm-hm correlations
    pval = pval.iloc[::5, :]

    # drop genes where no dPSI-HM correlations exist
    pval.dropna(subset=hms, how='all', inplace=True)

    # cleanup
    pval.reset_index(inplace=True)
    del pval['index']

    # Step 3: Adjust the p values using Benjamini-Hochberg method
    adj_pvals = adjust_pvalue(pval)

    # Step 4: Obtain genes where adjusted_pval < 0.05
    ## indiv-hms: find epigenes

    hm_cols = adj_pvals.columns.tolist()[:-1] 
    hm_epigenes = []
    i = 0
    for hm in hm_cols:
        tmp_df = adj_pvals[[hm, 'geneSymbol']]
        indiv_hm_epigenes = find_epigenes(tmp_df, coeff_genes[i])
        hm_epigenes.append(indiv_hm_epigenes)
        
        i += 1


    # get flanks of hm-specific epigenes
    i = 0
    true_epigenes = []
    for elem in hm_cols:
        hm = elem.split('_')[0]
        
        # get flanks of epispliced genes
        hm_df = flanks_meta[flanks_meta['geneSymbol'].isin(hm_epigenes[i])]

        # hm-specific corrplot
        true_epigenes.append(make_hm_plots(hm, hm_df))

        i+=1

    epigenes = list(set([item for items in true_epigenes for item in items]))
    print('Epigenes ', len(epigenes))
    print('Non-Epigenes ', len(non_epi))
    
    # get flanks of all epispliced genes
    flanks_meta[flanks_meta['geneSymbol'].isin(epigenes)].to_csv('0_Files/dPSI_Mval_epi_rmats.csv', sep='\t', index=False)

    # # get flanks of non-epispliced genes
    flanks_meta[flanks_meta['geneSymbol'].isin(non_epi)].to_csv('0_Files/dPSI_Mval_nonepi_rmats.csv', sep='\t', index=False)


    # remove unwnated
    os.remove('0_Files/pvals.csv')
    os.remove('0_Files/coeff.csv')

def plot_venn():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    tissue1 = d['tissue1'].capitalize()
    tissue2 = d['tissue2'].capitalize()

    hm_epigenes = []
    
    for hm in hms:
        with open(f'0_Files/RMATS/{hm}_truepos_epigenes.txt') as file:
            epigenes = [line.rstrip() for line in file]
            hm_epigenes.append(epigenes)
    
    # visualise overlap
    info = tissue1 + ' - ' + tissue2

    unique_epi =  set(i for j in hm_epigenes for i in j)
    title = info + ' : ' + str(len(unique_epi)) + ' Epigenes'

    labels = [elem.split('_adj')[0]+ ' : ' + str(len(hm_epigenes[ind])) for ind,elem in enumerate(hms)]
    values = [set(i) for i in hm_epigenes]

    data = dict(zip(labels, values))
    venn(data, cmap="plasma")

    plt.title(title)
    plt.savefig('0_Files/RMATS/hm_overlap_' + info + '.png')




if __name__ == "__main__":

    # # get epigenes and their correlation plots
    indiv_hms()

    # venn diagram of true epigenes
    plot_venn()