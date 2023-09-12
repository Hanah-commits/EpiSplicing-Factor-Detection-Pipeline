import os
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import rankdata
import matplotlib.pyplot as plt
from venn import venn


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
    df1['gene_id'] = df['gene_id'].values.tolist()
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
        strong_corr_genes.extend((df.loc[df[col] >= 0.5, 'gene_id']).values.tolist())

    strong_corr_genes = sorted(list(set(strong_corr_genes)))
    return strong_corr_genes

    

def find_epigenes(df, corr_genes):
    pval_cols = df.columns.tolist()[:-1]  # skipping gene-id column
    epi_genes = []
    for col in pval_cols:
        epi_genes.extend((df.loc[df[col] <= 0.05, 'gene_id']).values.tolist())

    epi_genes = sorted(list(set(epi_genes) & set(corr_genes)))
    return epi_genes


def indiv_hms():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    tissue1 = d['tissue1'].capitalize()
    tissue2 = d['tissue2'].capitalize()

    dPSI = pd.read_csv('0_Files/Filtered_dPSI.csv', delimiter='\t')
    peaks = pd.read_csv('0_Files/Filtered_MValues.csv', delimiter='\t')
    dPSI.drop_duplicates(inplace=True)
    peaks.drop_duplicates(inplace=True)
    flanks = pd.merge(dPSI, peaks, how="outer")


    # FILTER 1: drop genes with less than 3 flanks
    flanks = flanks[flanks.groupby('gene_id').gene_id.transform(len) > 2]

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    cols = ['gene_id'] + hms
    grouped = flanks[cols].groupby('gene_id')
    non_epi = []
    for gene, group in grouped:
        if group[hms].isnull().all().all():
            non_epi.append(gene)

    filtered_flanks = flanks[~flanks['gene_id'].isin(non_epi)].copy()
    filtered_flanks.set_index('flanks', inplace=True)

    # if one or more flanks of a gene have less than four peaks
    filtered_flanks.fillna(0, inplace=True)

    # Step 0: Find genes with strong DEU-DHM correlations
    coeff = filtered_flanks.groupby('gene_id').corr(method=pearsonr_coeff)

    # internal column filtering
    coeff.to_csv('0_Files/coeff.csv', sep='\t')
    coeff = pd.read_csv('0_Files/coeff.csv', delimiter='\t')
    coeff.drop(['Unnamed: 1', 'mean_dpsi_per_lsv_junction'], axis=1, inplace=True)
    
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
        tmp_df = coeff[['gene_id', col]]
        strong_corr_genes = strong_corr(coeff)
        coeff_genes.append(strong_corr_genes)
        

    # Step 1: Obtain p values
    pval = filtered_flanks.groupby('gene_id').corr(method=pearsonr_pval)

    # Step 2: Keep only relevant correlations

    # internal column filtering
    pval.to_csv('0_Files/pvals.csv', sep='\t')
    pval = pd.read_csv('0_Files/pvals.csv', delimiter='\t')
    pval.drop(['Unnamed: 1', 'mean_dpsi_per_lsv_junction'], axis=1, inplace=True)

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
        tmp_df = adj_pvals[[hm, 'gene_id']]
        indiv_hm_epigenes = find_epigenes(tmp_df, coeff_genes[i])
        hm_epigenes.append(indiv_hm_epigenes)
        
        i += 1

    print('Epigenes')
    i = 0
    for elem in hm_cols:
        hm = elem.split('_')[0]
        print(hm, '  ', len(hm_epigenes[i]))
        # get flanks of epispliced genes
        flanks[flanks['gene_id'].isin(hm_epigenes[i])].to_csv('0_Files/dPSI_Mval_epi_' + hm + '.csv', sep='\t', index=False)

        i+=1


    print('Non-Epigenes ', len(non_epi))

    # visualise overlap
    venn_plot(hm_epigenes, hm_cols, info = tissue1 + ' - ' + tissue2)

    # remove unwnated
    os.remove('0_Files/pvals.csv')
    os.remove('0_Files/coeff.csv')


def venn_plot(values, labels, info):

    # get epigenes common between h3k27ac and h3k27me3
    print(list(set(values[0]).intersection(values[1])))


    unique_epi =  set(i for j in values for i in j)
    title = info + ' : ' + str(len(unique_epi)) + ' Epigenes'

    labels = [elem.split('_adj')[0]+ ' : ' + str(len(values[ind])) for ind,elem in enumerate(labels)]
    values = [set(i) for i in values]

    data = dict(zip(labels, values))
    venn(data, cmap="plasma")

    plt.title(title)
    plt.savefig('0_Files/hm_overlap_' + info + '.png')

if __name__ == "__main__":

    indiv_hms()