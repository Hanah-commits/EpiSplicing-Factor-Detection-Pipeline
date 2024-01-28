import os
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import rankdata
import matplotlib.pyplot as plt


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


def plot_histogram(df, columns, status=0):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    col = 0
    for i in range(2):
        for j in range(2):
            ax = axes[i][j]
            ax.hist(df[columns[col]], bins=10, color='blue', alpha=0.5, label='{}'.format(columns[col]))
            ax.set_xlabel('P-values')
            ax.set_ylabel('count')
            ax.set_ylim([0, 250])
            leg = ax.legend(loc='upper left')
            leg.draw_frame(False)
            col += 1
    if status == 0:
        plt.suptitle('Non-adjusted P-value histogram')
    else:
        plt.suptitle('FDR adjusted P-value histogram')

    plt.show()


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



if __name__ == "__main__":

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    # read dPSI and M-values
    flanks = pd.read_csv('0_Files/Filtered_MValues_majiq.csv', delimiter='\t')

    # change geneID
    flanks['gene_id'] = flanks['gene_id'].str.split('.').str[0]  # ENSG00000116691.11 -> ENSG00000116691
    names = pd.read_csv('HelperFunctions/GeneID_Name.csv', delimiter='\t')
    names.columns = ['gene_id', 'gene']
    flanks = pd.merge(flanks, names, on='gene_id')
    del flanks['gene_id']
    flanks.rename(columns={'gene':'gene_id'}, inplace=True)

    flanks_meta = flanks.copy()
    flanks.drop(['seqid', 'strand'], axis=1, inplace=True)


    # FILTER 1: drop genes with less than 3 flanks
    flanks = flanks[flanks.groupby('gene_id').gene_id.transform(len) > 2]

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    AS_flanks = flanks[flanks.dPSI != 0]
    AS_flanks.replace(0, None, inplace=True) # to make comparison easier in next step

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    cols = ['gene_id'] + hms
    grouped = AS_flanks[cols].groupby('gene_id')
    non_epi = []
    for gene, group in grouped:
        if group[hms].isnull().all().all():
            non_epi.append(gene)

    filtered_flanks = flanks[~flanks['gene_id'].isin(non_epi)]
    filtered_flanks.set_index('flanks', inplace=True)

    # if one or more flanks of a gene have less than four peaks
    filtered_flanks.fillna(0, inplace=True)

    # Step 0: Find genes with strong DEU-DHM correlations
    coeff = filtered_flanks.groupby('gene_id').corr(method=pearsonr_coeff)

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

    strong_corr_genes = strong_corr(coeff)


    # Step 1: Obtain p values
    pval = filtered_flanks.groupby('gene_id').corr(method=pearsonr_pval)

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

    # # Step 3: Visualise p-value distribution
    # plot_histogram(pval, columns=["H3K4me3", "H3K27me3", "H3K9me3", "H3K27ac"])

    # Step 4: Adjust the p values using Benjamini-Hochberg method
    adj_pvals = adjust_pvalue(pval)

    # Step 5: Obtain genes where adjusted_pval < 0.05
    epigenes = find_epigenes(adj_pvals, strong_corr_genes)
    print('Epigenes ', len(epigenes))
    print('Non-Epigenes ', len(non_epi))


    # get flanks of epispliced genes
    flanks_meta[flanks_meta['gene_id'].isin(epigenes)].to_csv('0_Files/dPSI_Mval_epi_majiq.csv', sep='\t', index=False)

    # # get flanks of non-epispliced genes
    flanks_meta[flanks_meta['gene_id'].isin(non_epi)].to_csv('0_Files/dPSI_Mval_nonepi_majiq.csv', sep='\t', index=False)

    #clean-up
    os.remove('0_Files/pvals.csv')

    with open(f'0_Files/MAJIQ_epigenes.txt', 'w') as f:
        for line in list(set(epigenes)):
            f.write("%s\n" % line)

    with open(f'0_Files/MAJIQ_nonepigenes.txt', 'w') as f:
        for line in list(set(non_epi)):
            f.write("%s\n" % line)