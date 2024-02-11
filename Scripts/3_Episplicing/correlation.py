import os
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
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
    df1['gene_name'] = df['gene_name'].values.tolist()
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
        strong_corr_genes.extend((df.loc[df[col] >= 0.5, 'gene_name']).values.tolist())

    strong_corr_genes = sorted(list(set(strong_corr_genes)))
    return strong_corr_genes

    

def find_epigenes(df, corr_genes):
    pval_cols = df.columns.tolist()[:-1]  # skipping gene-id column
    epi_genes = []
    for col in pval_cols:
        epi_genes.extend((df.loc[df[col] <= 0.05, 'gene_name']).values.tolist())

    epi_genes = sorted(list(set(epi_genes) & set(corr_genes)))
    return epi_genes


def correlation(dir):

    print(f'Epigene Deection: {dir}')

    file = dir.lower()

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    # read dPSI and M-values
    flanks = pd.read_csv(f'0_Files/{dir}/DEU_DHM_{file}_flanks.tsv', delimiter='\t')
    del flanks['geneSymbol']

    # unique row index
    flanks['idx'] =  flanks['gene_name'] + flanks[['flank_start', 'flank_end']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)

    flanks_meta = flanks.copy()
    flanks.drop(['chr', 'strand'], axis=1, inplace=True)

    # FILTER 1: drop # genegenes with less than 3 flanks
    flanks = flanks[flanks.groupby('gene_name').gene_name.transform(len) > 2]

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    AS_flanks = flanks[flanks.dPSI != 0]
    AS_flanks = AS_flanks.copy()  # Make a copy to avoid the SettingWithCopyWarning
    AS_flanks.replace(0, None, inplace=True) # to make comparison easier in next step

    # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    cols = ['gene_name'] + hms
    grouped = AS_flanks[cols].groupby('gene_name')
    non_epi = []
    for gene, group in grouped:
        if group[hms].isnull().all().all():
            non_epi.append(gene)

    filtered_flanks = flanks[~flanks['gene_name'].isin(non_epi)]
    filtered_flanks.set_index('idx', inplace=True)

    # if one or more flanks of a gene have less than four peaks
    filtered_flanks.fillna(0, inplace=True)

    # Step 0: Find genes with strong DEU-DHM correlations

    ## remove unnecessary cols
    filtered_flanks.drop(columns=['feature', 'score', 'flank_start', 'flank_end'], inplace=True)
    coeff = filtered_flanks.groupby('gene_name').corr(method=pearsonr_coeff)

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
    pval = filtered_flanks.groupby('gene_name').corr(method=pearsonr_pval)

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

    # Step 5: Obtain genes where adjusted_pval < 0.05
    epigenes = find_epigenes(adj_pvals, strong_corr_genes)
    print('Epigenes ', len(epigenes))
    print('Non-Epigenes ', len(non_epi))


    # get flanks of epispliced genes
    flanks_meta[flanks_meta['gene_name'].isin(epigenes)].to_csv(f'0_Files/{dir}/dPSI_Mval_epi_{file}.csv', sep='\t', index=False)

    # # get flanks of non-epispliced genes
    flanks_meta[flanks_meta['gene_name'].isin(non_epi)].to_csv(f'0_Files/{dir}/dPSI_Mval_nonepi_{file}.csv', sep='\t', index=False)

    #clean-up
    os.remove('0_Files/pvals.csv')

    with open(f'0_Files/{dir}/{file}_epigenes.txt', 'w') as f:
        for line in list(set(epigenes)):
            f.write("%s\n" % line)

    with open(f'0_Files/{dir}/{file}_nonepigenes.txt', 'w') as f:
        for line in list(set(non_epi)):
            f.write("%s\n" % line)




if __name__ == "__main__":

    tools = ['MAJIQ', 'DEXSEQ', 'RMATS']

    for tool in tools:
        correlation(tool)