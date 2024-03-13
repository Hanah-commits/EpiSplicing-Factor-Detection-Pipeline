import sys
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from venn import venn
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import Counter


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

def make_hm_plots(hm, hm_flanks, hm_pvals, hm_coeff, dir):

    # create directory to save files
    Path(f'0_Files/{dir}/{hm}/plots/').mkdir(parents=True, exist_ok=True)

    if len(hm_flanks) == 0:
        true_genes = []
        print(hm, '  ', len(true_genes))

        with open(f'0_Files/{dir}/{hm}/{hm}_truepos_epigenes.txt', 'w') as f:
            for line in true_genes:
                f.write("%s\n" % line)


    else:

        file = dir.lower()

        hm_flanks = hm_flanks.copy() # Make a copy to avoid the SettingWithCopyWarning
        genes = list(set(hm_flanks['gene_name'].values.tolist()))
        hm_flanks["type"] = hm_flanks.apply(lambda row: 'deu' if row['dPSI'] != 0 else 'non-deu', axis=1)

        # impute missing data points
        hm_flanks.fillna(0,inplace=True)

        # get absolute DEU, DHM values
        hm_flanks.loc[:, 'dPSI'] = hm_flanks['dPSI'].abs()
        hm_flanks.loc[:, hm] = hm_flanks[hm].abs()
        
        # round off
        hm_flanks.loc[:, 'dPSI'] = hm_flanks['dPSI'].round(2)
        hm_flanks.loc[:, hm] = hm_flanks[hm].round(2)


        filter_out = []
        for gene in genes:
            ## FILTER 1: remve genes whee CS exons also have DHM peakss
            gene_df = hm_flanks[hm_flanks.gene_name == gene]
            if ((gene_df[hm] != 0) & (gene_df['type'] == 'non-deu')).any():
                filter_out.append(gene)

        true_genes = [gene for gene in genes if gene not in filter_out]

        with open(f'0_Files/{dir}/{hm}/{hm}_truepos_epigenes.txt', 'w') as f:
            for line in true_genes:
                f.write("%s\n" % line)

        hm_flanks = hm_flanks[hm_flanks['gene_name'].isin(true_genes)]
        hm_flanks.to_csv(f'0_Files/{dir}/{hm}/dPSI_Mval_epi_{hm}_{file}.csv', sep='\t', index=False)

        print(hm, '  ', len(true_genes))


        # get correlation plot of true epigenes
        for gene in true_genes:
            r = hm_coeff[hm_coeff.gene_name == gene][hm].values.tolist()[0]
            p = hm_pvals[hm_pvals.gene_name == gene][hm].values.tolist()[0]
            gene_df = hm_flanks[hm_flanks['gene_name'] == gene]
            corr_plot(gene_df, gene, hm, r, p, path=f'0_Files/{dir}/{hm}/plots/')
        
    return true_genes
 

def corr_plot(gene_df, gene, hm, r, p, path):
    
    # plot deu vs dhm with regression line
    sns.regplot(x='dPSI',y=hm,data=gene_df,fit_reg=True, scatter_kws={'alpha':0.3})

    # annotate the pearson correlation coefficient text to 2 decimal places
    ax = plt.gca()
    plt.text(0.45, 0.9, 'R = {:.2f}\np = {:.4f}'.format(r,p), transform=ax.transAxes)

    # num of flanks
    num = len(gene_df)
    # define titles and axes labels
    plt.title(f'{gene}  - {num} flanks')
    plt.xlabel('|DEU|')
    plt.ylabel(f'|DHM| - {hm}')

    plt.savefig(path+f'{gene}.png', bbox_inches='tight', dpi=300)
    plt.close()


def indiv_hms(dir):

    print(f'Epigene Detection: {dir}')

    file = dir.lower()

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    # read dPSI and M-values
    try:
        flanks = pd.read_csv(f'0_Files/{dir}/DEU_DHM_{file}_flanks.tsv', delimiter='\t')
        del flanks['geneSymbol']
    except:
        print(f'No DEUs, and consequently, no DHMs available for {dir}')
        return
    
    
    # unique row index
    flanks['idx'] =  flanks['gene_name'] + flanks[['flank_start', 'flank_end']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)

    flanks_meta = flanks.copy()
    flanks.drop(['chr', 'strand'], axis=1, inplace=True)

    # FILTER 1: drop genes with less than 3 flanks
    flanks = flanks[flanks.groupby('gene_name').gene_name.transform(len) > 2]

    print('# genes > 3 exon flanks:         ', len(set(flanks.gene_name.values.tolist()))) # log

    # # # FILTER 2: genes with dPSI values but no peak -> non-epigenes
    flanks_temp = flanks.copy() # Make a copy to avoid the SettingWithCopyWarning
    flanks_temp.replace(0, None, inplace=True) # to make comparison easier in next step
    
    cols = ['gene_name'] + hms
    grouped = flanks_temp[cols].groupby('gene_name')
    epigenes = []
    non_epi = []
    for gene, group in grouped:
        if group[hms].isnull().all().all():
            non_epi.append(gene)

    filtered_flanks = flanks[~flanks['gene_name'].isin(non_epi)].copy()
    filtered_flanks.set_index('idx', inplace=True)

    print('# correlation candidates: ', len(set(filtered_flanks.gene_name.values.tolist()))) # check num before and after noepi filtering

    if len(set(filtered_flanks.gene_name.values.tolist())) > 0:

        # if one or more flanks of a gene have less than four peaks
        filtered_flanks.fillna(0, inplace=True)
        
        ## STEP 0: Find genes with strong DEU-DHM correlations

        filtered_flanks.drop(columns=['feature', 'score', 'flank_start', 'flank_end'], inplace=True) #remove unnecessary cols
        coeff = filtered_flanks.groupby('gene_name').corr(method=pearsonr_coeff)
        coeff = coeff.fillna(0)

        # internal column filtering
        coeff.to_csv(f'0_Files/{dir}/coeff.csv', sep='\t')
        coeff = pd.read_csv(f'0_Files/{dir}/coeff.csv', delimiter='\t')
        coeff = process_dataframe(coeff, hms)
    
        ## STEP 1: Obtain p values of DEU-DHM correlations

        pval = filtered_flanks.groupby('gene_name').corr(method=pearsonr_pval)

        # internal column filtering
        pval.to_csv(f'0_Files/{dir}/pvals.csv', sep='\t')
        pval = pd.read_csv(f'0_Files/{dir}/pvals.csv', delimiter='\t')
        pval = process_dataframe(pval, hms)

        ## STEP 2: Adjust the p values using Benjamini-Hochberg method
        adj_pvals = adjust_pvalue(pval)

        # STEP 3: Find epigenes: genes where adjusted_pval <= 0.05, R  >= 0.5 

        #rename and rearrange columns
        adj_pvals.columns = hms + ['gene_name']
        adj_pvals = adj_pvals[['gene_name']+ hms]

        correlated_genes = list(set(adj_pvals.gene_name.values.tolist()))
        hm_epigenes_dict = {hm: [] for hm in hms}
        for gene in correlated_genes:
            for hm in hms:
                condition = (adj_pvals.gene_name == gene) & (adj_pvals[hm] <= 0.05) & (coeff.gene_name == gene) & (coeff[hm] >= 0.5)
                if condition.any():
                    hm_epigenes_dict[hm].append(gene)

        epigenes = list(set([item for sublist in hm_epigenes_dict.values() for item in sublist]))
        print('PRE-FILTERING:  ', len(epigenes)) # temp comment

        for gene in epigenes:

            with open(f'0_Files/{dir}/{file}_epigenes.txt', 'w') as f:
                for line in list(set(epigenes)):
                    f.write("%s\n" % line)

        ## STEP 4: Make corr plots of hm-specific epigenes:
                    
        if len(epigenes) > 0:            
            # get flanks of hm-specific epigenes
            i = 0
            true_epigenes = []
            for hm in hms:

                hm_epigenes = hm_epigenes_dict[hm]
            
                # get flanks, pvals and coeffs of hm-specific epispliced genes
                hm_flanks = flanks_meta[flanks_meta['gene_name'].isin(hm_epigenes)]
                hm_coeff = coeff[coeff['gene_name'].isin(hm_epigenes)]
                hm_pvals = adj_pvals[adj_pvals['gene_name'].isin(hm_epigenes)]

                # hm-specific corrplot
                true_epigenes.append(make_hm_plots(hm, hm_flanks, hm_pvals, hm_coeff, dir))

                i+=1

            epigenes = list(set([item for items in true_epigenes for item in items]))
    
    print('Epigenes     ', len(epigenes))
    print('Non-Epigenes ', len(non_epi))
    
    # get flanks of all epispliced genes
    flanks_meta[flanks_meta['gene_name'].isin(epigenes)].to_csv(f'0_Files/{dir}/dPSI_Mval_epi_{file}.csv', sep='\t', index=False)

    # # get flanks of non-epispliced genes
    flanks_meta[flanks_meta['gene_name'].isin(non_epi)].to_csv(f'0_Files/{dir}/dPSI_Mval_nonepi_{file}.csv', sep='\t', index=False)

    ## save epi and nonepigenes
    with open(f'0_Files/{dir}/{file}_filtered_epigenes.txt', 'w') as f:
        for line in list(set(epigenes)):
            f.write("%s\n" % line)

    with open(f'0_Files/{dir}/{file}_nonepigenes.txt', 'w') as f:
        for line in list(set(non_epi)):
            f.write("%s\n" % line)


def plot_venn(dir):

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]
    tissue1 = d['tissue1'].capitalize()
    tissue2 = d['tissue2'].capitalize()

    hm_epigenes = []
    
    try:
        for hm in hms:
            with open(f'0_Files/{dir}/{hm}/{hm}_truepos_epigenes.txt') as file:
                epigenes = [line.rstrip() for line in file]
                hm_epigenes.append(epigenes)
    except:
        print(f'No epigenes available for {dir}. Cannot make a venn diagram.')
        return
    
    # visualise overlap
    info = tissue1 + ' - ' + tissue2

    unique_epi =  set(i for j in hm_epigenes for i in j)
    title = info + ' : ' + str(len(unique_epi)) + ' Epigenes'

    labels = [elem.split('_adj')[0]+ ' : ' + str(len(hm_epigenes[ind])) for ind,elem in enumerate(hms)]
    values = [set(i) for i in hm_epigenes]

    data = dict(zip(labels, values))
    venn(data, cmap="plasma")

    plt.title(title)
    plt.savefig(f'0_Files/{dir}/hm_overlap_' + info + '.png')
    plt.close()


def common_genes():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    dirs = ['DEXSEQ', 'MAJIQ', 'RMATS']

    # get overlap from 2/3 tools : epigenes
    epigenes = {}
    for dir in dirs:
            hm_epigenes = []
            try:
                for hm in hms:
                    with open(f'0_Files/{dir}/{hm}/{hm}_truepos_epigenes.txt') as file:
                        genes = [line.rstrip() for line in file]
                        hm_epigenes.extend(genes)
            except:
                print(f'No epigenes available for {dir}. Cannot find common epigenes between this tool and the other two tools.')
                continue

            epigenes[dir] = list(set(hm_epigenes))        

    # Count occurrences of genes across all three tools
    epigenee_counts = Counter(value for sublist in epigenes.values() for value in sublist)

    # Filter genes that appear in all three tools
    overlap_epigenes = [gene for gene, count in epigenee_counts.items() if count >= 2]

    with open(f'0_Files/common_epigenes.txt', 'w') as f:
        for line in list(set(overlap_epigenes)):
            f.write("%s\n" % line)


    # get overlap from all three tools : nonepigenes
    nonepigenes = {}
    for dir in dirs:
        file = dir.lower()
        try:
            with open(f'0_Files/{dir}/{file}_nonepigenes.txt') as file:
                        nonepigenes[dir]  = [line.rstrip() for line in file]
        except:
            print(f'No DEUs available for {dir}. Cannot find common non-epigenes between this tool and the other two tools.')
            continue
                    
    # Count occurrences of genes across all three tools
    nonepigene_counts = Counter(value for sublist in nonepigenes.values() for value in sublist)

    # Filter genes that appear in all three tools
    overlap_nonepigenes = [gene for gene, count in nonepigene_counts.items() if count > 2]

    with open(f'0_Files/common_nonepigenes.txt', 'w') as f:
        for line in list(set(overlap_nonepigenes)):
            f.write("%s\n" % line)

    
    # outputlog
    print('# Epigenes in 2/3 tools:    ', len(overlap_epigenes))
    print('# NonEpigenes in 3/3 tools: ', len(overlap_nonepigenes))




if __name__ == "__main__":

    tools = ['DEXSEQ', 'MAJIQ', 'RMATS']

    for tool in tools:
        # get epigenes and their correlation plots
        indiv_hms(tool)

        # venn diagram of true epigenes
        plot_venn(tool)

    common_genes()