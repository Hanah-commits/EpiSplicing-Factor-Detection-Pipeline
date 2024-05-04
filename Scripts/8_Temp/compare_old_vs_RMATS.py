import pandas as pd
import json
from pathlib import Path
import os
import sys
import numpy as np
import warnings
from scipy.stats import pearsonr


def SE(path, genes):


    # STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

    # Keep relevant columns
    file = path + '/RMATS/SE.MATS.JC.txt'
    rmats = pd.read_csv(file, delimiter='\t')
    col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
    rmats = rmats[col_list]

    print('Processing RMATS output: Skipped Exons \n')
    print('Originally: ', len(genes))
    print('# genes reported:                ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log

    # use | dPSI | and only true values
    rmats['IncLevelDifference'] = rmats['IncLevelDifference'].abs()
    rmats = rmats[rmats['FDR'] <=0.05]

    # print('FDR-adj pvalue <= 0.05:          ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log

    if len(rmats) == 0:
        print(' No skipped exons to process \n')
        sys.exit(0)

    # FILTER 1: Get AS ( |dPSI| > 0.2, FDR < 0.05)
    rmats_AS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() >= 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]

    if len(rmats_AS) == 0:
        print(' No skipped exons to process \n')
        sys.exit(0)

    # FILTER 2: If skipped exon is reported many times,  pick single dPSI score (can happen if down/upstream exons vary)

    ## get the largest dPSI value for AS exons (most differentially used score)
    df = rmats_AS.copy()
    # Create 'dPSI' and 'dPSI' columns
    df['dPSI'] = df.groupby('exonStart_0base')['IncLevelDifference'].transform(lambda x: ','.join(x.astype(str)))
    df['dPSI'] = df['dPSI'].str.split(',').apply(lambda x: max(map(float, x)) if x[0] else None)

    # Keep only rows where 'IncLevelDifference' is equal to 'dPSI'
    df = df[df['IncLevelDifference'] == df['dPSI']]

    # FILTER 3: Drop duplicate exon entries
    df = df.drop_duplicates(subset=["GeneID", "strand", "exonStart_0base", "exonEnd"], keep='first')

    # Assign the modified DataFrame back to the original variable
    rmats_AS = df

    # FILTER 4: Keep coords of single version of exon if A3SS/A5SS events exist (to prevent 2+ flanks per exon)

    # #                   GeneID geneSymbol    chr strand  IncLevelDifference       FDR  exonStart_0base   exonEnd   dPSI
    # # 4634  ENSG00000126456.15       IRF3  chr19      -               0.449  0.000202         49664442  49664673  0.449
    # # 4636  ENSG00000126456.15       IRF3  chr19      -               0.584  0.000012         49664552  49664673  0.584
    # # 4640  ENSG00000126456.15       IRF3  chr19      -               0.363  0.015466         49664586  49664673  0.363

    def A3SS_A5SS_filter(group, subset_column):
        group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
        group.drop_duplicates(subset=['exonStart_0base'], keep='first', inplace=True)
        group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
        group.drop_duplicates(subset=['exonEnd'], keep='first', inplace=True)
        return group

    rmats_AS = rmats_AS.groupby('GeneID').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
    rmats_AS.reset_index(drop=True, inplace=True)

    print('| IncLevelDifference | > 0.2:    ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log

    ## STEP 2: Prepare bedtools input

    # temp output fiilee
    df = rmats_AS.copy()
    df['feature'] = "Exon"
    df['score'] = "."
    df['exonStart_0base'] = pd.to_numeric(df['exonStart_0base']) + 1
    df[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "GeneID", "dPSI"]].to_csv(f'0_Files/RMATS/SE_exons.tsv', index=False, sep='\t', header=True)

    # save found epigenes
    geneids = list(set(df[df['geneSymbol'].isin(genes)]['GeneID'].values.tolist()))
    with open(f'0_Files/geneids.txt', 'w') as f:
        for line in list(set(geneids)):
            f.write("%s\n" % line)

    # save geensymbol-geneids of found epigenes:
    filtered_df = df[df['geneSymbol'].isin(genes)]
    gene_id_list = list(set(filtered_df['GeneID'].values.tolist()))
    gene_symbols_list = list(set(filtered_df['geneSymbol'].values.tolist()))

    # Create a dictionary mapping gene symbol to gene ID
    gene_symbol_gene_id_dict = dict(zip(gene_id_list, gene_symbols_list))


    df_temp = df.copy()
    del(df_temp['exonStart_0base'])
    del(df['exonEnd'])

    df.rename(columns={'exonStart_0base': 'exon_coord0'}, inplace = True)
    df_temp.rename(columns={'exonEnd': 'exon_coord0'}, inplace=True)

    df = pd.concat([df_temp, df]).sort_index(kind='merge')

    keep_cols = ['chr', 'exon_coord0', 'strand']
    df_bed = df[keep_cols]
    df_bed = df_bed.drop_duplicates()
    # to fit bedtools input requirements
    df_bed['exon_coord1'] = pd.to_numeric(df_bed['exon_coord0']) + 1
    df_bed['feature'] = "flank"
    df_bed['score'] = "."


    df_bed = df_bed[['chr', "exon_coord0", "exon_coord1", "feature", "score", "strand"]]
    df_bed.to_csv(f'0_Files/RMATS/SE.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
    df.to_csv(f'0_Files/RMATS/SE_exons.csv', index=False, sep='\t', header=True)

    return gene_symbol_gene_id_dict

    
def MXE(path, genes):


    # Keep relevant columns
    file = path+ '/RMATS/MXE.MATS.JC.txt'
    rmats = pd.read_csv(file, delimiter='\t')

    col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd']
    rmats = rmats[col_list]

    print('Processing RMATS output: Mutually Exclusive Exons \n')
    print('Originally: ', len(genes))
    print('# genes reported:                ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log

    # use | dPSI | and only true values
    rmats['IncLevelDifference'] = rmats['IncLevelDifference'].abs()
    rmats = rmats[rmats['FDR'] <=0.05]

    # print('FDR-adj pvalue <= 0.05:          ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log

    if len(rmats) == 0:
        print(' No mutually exclusive exons to process \n')
        sys.exit(0)


    # STEP 2 : Split into multiple rows, keeping one exon coord in one row.

    # Create two DataFrames, one for each row
    row1 = rmats[['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', '1stExonStart_0base', '1stExonEnd']].copy()
    row2 = rmats[['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', '2ndExonStart_0base', '2ndExonEnd']].copy()

    # Rename columns 
    row1.columns = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
    row2.columns = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']

    # mark exon order
    row1['exon_order'] = 1
    row2['exon_order'] = 2

    # Concatenate the two DataFrames to get the output
    rmats = pd.concat([row1, row2], ignore_index=True)

    # # housekeeping
    # os.system('rm 0_Files/RMATS/rmats_*.bed')

    # STEP 3: Get dPSI scores based on inclusion exon
    # NOTE: the inclusion isoform includes the exon that is “earlier” in the transcript.

    rmats['dPSI'] = np.where(
        (rmats['strand'] == '+') & (rmats['exon_order'] == 1),
        rmats['IncLevelDifference'],
        np.where(
            (rmats['strand'] == '+') & (rmats['exon_order'] == 2),
            1 - rmats['IncLevelDifference'],
            np.where(
                (rmats['strand'] == '-') & (rmats['exon_order'] == 2),
                rmats['IncLevelDifference'],
                np.where(
                    (rmats['strand'] == '-') & (rmats['exon_order'] == 1),
                    1 - rmats['IncLevelDifference'],
                    np.nan  # default value for other cases
                )
            )
        )
    )

    # FILTER 1: Get true MXE events
    try:
        SE_exons = pd.read_csv("0_Files/RMATS/SE_exons.csv", delimiter='\t')
        SE_exons = list(set(SE_exons[SE_exons.dPSI > 0.2].exon_coord0.values.tolist()))
        rmats = rmats[(~rmats['exonStart_0base'].isin(SE_exons)) & (~rmats['exonEnd'].isin(SE_exons))] # covers A3SS,A5SS versions of skipped exons

        # print('True MXE:                        ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log
    except:
        print('No skipped exons available: All are considered True MXE')

    # FILTER 2: Get AS ( |dPSI| > 0.2, FDR < 0.05) and CS exons ( |dPSI| < 0.2, FDR < 0.05)
    rmats_AS = rmats[(pd.to_numeric(rmats['dPSI'] >= 0.2)) & (pd.to_numeric(rmats['FDR']) <= 0.05)]

    # FILTER 3: Keep coords of single version of exon if A3SS/A5SS events exist (to prevent 2+ flanks per exon)

    # #                   GeneID geneSymbol    chr strand  IncLevelDifference       FDR  exonStart_0base   exonEnd   dPSI
    # # 4634  ENSG00000126456.15       IRF3  chr19      -               0.449  0.000202         49664442  49664673  0.449
    # # 4636  ENSG00000126456.15       IRF3  chr19      -               0.584  0.000012         49664552  49664673  0.584
    # # 4640  ENSG00000126456.15       IRF3  chr19      -               0.363  0.015466         49664586  49664673  0.363

    def A3SS_A5SS_filter(group, subset_column):
        group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
        group.drop_duplicates(subset=['exonStart_0base'], keep='first', inplace=True)
        group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
        group.drop_duplicates(subset=['exonEnd'], keep='first', inplace=True)
        return group

    rmats_AS = rmats_AS.groupby('GeneID').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
    rmats_AS.reset_index(drop=True, inplace=True)

    print('| IncLevelDifference | > 0.2:    ', len(list(set(rmats.geneSymbol.values.tolist()) & set(genes)))) # log

    df = rmats_AS.copy()
    # temp output fiilee
    df['feature'] = "Exon"
    df['score'] = "."
    df['exonStart_0base'] = pd.to_numeric(df['exonStart_0base']) + 1
    df[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "GeneID", "dPSI"]].to_csv(f'0_Files/RMATS/MXE_exons.tsv', index=False, sep='\t', header=True)
    df_temp = df.copy()
    del(df_temp['exonStart_0base'])
    del(df['exonEnd'])

    # save found epigenes
    geneids = list(set(df[df['geneSymbol'].isin(genes)]['GeneID'].values.tolist()))
    genesymbols = list(set(df[df['geneSymbol'].isin(genes)]['geneSymbol'].values.tolist()))
    with open(f'0_Files/geneids.txt', 'a') as f:
        for line in list(set(geneids)):
            f.write("%s\n" % line)

    # Create a dictionary mapping gene symbol to gene ID
    gene_symbol_gene_id_dict = dict(zip(geneids, genesymbols))

    df.rename(columns={'exonStart_0base': 'exon_coord0'}, inplace = True)
    df_temp.rename(columns={'exonEnd': 'exon_coord0'}, inplace=True)

    df = pd.concat([df_temp, df]).sort_index(kind='merge')

    keep_cols = ['chr', 'exon_coord0', 'strand']
    df_bed = df[keep_cols]
    df_bed = df_bed.drop_duplicates()
    # to fit bedtools input requirements
    df_bed['exon_coord1'] = pd.to_numeric(df_bed['exon_coord0']) + 1
    df_bed['feature'] = "flank"
    df_bed['score'] = "."


    df_bed = df_bed[['chr', "exon_coord0", "exon_coord1", "feature", "score", "strand"]]
    df_bed.to_csv(f'0_Files/RMATS/MXE.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
    df.to_csv(f'0_Files/RMATS/MXE_exons.csv', index=False, sep='\t', header=True)

    return gene_symbol_gene_id_dict


def SE_MXE(geneids):

    with open('paths.json') as f:
        d = json.load(f)

    fasta = d['Reference fasta']
    ref_genome= fasta+".fai"

    e = 0
    try:
        SE = pd.read_csv(f'0_Files/RMATS/SE_exons.tsv', delimiter='\t', names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI" ], skiprows=1)
    except:
        SE = pd.DataFrame()
        e +=1

    try:
        MXE = pd.read_csv(f'0_Files/RMATS/MXE_exons.tsv', delimiter='\t',  names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI"], skiprows=1)
    except:
        MXE = pd.DataFrame()
        e +=1

    if e == 2:
        print('No RMATS exons available \n')
        sys.exit()


    def A3SS_A5SS_filter(group, subset_column):
        group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
        # Drop duplicates based on 'exonStart_0base' and keep the one with the largest 'dPSI'
        group = group.groupby('exonStart_0base').apply(lambda x: x.loc[x[subset_column].idxmax()])
        

        group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
        # Drop duplicates based on 'exonEnd' and keep the one with the largest 'dPSI'
        group = group.groupby('exonEnd').apply(lambda x: x.loc[x[subset_column].idxmax()])
        
        return group


    ## STEP 1: get combined dPSI scores of AS exons

    if len(SE) > 0 and len(MXE) >0 : 
        SE_MXE_exons =  pd.concat([SE, MXE], ignore_index=True)
    elif len(SE) > 0:
        SE_MXE_exons = SE
    else:
        SE_MXE_exons = MXE

    ## FILTER 1: A3SS/A5SS -  Avoid many A3SS/A5SS versions of AS exons
    SE_MXE_exons = SE_MXE_exons.groupby('geneSymbol').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
    SE_MXE_exons = SE_MXE_exons.reset_index(drop=True)

    ## FILTER 2: Keep only AS exons
    SE_MXE_exons = SE_MXE_exons[SE_MXE_exons.dPSI > 0.2]

    print('# genes with SE and/or MXE exons:   ', len(list(set(SE_MXE_exons.geneSymbol.values.tolist()) & set(geneids)))) # log

    ## STEP 2: Get exon flanks

    SE_MXE_exons[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI"]].to_csv('0_Files/RMATS/rmats_exons_coords.bed', index=False, sep='\t', header=False)

    # exon boundary external flanks
    os.system("bedtools flank -i 0_Files/RMATS/rmats_exons_coords.bed -g " + ref_genome + " -b 200 > 0_Files/flanks.bed" )

    # separate start,stop flank coords
    os.system("sed -n 'n;p' 0_Files/flanks.bed > 0_Files/stop.bed")
    os.system("sed -n 'p;n' 0_Files/flanks.bed > 0_Files/start.bed")

    # exon boundary internal flanks
    os.system("bedtools slop -i 0_Files/start.bed -g " + ref_genome + " -l 0 -r 200 > 0_Files/start_flanks.bed")
    os.system("bedtools slop -i 0_Files/stop.bed -g " + ref_genome +" -l 200 -r 0 > 0_Files/stop_flanks.bed")

    # combine start,stop flank coords
    os.system("paste -d'\n' 0_Files/start_flanks.bed 0_Files/stop_flanks.bed | sort -k1,1 -k2,2n > 0_Files/RMATS/rmats_flanks200.bed")

    # remove intermediate files
    os.system("rm 0_Files/start*.bed")
    os.system("rm  0_Files/stop*.bed")
    os.system("rm 0_Files/flanks.bed")


    ## FILTER 3: Drop flanked AS exns overlapping with TSS regions. CS exons are TSS-free since exon_coords.bed alreeady has TSS-filtered exons

    os.system('bedtools intersect -wa -a 0_Files/RMATS/rmats_flanks200.bed -b 0_Files/TSS.bed -s -v > 0_Files/rmats_flanks200_temp.bed && mv 0_Files/rmats_flanks200_temp.bed 0_Files/RMATS/rmats_flanks200.bed')

    #%# Note: Exons have varying lengths. Flanks can overlap.


def MANorm(genes):

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    try:
        deu_flanks = pd.read_csv('0_Files/RMATS/rmats_flanks200.bed', delimiter='\t', header=None)
    except:
        print('No RMATS exons to annotate \n \n')
        sys.exit()

    deu_flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'geneSymbol', 'dPSI']
    deu_flanks.drop(columns=['feature', 'score'], inplace=True)
    deu_flanks.drop_duplicates(inplace=True)

    print ('Annotating RMATS exons with HM peaks \n')
    print('TSS Filtering:                   ', len(list(set(deu_flanks.geneSymbol.values.tolist()) & set(genes)))) # log

    dhm_flanks = pd.read_csv('0_Files/MANorm/DHM_peaks_annotation.tsv', delimiter='\t')

    ## STEP 1: combine DEU and DHM scores
    dhm_flanks = dhm_flanks.merge(deu_flanks, on=['chr', 'flank_start', 'flank_end', 'strand', 'geneSymbol'], how='left')

    # Fill NaN values with 0 
    dhm_flanks['dPSI'] = dhm_flanks['dPSI'].fillna(0)

    ## STEP 2:  keep only genes where DEU scores are available
    deu_genes = list(set(deu_flanks.geneSymbol.values.tolist()))
    dhm_flanks = dhm_flanks[dhm_flanks.geneSymbol.isin(deu_genes)]

    print('TSL Filtering:                   ', len(list(set(dhm_flanks.geneSymbol.values.tolist()) & set(genes))), '\n \n') # log

    return list(set(dhm_flanks.geneSymbol.values.tolist()) & set(genes))


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
        # print(hm, '  ', len(true_genes))

        # with open(f'0_Files/{dir}/{hm}/{hm}_truepos_epigenes.txt', 'w') as f:
        #     for line in true_genes:
        #         f.write("%s\n" % line)


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

        # with open(f'0_Files/{dir}/{hm}/{hm}_truepos_epigenes.txt', 'w') as f:
        #     for line in true_genes:
        #         f.write("%s\n" % line)

        hm_flanks = hm_flanks[hm_flanks['gene_name'].isin(true_genes)]
        # hm_flanks.to_csv(f'0_Files/{dir}/{hm}/dPSI_Mval_epi_{hm}_{file}.csv', sep='\t', index=False)

        # print(hm, '  ', len(true_genes))


        # get correlation plot of true epigenes
        for gene in true_genes:
            r = hm_coeff[hm_coeff.gene_name == gene][hm].values.tolist()[0]
            p = hm_pvals[hm_pvals.gene_name == gene][hm].values.tolist()[0]
            gene_df = hm_flanks[hm_flanks['gene_name'] == gene]
        
    return true_genes
 

def indiv_hms(ogdir, newdir, genes, nonepi_og, check_genes, i, dir = "RMATS"):

    file = dir.lower()

    with open(f'{ogdir}/../paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    # read dPSI and M-values
    try:
        flanks = pd.read_csv(f'{newdir}/0_Files/{dir}/DEU_DHM_{file}_flanks.tsv', delimiter='\t')
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

    # print('# genes > 3 exon flanks:         ', len(set(flanks.gene_name.values.tolist()))) # log

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

    # print('# correlation candidates:        ', len(set(filtered_flanks.gene_name.values.tolist()))) # check num before and after noepi filtering

    if len(set(filtered_flanks.gene_name.values.tolist())) > 0:

        # if one or more flanks of a gene have less than four peaks
        filtered_flanks.fillna(0, inplace=True)
        
        ## STEP 0: Find genes with strong DEU-DHM correlations

        filtered_flanks.drop(columns=['feature', 'score', 'flank_start', 'flank_end'], inplace=True) #remove unnecessary cols
        warnings.filterwarnings("ignore", message="An input array is constant; the correlation coefficient is not defined.") # Suppress warning when corr is NaN
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
        print('Epigenes detcted before filtering: ', len(epigenes))
        print('Overlap before filtering:                   ', len(list(set(epigenes) & set(genes)))) # temp comment

        print('\n\n')
        filtered_out = [g for g in check_genes if g not in epigenes]
        flanks[flanks.gene_name.isin(filtered_out)].to_csv(f'0_Files/dropped_epigenes_{i}.tsv', sep='\t', index=False)
        print('\n\n')

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
    
    print('Epigenes detcted after filtering: ', len(epigenes))
    print('Overlap after filtering:                         ', len(list(set(epigenes)& set(genes))))
    print('Overlap Non-Epigenes                     ', len(list(set(non_epi) & set(nonepi_og))), '\n\n')
    
    # # get flanks of all epispliced genes
    # flanks_meta[flanks_meta['gene_name'].isin(epigenes)].to_csv(f'0_Files/{dir}/dPSI_Mval_epi_{file}.csv', sep='\t', index=False)

    # # # get flanks of non-epispliced genes
    # flanks_meta[flanks_meta['gene_name'].isin(non_epi)].to_csv(f'0_Files/{dir}/dPSI_Mval_nonepi_{file}.csv', sep='\t', index=False)

    # ## save epi and nonepigenes
    # with open(f'0_Files/{dir}/{file}_filtered_epigenes.txt', 'w') as f:
    #     for line in list(set(epigenes)):
    #         f.write("%s\n" % line)

    # with open(f'0_Files/{dir}/{file}_nonepigenes.txt', 'w') as f:
    #     for line in list(set(non_epi)):
    #         f.write("%s\n" % line)



# STEP 1: Get features of DEU flanks of dropped genes
old_dirs = ['/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/ectodermalcell_mesodermalcell_1685435254.5447733/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_ectodermalcell_1685402322.3652477/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/neuronalcell_H1_1685391226.3440962/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/mesodermalcell_neuronalcell_1685428884.2537532/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/mesodermalcell_H1_1685398815.2181895/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_mesodermalcell_1685404937.8074098/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_neuronalcell_1685416269.2483714/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/ectodermalcell_neuronalcell_1685413355.7715614/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_H1_1685422988.3170376/0_Files', '/home/hanah/EpiSplicing-Factor-Detection-Pipeline/Output/ectodermalcell_H1_1685437523.1683042/0_Files']
new_dirs = ['/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/ectodermalcell_mesodermalcell_1685435254.5447733', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_ectodermalcell_1685402322.3652477', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/neuronalcell_H1_1685391226.3440962', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/mesodermalcell_neuronalcell_1685428884.2537532', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/mesodermalcell_H1_1685398815.2181895', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_mesodermalcell_1685404937.8074098', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_neuronalcell_1685416269.2483714', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/ectodermalcell_neuronalcell_1685413355.7715614', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/endodermalcell_H1_1685422988.3170376', '/home/hanah/Epi_new/EpiSplicing-Factor-Detection-Pipeline/Output/ectodermalcell_H1_1685437523.1683042']

i = 0
epigenes_dropped = []
non_epigenes_dropped = []
for dir in old_dirs:

    print(f'\n {dir} \n')


    with open(f'{dir}/epigenes.txt', 'r') as file:
        epigenes = list(set([line.strip() for line in file.readlines()]))

    with open(f'{dir}/nonepigenes.txt', 'r') as file:
        non_epigenes = list(set([line.strip() for line in file.readlines()]))

    print('\n EPIGENES \n')

    geneids_dict_1 = SE(new_dirs[i],epigenes)
    print('\n')
    geneids_dict_2 = MXE(new_dirs[i],epigenes)
    print('\n')

    with open('0_Files/geneids.txt', 'r') as file:
        geneids = list(set([line.strip() for line in file.readlines()]))

    ## TODO: Get genes that were dropped as a result of ONLY the peak-leak filter (DEU_DHM)

    ## get ensembl version of these genes
    SE_MXE(geneids)
    print('\n')
    filtered_ids = MANorm(geneids)
    print('\n')

    # get genesymbol of overlapping epigenes from previous step
    geneids_dict = {**geneids_dict_1, **geneids_dict_2}
    check_genes = [geneids_dict[id] for id in filtered_ids]
    indiv_hms(dir, new_dirs[i], epigenes, non_epigenes, check_genes, i)

    print('\n')
    
    print ('#################')

    i += 1
