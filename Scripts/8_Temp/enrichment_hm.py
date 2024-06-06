import pandas as pd
import math
import json
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def prep():
    epi = pd.read_csv('0_Files/Post-processing/features_epi.csv', delimiter='\t')
    nonepi = pd.read_csv('0_Files/Post-processing/features_nonepi.csv', delimiter='\t')

    epi["label"] = "epi"
    nonepi["label"] = "nonepi"

    # remove genes with both labels
    common_genes = list(set(epi.gene_name.values.tolist()) & set(nonepi.gene_name.values.tolist()))
    epi = epi[~(epi.gene_name.isin(common_genes))]
    nonepi = nonepi[~(nonepi.gene_name.isin(common_genes))]

    sfs = pd.read_csv('0_Files/Post-processing/impt_features.csv', delimiter='\t')['Unnamed: 0'].values.tolist()

    # keep only strong binding events
    for df in [epi, nonepi]:
        df.loc[:, sfs] = df.loc[:, sfs].applymap(lambda val: 0 if val < 2 else val)

    # RBPs with no binding site in any flank
    all_zero = []
    for column in epi:  # iterates by-name
        if epi[column].isna().all() or (epi[column] == 0).all():
            all_zero.append(column)

    for column in nonepi:  # iterates by-name
        if nonepi[column].isna().all() or (nonepi[column] == 0).all():
            all_zero.append(column)

    sfs = [sf for sf in sfs if sf not in all_zero]

    return epi, nonepi, sfs


def significane_hms(hm):

    prob_dict = {}
    # sfs = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    enriched_epi = []
    enriched_nonepi = []

    epi, nonepi, sfs = prep()

    epi = epi[epi['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    nonepi = nonepi[nonepi['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]

    for sf in sfs:
        prob = stats.mannwhitneyu(epi[sf], nonepi[sf], alternative='greater')
        prob_dict[sf] = prob[1]
    
    # adjust pvalues
    adjusted_pvalues = multipletests(list(prob_dict.values()), method='fdr_bh')[1]

    # Filter features based on adjusted p-values
    for sf, adj_pval in zip(prob_dict.keys(), adjusted_pvalues):
        if adj_pval < 0.05:
            enriched_epi.append(sf)

    prob_dict = {}
    for sf in sfs:
        prob = stats.mannwhitneyu(nonepi[sf], epi[sf], alternative='greater')
        prob_dict[sf] = prob[1]

    # adjust pvalues
    adjusted_pvalues = multipletests(list(prob_dict.values()), method='fdr_bh')[1]

    # Filter features based on adjusted p-values
    for sf, adj_pval in zip(prob_dict.keys(), adjusted_pvalues):
        if adj_pval < 0.05:
            enriched_nonepi.append(sf)


    print('RBPS enriched in epigene flanks: ' ,len(enriched_epi))
    print('RBPs enriched in non-epigene flanks ', len(enriched_nonepi))
    print('###############################')

    with open(f'0_Files/Post-processing/enriched_epi_{hm}.txt', 'w') as f:
        for rbp in enriched_epi:
            f.write(f"{rbp}\n")

    with open(f'0_Files/Post-processing/enriched_nonepi_{hm}.txt', 'w') as f:
        for rbp in enriched_nonepi:
            f.write(f"{rbp}\n")
            
    return enriched_epi, enriched_nonepi


if __name__ == "__main__":


    # HM-specific epi vs nonepi
    with open('paths.json') as f:
            d = json.load(f)

    hms = d["Histone modifications"]

    for hm in hms:
        significane_hms(hm)
