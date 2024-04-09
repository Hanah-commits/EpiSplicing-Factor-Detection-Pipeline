import pandas as pd
import math
import json
import numpy as np
from scipy import stats

def prep():
    epi = pd.read_csv('0_Files/Post-processing/features_epi.csv', delimiter='\t')
    nonepi = pd.read_csv('0_Files/Post-processing/features_nonepi.csv', delimiter='\t')

    epi["label"] = "epi"
    nonepi["label"] = "nonepi"

    sfs = pd.read_csv('0_Files/Post-processing/impt_features.csv', delimiter='\t')['Unnamed: 0'].values.tolist()

    # RBPs with no binding site in any flank
    all_zero = []
    for column in epi:  # iterates by-name
        if epi[column].isna().all() or (epi[column] == 0).all():
            all_zero.append(column)

    for column in nonepi:  # iterates by-name
        if nonepi[column].isna().all() or (nonepi[column] == 0).all():
            all_zero.append(column)

    sfs = [sf for sf in sfs if sf not in all_zero]

    return epi, nonepi

def adjust_pvalue(df):
    pval_cols = df.columns.tolist()
    new_cols = []
    col_names = []
    for col in pval_cols:

        # get indices of null values
        na_idx = df[df[col].isnull()].index.tolist()

        # adjust non-null p values
        pvals = df[col].values.tolist()
        pvals = [x for x in pvals if not math.isnan(x)]
        adj_pval = p_adjust_bh(pvals).tolist()

        # insert null at original indices
        for idx in na_idx:
            adj_pval.insert(idx, None)

        new_cols.append(adj_pval)
        col_names.append(col)

    # adjusted p values as new df
    df1 = pd.DataFrame(columns=col_names)
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


def significane_hms(hm):

    prob_dict = {}
    # sfs = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    enriched_epi = []
    enriched_nonepi = []

    epi, nonepi = prep()

    epi = epi[epi['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    nonepi = nonepi[nonepi['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]

    for sf in sfs:
        prob = stats.ttest_ind(epi[sf], nonepi[sf], equal_var = False, alternative='greater')
        prob_dict[sf] = prob[1]
        if prob[1] < 0.05 :
            enriched_epi.append(sf)

    for sf in sfs:
        prob = stats.ttest_ind(nonepi[sf], epi[sf], equal_var = False, alternative='greater')
        prob_dict[sf] = prob[1]
        if prob[1] < 0.05 :
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
