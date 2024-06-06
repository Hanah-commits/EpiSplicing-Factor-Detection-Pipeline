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

    sfs = ['A1CF', 'ANKHD1', 'BOLL', 'BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'CELF1', 'CNOT4', 'CPEB1', 'CPEB2', 'CPEB4', 'DAZ3', 'DAZAP1', 'EIF4G2', 'ELAVL4', 'ENOX1', 'ESRP1', 'ESRP2', 'EWSR1', 'FMR1', 'FUBP1', 'FUBP3', 'FUS', 'FXR1', 'FXR2', 'G3BP2', 'HNRNPA0', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPCL1', 'HNRNPD', 'HNRNPDL', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRPLL', 'HuR', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'ILF2', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'KHSRP', 'LIN28A', 'MATR3', 'MBNL1', 'MSI1', 'NOVA1', 'NUPL2', 'PABPC1', 'PABPC3', 'PABPC4', 'PABPC5', 'PABPN1', 'PABPN1L', 'PCBP1', 'PCBP2', 'PCBP3', 'PCBP4', 'PPRC1', 'PRR3', 'PTB3', 'PTBP3', 'PUF60', 'PUM1', 'PUM2', 'QKI', 'RALY', 'RBFOX1', 'RBFOX2', 'RBFOX3', 'RBM15B', 'RBM22', 'RBM23', 'RBM24', 'RBM25', 'RBM28', 'RBM3', 'RBM38', 'RBM4', 'RBM41', 'RBM42', 'RBM45', 'RBM46', 'RBM47', 'RBM4B', 'RBM5', 'RBM6', 'RBM8A', 'RBMS1', 'RBMS2', 'RBMS3', 'RC3H1', 'SAMD4A', 'SART3', 'SF1', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF11', 'SRSF2', 'SRSF4', 'SRSF5', 'SRSF7', 'SRSF8', 'SRSF9', 'TAF15', 'TARDBP', 'TIA1', 'TRA2A', 'TRNAU1AP', 'TUT1', 'U2AF2', 'UNK', 'YBX1', 'YBX2', 'ZC3H10', 'ZC3H14', 'ZCRB1', 'ZFP36', 'ZNF326', 'ZNF638']

    # # keep only strong binding events
    # for df in [epi, nonepi]:
    #     df.loc[:, sfs] = df.loc[:, sfs].applymap(lambda val: 0 if val < 2 else val)

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
    hms = [ "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]

    for hm in hms:
        print(hm)
        significane_hms(hm)
