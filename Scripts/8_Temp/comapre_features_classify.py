import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve
from collections import defaultdict
from functools import reduce



path = '0_Files/test_neuro_H1'

def classify_v1_vs_SE_MXE():

    features_epi = pd.read_csv(f'{path}/features_nonepi_aktuell.csv', delimiter='\t')
    features_nonepi = pd.read_csv(f'{path}/features_nonepi_1.1.csv', delimiter= '\t')

    common_genes = list(set(features_epi.gene_name.values.tolist()) & set(features_nonepi.gene_name.values.tolist()))
    features_epi = features_epi[~(features_epi.gene_name.isin(common_genes))]
    features_nonepi = features_nonepi[~(features_nonepi.gene_name.isin(common_genes))]


    # sf = ['A1CF', 'ANKHD1', 'BOLL', 'BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'CELF1', 'CNOT4', 'CPEB1', 'CPEB2', 'CPEB4', 'DAZ3', 'DAZAP1', 'EIF4G2', 'ELAVL4', 'ENOX1', 'ESRP1', 'ESRP2', 'EWSR1', 'FMR1', 'FUBP1', 'FUBP3', 'FUS', 'FXR1', 'FXR2', 'G3BP2', 'HNRNPA0', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPCL1', 'HNRNPD', 'HNRNPDL', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRPLL', 'HuR', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'ILF2', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'KHSRP', 'LIN28A', 'MATR3', 'MBNL1', 'MSI1', 'NOVA1', 'NUPL2', 'PABPC1', 'PABPC3', 'PABPC4', 'PABPC5', 'PABPN1', 'PABPN1L', 'PCBP1', 'PCBP2', 'PCBP3', 'PCBP4', 'PPRC1', 'PRR3', 'PTB3', 'PTBP3', 'PUF60', 'PUM1', 'PUM2', 'QKI', 'RALY', 'RBFOX1', 'RBFOX2', 'RBFOX3', 'RBM15B', 'RBM22', 'RBM23', 'RBM24', 'RBM25', 'RBM28', 'RBM3', 'RBM38', 'RBM4', 'RBM41', 'RBM42', 'RBM45', 'RBM46', 'RBM47', 'RBM4B', 'RBM5', 'RBM6', 'RBM8A', 'RBMS1', 'RBMS2', 'RBMS3', 'RC3H1', 'SAMD4A', 'SART3', 'SF1', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF11', 'SRSF2', 'SRSF4', 'SRSF5', 'SRSF7', 'SRSF8', 'SRSF9', 'TAF15', 'TARDBP', 'TIA1', 'TRA2A', 'TRNAU1AP', 'TUT1', 'U2AF2', 'UNK', 'YBX1', 'YBX2', 'ZC3H10', 'ZC3H14', 'ZCRB1', 'ZFP36', 'ZNF326', 'ZNF638']

    sf = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    features_epi['label'] = 'epigene'
    features_nonepi['label'] = 'nonepigene'

    features_epi = features_epi.set_index('label')
    features_nonepi = features_nonepi.set_index('label')

    features_epi = features_epi[sf]
    features_nonepi = features_nonepi[sf]

    features = pd.concat([features_epi, features_nonepi], axis=0)
    features.fillna(0,inplace=True)

    sf_data = features[sf]
    features['label'] = features.index.map({'epigene': 1, 'nonepigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    clf = RandomForestClassifier(n_estimators=100, criterion='gini')
    kf = StratifiedKFold(n_splits=10)

    tprs = []
    aucs= []
    base_fpr = np.linspace(0, 1, 101)

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')
    gini_scores = []

    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        fpr, tpr, _ = roc_curve(y[test], y_score[:, 1])
        roc_auc = metrics.auc(fpr, tpr)
        aucs.append(roc_auc)

        plt.plot(fpr, tpr, 'b', alpha=0.15)
        tpr = np.interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)
        gini_scores.append(dict(zip(sf,model.feature_importances_)))

    # obtain range of AUCs
    auc_range = np.percentile(aucs, (2.5, 97.5))
    ci = float("%.2f" % (auc_range[1] - auc_range[0]))/2
    mean_auc = "%.2f" % (auc_range[1] - ci)

    # mean gini_impurity across all folds
    def foo(r, d):
        for k in d:
            r[k].append(d[k])
    gini_scores = reduce(lambda r, d: foo(r, d) or r, gini_scores, defaultdict(list))
    gini_scores = pd.DataFrame(gini_scores)
    mean_gini = gini_scores.mean(axis=0).sort_values(ascending=False)
    # print('Gini Scores: ')
    # print(mean_gini)
    mean_gini.to_csv(f'{path}/impt_features.csv', sep='\t')


    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)
    aucs = np.array(aucs)
    mean_auc = aucs.mean(axis=0)
    print('EPI vs NONEPI ', mean_auc)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b', label = 'Mean AUC = ' + str(mean_auc) + ' '+ r'$\pm$' + ' ' + str(ci))
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.title="Receiver Operating Characteristic"
    plt.show()


classify_v1_vs_SE_MXE()