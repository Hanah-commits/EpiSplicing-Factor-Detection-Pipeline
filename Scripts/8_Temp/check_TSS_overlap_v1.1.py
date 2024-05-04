import pandas as pd
import numpy as np
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve
from collections import defaultdict
from functools import reduce
import seaborn as sns


# awk -F'[:-]' -v OFS='\t' '{print $1, $2, $3, "Flank", ".", $4}' rbp_input_nonepi1.csv | awk -F'\t' -v OFS='\t' '$NF == "" {$NF = "-"} {print}' > nonepi.tsv # get flanks of epi and non-pi features
# bedtools intersect -wa -s -v -a epi.tsv -b TSS_exons.bed | sort | uniq > filtered_epi.tsv
# bedtools intersect -wa -s -a epi.tsv -b TSS_exons.bed | sort | uniq > TSS_epi.tsv

path = '0_Files/test_neuro_H1'

# TODO: check mode and extract features based on that
mode = sys.argv[1]

features_dfs = []
if mode == 'TSS11': # TSS overlap flanks epi & nonepi

    for type in ['epi', 'nonepi']:

        flanks = pd.read_csv(f'{path}/TSS_{type}.tsv', delimiter='\t', header=None)
        flanks.columns = ['chr', 'start', 'stop', 'feature', 'score', 'strand']
        features = pd.read_csv(f'{path}/features_{type}.csv', delimiter='\t')
        features[['start', 'stop']] = features['flanks'].str.split('-', n=1, expand=True)
        for col in ['start', 'stop']:
            features[col] = features[col].astype('int64')

        # TODO: check mode and extract features based on that

        # Assuming df1 and df2 are your DataFrames
        merged_df = features.merge(flanks, on=['start', 'stop'], how='left')

        # Filter out rows where values are not present in flanks
        filtered_df = merged_df[merged_df['chr'].notna() & merged_df['feature'].notna()]

        # Drop the columns from flanks that were added during the merge
        filtered_df.drop(['start', 'stop'], axis=1, inplace=True)

        filtered_df = filtered_df[[col for col in features.columns if col not in ['start', 'stop']]]

        filtered_df['label'] = type + 'gene'
        features_dfs.append(filtered_df)

elif mode == 'TSS00': # non-TSS overlap flanks epi & nonepi

    for type in ['epi', 'nonepi']:

        flanks = pd.read_csv(f'{path}/TSS_{type}.tsv', delimiter='\t', header=None)
        flanks.columns = ['chr', 'start', 'stop', 'feature', 'score', 'strand']
        features = pd.read_csv(f'{path}/features_{type}.csv', delimiter='\t')
        features[['start', 'stop']] = features['flanks'].str.split('-', n=1, expand=True)
        for col in ['start', 'stop']:
            features[col] = features[col].astype('int64')

        # TODO: check mode and extract features based on that

        # Assuming df1 and df2 are your DataFrames
        merged_df = features.merge(flanks, on=['start', 'stop'], how='left')

        # Filter out rows where values are not present in flanks
        filtered_df = merged_df[merged_df['chr'].notna() & merged_df['feature'].notna()]

        # Drop the columns from flanks that were added during the merge
        filtered_df.drop(['start', 'stop'], axis=1, inplace=True)

        filtered_df = filtered_df[[col for col in features.columns if col not in ['start', 'stop']]]

        filtered_df['label'] = type + 'gene'
        features_dfs.append(filtered_df)

elif mode == 'TSS10': # TSS overlap epi & non-TSS overlap nonepi

    suffices =  ['TSS', 'filtered']

    i = 0
    for type in ['epi', 'nonepi']:

        suffix = suffices[i]
        flanks = pd.read_csv(f'{path}/{suffix}_{type}.tsv', delimiter='\t', header=None)
        flanks.columns = ['chr', 'start', 'stop', 'feature', 'score', 'strand']
        features = pd.read_csv(f'{path}/features_{type}.csv', delimiter='\t')
        features[['start', 'stop']] = features['flanks'].str.split('-', n=1, expand=True)
        for col in ['start', 'stop']:
            features[col] = features[col].astype('int64')

        # TODO: check mode and extract features based on that

        # Assuming df1 and df2 are your DataFrames
        merged_df = features.merge(flanks, on=['start', 'stop'], how='left')

        # Filter out rows where values are not present in flanks
        filtered_df = merged_df[merged_df['chr'].notna() & merged_df['feature'].notna()]

        # Drop the columns from flanks that were added during the merge
        filtered_df.drop(['start', 'stop'], axis=1, inplace=True)

        filtered_df = filtered_df[[col for col in features.columns if col not in ['start', 'stop']]]

        filtered_df['label'] = type + 'gene'

        print(filtered_df)

        i += 1

elif mode == 'TSS01': # non-TSS overlap epi & TSS overlap epi

    suffices =  ['filtered', 'TSS']

    i = 0
    for type in ['epi', 'nonepi']:

        suffix = suffices[i]
        flanks = pd.read_csv(f'{path}/{suffix}_{type}.tsv', delimiter='\t', header=None)
        flanks.columns = ['chr', 'start', 'stop', 'feature', 'score', 'strand']
        features = pd.read_csv(f'{path}/features_{type}.csv', delimiter='\t')
        features[['start', 'stop']] = features['flanks'].str.split('-', n=1, expand=True)
        for col in ['start', 'stop']:
            features[col] = features[col].astype('int64')

        # TODO: check mode and extract features based on that

        # Assuming df1 and df2 are your DataFrames
        merged_df = features.merge(flanks, on=['start', 'stop'], how='left')

        # Filter out rows where values are not present in flanks
        filtered_df = merged_df[merged_df['chr'].notna() & merged_df['feature'].notna()]

        # Drop the columns from flanks that were added during the merge
        filtered_df.drop(['start', 'stop'], axis=1, inplace=True)

        filtered_df = filtered_df[[col for col in features.columns if col not in ['start', 'stop']]]

        filtered_df['label'] = type + 'gene'
        features_dfs.append(filtered_df)

        i += 1


def classify(): 
    features = pd.concat(features_dfs,axis=0,sort=False)
    # TODO: save with mode

    features.fillna(0,inplace=True)


    sf = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    sf_data = features[sf]
    features['label'] = features['label'].map({'epigene': 1, 'nonepigene': 0}).astype(int)
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


def heatmap():

    features = pd.concat(features_dfs,axis=0,sort=False)
    features.fillna(0,inplace=True)
    features = features.set_index('label')

    ## EPI-ENRICHED RBPS
    with open(f'{path}/enriched_epi.txt', 'r') as file:
        sfs = [line.strip() for line in file.readlines()]

    #     # sfs = [ val for val in list(features.columns) if val not in ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene', 'type']]

    features.loc[:, sfs] = features.loc[:, sfs].applymap(lambda val: 0 if val < 3 else val)

    plt.close('all')
    # plt.figure(figsize=(10, 8))

    label_colors = dict(zip(['epigene', 'nonepigene'], ["green", "grey"]))
    row_colors = features.index.map(label_colors)
    cluster = sns.clustermap(features[sfs], annot=False, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='mako',  cbar_pos=(0.9, .2, .03, .4))
    heatmap = cluster.ax_heatmap
    cbar = heatmap.collections[0].colorbar # custom y tick colorbar

    # Remove y-axis ticks and tick labels
    heatmap.set_yticks([])
    heatmap.set_yticklabels([])

    # Set x-axis tick labels
    heatmap.set_xticks(range(len(sfs)))
    heatmap.set_xticklabels(sfs, size=5, rotation=90)

    # add title and axes labels
    heatmap.set_title(f'Epi-enriched RBPs')
    heatmap.set_ylabel('')

    # Position the legend next to the plot
    plt.show()


    ## NONEPI ENRICHED RBPS
    with open(f'{path}/enriched_nonepi.txt', 'r') as file:
        sfs = [line.strip() for line in file.readlines()]

    plt.close('all')

    label_colors = dict(zip(['epigene', 'nonepigene'], ["green", "grey"]))
    row_colors = features.index.map(label_colors)
    cluster = sns.clustermap(features[sfs], annot=False, linewidth=0, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap='mako',  cbar_pos=(0.9, .2, .03, .4))
    heatmap = cluster.ax_heatmap
    cbar = heatmap.collections[0].colorbar # custom y tick colorbar

    # Remove y-axis ticks and tick labels
    heatmap.set_yticks([])
    heatmap.set_yticklabels([])

    # Set x-axis tick labels
    heatmap.set_xticks(range(len(sfs)))
    heatmap.set_xticklabels(sfs, size=5, rotation=90)

    # add title and axes labels
    heatmap.set_title(f'Nonepi-enriched RBPs')
    heatmap.set_ylabel('')

    # Position the legend next to the plot
    plt.show()
