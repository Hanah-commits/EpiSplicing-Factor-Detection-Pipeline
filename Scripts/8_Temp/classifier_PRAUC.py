import pandas as pd
import numpy as np
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, auc
from collections import defaultdict
from functools import reduce
from sklearn.utils import shuffle
import json


def stratified_classifier(output_dir):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features = features.drop('type', axis=1) # drop hm info
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)


    sf = [val for val in features.columns if val != 'label']

    sf_data = features[sf]
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    # Initialize classifier and cross-validation
    clf = RandomForestClassifier(n_estimators=100, max_features= "sqrt", class_weight='balanced', n_jobs = -1, random_state=0)
    kf = StratifiedKFold(n_splits=10)

    # Initialize arrays to store PR-AUC values
    pr_aucs = []

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')

    # Loop over folds
    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        precision, recall, _ = precision_recall_curve(y[test], y_score[:, 1])
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)
        plt.plot(recall, precision, 'b', alpha=0.15)

    # Calculate mean PR-AUC and confidence interval
    pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
    pr_auc_ci = float("%.2f" % (pr_auc_range[1] - pr_auc_range[0])) / 2
    mean_pr_auc = "%.2f" % (pr_auc_range[1] - pr_auc_ci)
    print('EPI vs NONEPI PR-AUC:', mean_pr_auc)


    # Calculate mean PR-AUC
    pr_auc_mean = np.mean(pr_aucs)
    pr_auc_std = np.std(pr_aucs)
    ci = pr_auc_std / 2

    # Plot single point for mean PR-AUC
    plt.plot(1, pr_auc_mean, 'bo', label='Mean PR-AUC = {:.2f} $\pm$ {:.2f}'.format(pr_auc_mean, ci))

    # Set plot properties
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve')
    plt.legend(loc='lower right')

    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.savefig(output_dir + 'PR_curve.png')

    
def stratified_hms_classifier(output_dir, hm):

    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info

    if len(features) == 0:
        return

    sf = [val for val in features.columns if val != 'label']

    # keep only strong binding events
    sf_data = features[sf]
    sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    # Initialize classifier and cross-validation
    clf = RandomForestClassifier(n_estimators=100, max_features= "sqrt", class_weight='balanced', n_jobs = -1, random_state=0)
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    # Initialize arrays to store PR-AUC values
    pr_aucs = []

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')

    # Loop over folds
    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        precision, recall, _ = precision_recall_curve(y[test], y_score[:, 1])
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)
        plt.plot(recall, precision, color_dict[hm])

    # Calculate mean PR-AUC and confidence interval
    pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
    pr_auc_ci = float("%.2f" % (pr_auc_range[1] - pr_auc_range[0])) / 2
    mean_pr_auc = "%.2f" % (pr_auc_range[1] - pr_auc_ci)
    print(f'EPI vs NONEPI - {hm}:', mean_pr_auc)

    # Calculate mean PR-AUC
    pr_auc_mean = np.mean(pr_aucs)
    pr_auc_std = np.std(pr_aucs)
    ci = pr_auc_std / 2

    # Plot single point for mean PR-AUC
    plt.plot(1, pr_auc_mean, 'o', color=color_dict[hm], label='Mean PR-AUC = {:.2f} $\pm$ {:.2f}'.format(pr_auc_mean, ci))

    # Set plot properties
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve - {hm}')
    plt.legend(loc='lower right')

    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.savefig(output_dir + f'{hm}_PR_curve.png')


def stratified_classifier_2(output_dir):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features = features.drop('type', axis=1) # drop hm info
    features.fillna(0,inplace=True)


    sf = [val for val in features.columns if val != 'label']

    sf_data = features[sf]
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    clf = RandomForestClassifier(n_estimators=100, criterion='gini')
    kf = StratifiedKFold(n_splits=10)

    # Initialize arrays to store PR-AUC values
    pr_aucs = []
    precisions = []
    recalls = []

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')

    # Loop over folds
    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        precision, recall, _ = precision_recall_curve(y[test], y_score[:, 1])
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)
        precisions.append(precision)
        recalls.append(recall)
        plt.plot(recall, precision, 'b', alpha=0.15)

    # Calculate mean precision-recall curve
    mean_precision = np.mean(precisions, axis=0)
    mean_recall = np.mean(recalls, axis=0)
    std_precision = np.std(precisions, axis=0)
    precisions_upper = np.minimum(mean_precision + std_precision, 1)
    precisions_lower = mean_precision - std_precision
    pr_auc_mean = np.mean(pr_aucs)
    pr_auc_std = np.std(pr_aucs)
    ci = pr_auc_std / 2

    # Plot mean precision-recall curve
    plt.plot(mean_recall, mean_precision, 'b', label='Mean PR-AUC = {:.2f} $\pm$ {:.2f}'.format(pr_auc_mean, ci))
    plt.fill_between(mean_recall, precisions_lower, precisions_upper, color='grey', alpha=0.3)
    plt.legend(loc='lower right')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.title('Precision-Recall Curve')
    plt.savefig(output_dir + 'PR_curve.png')



if __name__ == "__main__":
    # epi vs nonepi (all marks)
    # stratified_classifier(output_dir=sys.argv[1])
    
    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1']

    for hm in hms:
        stratified_hms_classifier(output_dir=sys.argv[1], hm=hm)