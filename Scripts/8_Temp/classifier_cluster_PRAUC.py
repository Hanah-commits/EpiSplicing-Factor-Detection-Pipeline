import pandas as pd
import numpy as np
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from collections import defaultdict
from sklearn.metrics import precision_recall_curve, auc, average_precision_score


def hms_hyperparamter_tuning(output_dir,hm):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    ## STEP 1: extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info


    ## STEP 2: Handling Multicollinear Features
    def cluster_corr_features(X, threshold):
        
        X = X.applymap(lambda val: 0 if val < 2 else val)

        if threshold == 0:
            return X
        else:

            corr = spearmanr(X).correlation
            corr = (corr + corr.T) / 2
            np.fill_diagonal(corr, 1)
            distance_matrix = 1 - np.abs(corr)
            dist_linkage = hierarchy.ward(squareform(distance_matrix))
            cluster_ids = hierarchy.fcluster(dist_linkage, threshold, criterion="distance")
            cluster_id_to_feature_ids = defaultdict(list)
            for idx, cluster_id in enumerate(cluster_ids):
                cluster_id_to_feature_ids[cluster_id].append(idx)
            
            ## Extract clusters
            clusters = hierarchy.fcluster(dist_linkage, t=threshold, criterion='distance')
            labels = X.columns.to_list() # Get the cluster labels

            # Create a dictionary to hold the clusters
            cluster_dict = {}
            for label, cluster_id in zip(labels, clusters):
                if cluster_id not in cluster_dict:
                    cluster_dict[cluster_id] = []
                cluster_dict[cluster_id].append(label)

            ## Get mean of cluster as cluster representative
            keys = list(cluster_dict.keys())
            keys.sort()
            cluster_dict = {i: cluster_dict[i] for i in keys}
            for key in cluster_dict.keys():
                X['cluster_' +str(key)] = X[cluster_dict[key]].mean(axis=1)
            clustered_features = ['cluster_'+ str(key) for key in list(cluster_dict.keys())]

            return X[clustered_features]


    ## STEP 3: Build Model, intialize CV
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    clf =  RandomForestClassifier(n_estimators=100, max_features= "sqrt", bootstrap=True, class_weight='balanced', n_jobs = -1, random_state=0)
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42) 

    ##  STEP 4: RFECV implementation
    # Thresholds to iterate over
    thresholds = [0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5]

    # set color palette
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))

    # initialize lists to store prauc of every iteration
    num_clusters = []
    mean_pr_aucs = []

    for i, threshold in enumerate(thresholds):
        sf = [val for val in features.columns if val != 'label']
        clustered_features = cluster_corr_features(X=features[sf], threshold=threshold)
        X, y = clustered_features.values, features['label'].values

        # Initialize arrays to store PR-AUC values
        pr_aucs = []

        # Loop over folds
        for i, (train, test) in enumerate(kf.split(X, y)):
            model = clf.fit(X[train], y[train])
            y_score = model.predict_proba(X[test])
            avg_precision = average_precision_score(y[test], y_score[:, 1])
            pr_aucs.append(avg_precision)
        
        # Calculate mean PR-AUC and confidence interval
        pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
        pr_auc_ci = (pr_auc_range[1] - pr_auc_range[0]) / 2
        mean_pr_auc = pr_auc_range[1] - pr_auc_ci
        mean_pr_auc = float("%.2f" % mean_pr_auc)

        mean_pr_aucs.append(mean_pr_auc)
        num_clusters.append(len(clustered_features.columns))

    # Plotting after the loop
    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(thresholds)), mean_pr_aucs, color=color_dict[hm])
    plt.xlabel('Clustering Threshold')
    plt.ylabel('Mean PR AUC')
    plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.xticks(range(len(thresholds)), thresholds)
    plt.title(f'Effect of RBP Clustering on Classifier Performance - {hm}')
    
    # Add text on each bar
    for i, bar in enumerate(bars):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2, f'{num_clusters[i]} clusters', 
                ha='center', va='center', color='white', fontsize=8)

    # Save the plot
    plt.savefig(output_dir + f'Cluster_PRAUC_{hm}.png')
    plt.close()

    print(f"Plot saved as 'Cluster_PRAUC_{hm}.png'")


if __name__ == "__main__":

    hms = ['H3K27ac', 'H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1', 'H3K27me3']

    for hm in hms:
        print(hm)
        hms_hyperparamter_tuning(output_dir=sys.argv[1], hm=hm)
