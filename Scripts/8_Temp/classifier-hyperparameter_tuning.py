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


def hms_hyperparamter_tuning(output_dir,hm):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    ## STEP 1: extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info


    ## STEP 2: Handling Multicollinear Features
    def cluster_corr_features(X, threshold):

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
    thresholds = [0, 0.25, 0.5, 1.0, 1.5]
    plt.figure(figsize=(10, 6))

    # set color palette
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))
    alphas = [0.05, 0.25, 0.50, 0.75, 1.0]

    for i, threshold in enumerate(thresholds):
        sf = [val for val in features.columns if val != 'label']
        clustered_features = cluster_corr_features(X=features[sf], threshold=threshold)
        X, y = clustered_features.values, features['label'].values

        # RFECV implementation
        rfecv = RFECV(estimator=clf, step=1, cv=kf, scoring='average_precision', n_jobs=-1)
        rfecv.fit(X, y)

        # Get the grid of cross-validated scores
        cv_scores = rfecv.cv_results_['mean_test_score']

        # Create a feature scores array
        feature_scores = np.array([(i + 1, score) for i, score in enumerate(cv_scores)])

        # Find the index of the maximum score
        max_index = np.argmax(cv_scores)
        max_num_features = feature_scores[max_index, 0]

        # Plot the scores
        plt.plot(feature_scores[:, 0], feature_scores[:, 1], marker='o', alpha = alphas[i],color=color_dict[hm], label=f'Threshold {threshold} ')
        plt.axvline(x=max_num_features, color='gray', linestyle='--', label=f'Optimal # Features : {int(max_num_features)} / {len(clustered_features.columns)}')

    plt.xlabel('Number of Features')
    plt.ylabel('Mean PR AUC Score')
    plt.title(f'Recursive Feature Elimination CV - {hm}')
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_dir + f'Feature_Selection_{hm}.png')
    plt.close()

    print("Plot saved as 'Feature_Selection_{hm}.png'")


if __name__ == "__main__":

    hms = ['H3K27ac', 'H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1', 'H3K27me3']

    for hm in hms:
        print(hm)
        hms_hyperparamter_tuning(output_dir=sys.argv[1], hm=hm)
