import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from scipy.cluster import hierarchy
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier
import json
from sklearn.model_selection import StratifiedKFold
from sklearn.inspection import permutation_importance
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests
import sys
import math
import json
from scipy import stats


def prep(rbps):
    epi = pd.read_csv('0_Files/Post-processing/features_epi.csv', delimiter='\t')
    nonepi = pd.read_csv('0_Files/Post-processing/features_nonepi.csv', delimiter='\t')

    epi["label"] = "epi"
    nonepi["label"] = "nonepi"

    # remove genes with both labels
    common_genes = list(set(epi.gene_name.values.tolist()) & set(nonepi.gene_name.values.tolist()))
    epi = epi[~(epi.gene_name.isin(common_genes))]
    nonepi = nonepi[~(nonepi.gene_name.isin(common_genes))]

    sfs = rbps

    # keep only strong binding events
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


def significane_hms(hm, rbps):

    enriched_epi = []
    enriched_nonepi = []

    epi, nonepi, sfs = prep(rbps)

    epi = epi[epi['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    nonepi = nonepi[nonepi['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]

    prob_dict = {}
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


def cluster_feature_impt(hm, output_dir = sys.argv[1]):

    # set color palette
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))

    ##################################### STEP 1: Load and preprocess the data. #####################################
    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]  # Extract features for hm
    features = features.drop('type', axis=1)  # Drop hm info
    features.fillna(0, inplace=True)  # Fill NaNs with 0
    features = shuffle(features, random_state=42)  # Shuffle the data for randomness

    sf = [val for val in features.columns if val != 'label']
    X = features[sf]

    ##################################### STEP 2: Cluster data.  #####################################
    # Compute the correlation matrix
    corr = spearmanr(X).correlation
    corr = (corr + corr.T) / 2  # Ensure the matrix is symmetric
    np.fill_diagonal(corr, 1)  # Fill the diagonal with 1s
    distance_matrix = 1 - np.abs(corr)  # Conver corr coeff matrix to distance matrix

    # Perform hierarchical clustering
    dist_linkage = hierarchy.ward(squareform(distance_matrix))

    ##################################### STEP 3: Extract clusters. #####################################
    threshold = 1.5
    
    clusters = hierarchy.fcluster(dist_linkage, t=threshold, criterion='distance')

    # Get the cluster labels
    labels = X.columns.to_list()

    # Create a dictionary to hold the clusters
    cluster_dict = {}
    for label, cluster_id in zip(labels, clusters):
        if cluster_id not in cluster_dict:
            cluster_dict[cluster_id] = []
        cluster_dict[cluster_id].append(label)

    ##################################### STEP 4: Visualize clusters. #####################################

    # Plot the dendrogram initially to extract  information
    dendro = hierarchy.dendrogram(
        dist_linkage,
        labels=X.columns.to_list(),
        leaf_rotation=90,
        color_threshold=threshold,
        no_plot=True  # Do not plot the dendrogram yet
    )

    # Get # cluster num
    cluster_colors = dendro['leaves_color_list']
    unique_colors = []
    cluster_counter = 1
    for color in cluster_colors:
        if color not in unique_colors:
            unique_colors.append(color)
            cluster_counter += 1

    # Define a custom color palette for the clusters using color palette
    cmap = plt.colormaps.get_cmap("jet")
    num_clusters = cluster_counter  # Number of unique clusters to color
    colors = [to_hex(cmap(i / num_clusters)) for i in range(num_clusters)]
    hierarchy.set_link_color_palette(colors)    

    # Plot dendrogram
    fig, ax1 = plt.subplots(figsize=(10, 7))
    dendro = hierarchy.dendrogram(
        dist_linkage,
        labels=X.columns.to_list(),
        leaf_rotation=90,
        color_threshold=threshold,
        above_threshold_color=color_dict[hm]
    )

    # Extract color information from the dendrogram
    colormap = {label: color for label, color in zip(dendro['ivl'], dendro['leaves_color_list'])}

    # Extract unique cluster colors in the order they are plotted
    cluster_colors = dendro['leaves_color_list']
    unique_colors = []
    unique_labels = {}
    cluster_counter = 1
    for color in cluster_colors:
        if color not in unique_colors:
            unique_colors.append(color)
            unique_labels[color] = f'Cluster {cluster_counter}'
            cluster_counter += 1


    # Set tick label colors and positions based on the colormap and cluster_dict
    for tick in ax1.get_xticklabels():
        tick_text = tick.get_text()
        tick.set_color(colormap[tick_text])
        # Adjust the position of every alternate cluster group
        for cluster_id, labels in cluster_dict.items():
            if tick_text in labels:
                offset = -0.01 if cluster_id %2 == 0 else -0.06  # Adjust the offset values as needed
                tick.set_y(offset)  # Adjust the position closer or further from the axis
                break

    # title
    ax1.set_title('Clustering RBPs by Binding Sequence Similarity - ' + f'{hm}')
    ax1.set_ylabel('Distance')
    
    plt.tight_layout()
    plt.savefig(output_dir + f'Clustering_{hm}.png')
    plt.close()


    ##################################### STEP 5: Select cluster-representative RBP and get feature importance.  #####################################
    keys = list(cluster_dict.keys())
    keys.sort()
    cluster_dict = {i: cluster_dict[i] for i in keys}
    
    ## get mean of cluster as cluster representative
    for key in cluster_dict.keys():
        features['cluster_' +str(key)] = features[cluster_dict[key]].mean(axis=1)
    selected_features = ['cluster_'+ str(key) for key in list(cluster_dict.keys())]
    cluster_names = dict(zip(list(cluster_dict.keys()), selected_features))  

    # prep data 
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = features[selected_features].values, features['label'].values

    clf =  RandomForestClassifier(n_estimators=100, max_features= "sqrt", bootstrap=True, class_weight='balanced', n_jobs = -1, random_state=0)
    kf = StratifiedKFold(n_splits=10) 

    perm_importances = []
    for i, (train, test) in enumerate(kf.split(X, y)):
            model = clf.fit(X[train], y[train])
            result = permutation_importance(model, X[test], y[test], n_repeats=10, n_jobs=-1, random_state=42, scoring='average_precision')  # Calculate permutation importance
            importance = result.importances_mean
            perm_importances.append(importance)

    ##################################### STEP 6: Visualize feature importance:    #####################################

    # Convert to numpy array
    mean_importance = np.array(perm_importances)

    # Calculate mean and standard deviation of importances
    mean_importance = np.mean(perm_importances, axis=0)
    std_importance = np.std(perm_importances, axis=0)

    # Sort features by mean importance from highest to lowest
    sorted_indices = np.argsort(mean_importance)[::1]
    mean_importance_sorted = mean_importance[sorted_indices]
    std_importance_sorted = std_importance[sorted_indices]
    sf_sorted = np.array(selected_features)[sorted_indices]

    # Set the figure size
    plt.figure(figsize=(20, len(sf_sorted)*0.5))  # Adjust the width and height of the figure as needed

    plt.errorbar(mean_importance_sorted, np.arange(len(sf_sorted)), xerr=std_importance_sorted, fmt='o', color=color_dict[hm], ecolor=color_dict[hm], elinewidth=3, markersize =20, capsize=5)
    plt.barh(sf_sorted, mean_importance_sorted, height=0.5, color='lightgray')

    # get cluster names in order of importance
    value_to_key = {v: k for k, v in cluster_names.items()} ## Create a reverse mapping from values to keys
    cluster_sorted = [value_to_key[value] for value in sf_sorted if value in value_to_key]

    plt.yticks(np.arange(len(sf_sorted)), cluster_sorted, fontsize=40, rotation=45) # Increase space between y-axis ticks

    plt.xlabel(f'Mean Permutation Importance Value', fontsize=15)
    plt.ylabel('Clusters', fontsize=15)
    plt.title(f'Feature Importance - {hm}', fontsize=20)
    plt.tick_params(axis='y', labelsize=8)  # Adjust y-axis tick label size


    plt.tight_layout()
    plt.savefig(output_dir + f'Permutation_Importance_corr_{hm}.png')
    plt.close()

    ##################################### STEP 7: Enrichment Analyses of members of Impt Clusters:    #####################################
    # cluster_sel = {"H3K27ac" : 2,
    #                "H3K27me3" : 2,
    #                "H3K4me3" : 2,
    #                "H3K9me3" : 2,
    #                 "H3K36me3" : 3, 
    #                 "H3K4me1": 2}
    
    # num_clusters = cluster_sel[hm]
    # impt_clusters = cluster_sorted[-num_clusters:]
    # candidate_rbps = []
    # for i in impt_clusters:
    #     candidate_rbps.extend(cluster_dict[i])
    # candidate_rbps = list(set(candidate_rbps))
    # significane_hms(hm, candidate_rbps)


for hm in [ "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]:
     print(hm)
     cluster_feature_impt(hm)