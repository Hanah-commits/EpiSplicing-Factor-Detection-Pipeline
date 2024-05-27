import pandas as pd
import numpy as np
import sys
import json
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.model_selection import StratifiedKFold
from sklearn.inspection import permutation_importance


def stratified_classifier(output_dir, rbp_classes):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features = features.drop('type', axis=1) # drop hm info
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    sf = [val for val in features.columns if val != 'label']
    
    sf_data = features[sf]
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    clf = RandomForestClassifier(n_estimators=100, max_features= "sqrt", class_weight='balanced', bootstrap=True, oob_score=True, n_jobs = -1, random_state=0) 
    kf = StratifiedKFold(n_splits=10)

    impt_scores = []

    for i, (train, test) in enumerate(kf.split(X, y)):
        
        model = clf.fit(X[train], y[train])
        result = permutation_importance(model, X[test], y[test], n_repeats=10, n_jobs=-1, random_state=42)  # Calculate permutation importance
        importance = result.importances_mean
        impt_scores.append(importance)

    # mean permutation importance across all folds
    mean_importance = np.mean(impt_scores, axis=0)

    # Sort features by mean importance from highest to lowest
    sorted_indices = np.argsort(mean_importance)[::1]
    mean_importance_sorted = mean_importance[sorted_indices]
    sf_sorted = np.array(sf)[sorted_indices]

    with open('0_Files/impt_features/scores.txt', 'w') as f:
        for item in mean_importance_sorted:
            f.write("%s\n" % item)

    with open('0_Files/impt_features/sfs.txt', 'w') as f:
        for item in sf_sorted:
            f.write("%s\n" % item)


    ### Plot permutation importance

    # Define a color map
    color_map = {'RNA metabolism protein(PC00031)': '#FFA500', 'RNA splicing factor(PC00148)': '#C20078', 'RNA processing factor(PC00147)': '#029386', 'Misc': '#00FFFF', 'Unknown': '#FFC0CB'}
    
    # Set the figure size
    plt.figure(figsize=(50, len(sf_sorted)*0.2))  # Adjust the width and height of the figure as needed

    plt.yticks(np.arange(len(sf_sorted)), sf_sorted, fontsize=15, rotation=45) # Increase space between y-axis ticks
    bars = plt.barh(sf_sorted, mean_importance_sorted, height=0.5)


    # mark SF RBPs in red
    for label in sf:
        for key, values in rbp_classes.items():
            if label in values:
                index = np.where(sf_sorted == label)[0][0]
                color = color_map.get(key, 'black')  # Get color from color_map based on key, default to black if key not found
                bars[index].set_color(color)  # Set color based on the color map


    plt.xlabel(f'Mean Permutation Importance (K=10 Folds)', fontsize=20)
    plt.ylabel('RBPs', fontsize=25)
    plt.title('Permutation Importance', fontsize=25)
    plt.tick_params(axis='y', labelsize=8)  # Adjust y-axis tick label size

    # Create legend based on color_map
    legend_elements = [Patch(facecolor=color, edgecolor='black', label=key) for key, color in color_map.items()]
    plt.legend(handles=legend_elements, loc='lower right', fontsize='30')

    plt.tight_layout()
    plt.savefig(output_dir + '0_Files/impt_features/Permutation_Importance.png')
    plt.close()


def stratified_hms_classifier(output_dir, hm, rbp_classes):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info

    sf = [val for val in features.columns if val != 'label']
    sf_rbps = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTB3', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']


    sf_data = features[sf]
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    clf = RandomForestClassifier(n_estimators=100, max_features= "sqrt", class_weight='balanced', bootstrap=True, oob_score=True, n_jobs = -1, random_state=0) 
    kf = StratifiedKFold(n_splits=10)

    impt_scores = []

    for i, (train, test) in enumerate(kf.split(X, y)):
        
        model = clf.fit(X[train], y[train])
        result = permutation_importance(model, X[test], y[test], n_repeats=10, n_jobs=-1, random_state=42)  # Calculate permutation importance
        importance = result.importances_mean
        impt_scores.append(importance)

    # mean permutation importance across all folds
    mean_importance = np.mean(impt_scores, axis=0)

    # Sort features by mean importance from highest to lowest
    sorted_indices = np.argsort(mean_importance)[::1]
    mean_importance_sorted = mean_importance[sorted_indices]
    sf_sorted = np.array(sf)[sorted_indices]

    with open(f'0_Files/impt_features/scores_{hm}.txt', 'w') as f:
        for item in mean_importance_sorted:
            f.write("%s\n" % item)

    with open(f'0_Files/impt_features/sfs_{hm}.txt', 'w') as f:
        for item in sf_sorted:
            f.write("%s\n" % item)

    # Plot permutation importance
    # Define a color map
    color_map = {'RNA metabolism protein(PC00031)': '#FFA500', 'RNA splicing factor(PC00148)': '#C20078', 'RNA processing factor(PC00147)': '#029386', 'Misc': '#00FFFF', 'Unknown': '#FFC0CB'}
    
    # Set the figure size
    plt.figure(figsize=(50, len(sf_sorted)*0.2))  # Adjust the width and height of the figure as needed

    plt.yticks(np.arange(len(sf_sorted)), sf_sorted, fontsize=15, rotation=45) # Increase space between y-axis ticks
    bars = plt.barh(sf_sorted, mean_importance_sorted, height=0.5)


    # mark SF RBPs in red
    for label in sf:
        for key, values in rbp_classes.items():
            if label in values:
                index = np.where(sf_sorted == label)[0][0]
                color = color_map.get(key, 'black')  # Get color from color_map based on key, default to black if key not found
                bars[index].set_color(color)  # Set color based on the color map


    plt.xlabel(f'Mean Permutation Importance (K=10 Folds)', fontsize=25)
    plt.ylabel('RBPs', fontsize=25)
    plt.title(f'Permutation Importance - {hm}', fontsize=25)
    plt.tick_params(axis='y', labelsize=8)  # Adjust y-axis tick label size

    # Create legend based on color_map
    legend_elements = [Patch(facecolor=color, edgecolor='black', label=key) for key, color in color_map.items()]
    plt.legend(handles=legend_elements, loc='lower right', fontsize=30)

    plt.tight_layout()
    plt.savefig(output_dir + f'0_Files/impt_features/Permutation_Importance_{hm}.png')
    plt.close()


if __name__ == "__main__":

    with open('HelperFunctions/RBP_Classes.json') as f:
            RBP_classes = json.load(f)

    stratified_classifier(output_dir=sys.argv[1], rbp_classes=RBP_classes)

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1']

    for hm in hms:
        stratified_hms_classifier(output_dir=sys.argv[1], hm=hm, rbp_classes=RBP_classes)
