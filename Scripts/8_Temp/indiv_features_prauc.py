import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, auc
from sklearn.utils import shuffle
import json
import seaborn as sns
from matplotlib.ticker import MaxNLocator



def stratified_hms_classifier_indiv(hm):

    features = pd.read_csv('0_Files/Post-processing/features_all.csv', delimiter='\t')
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info

    if len(features) == 0:
        return
    
     ## EPI-ENRICHED RBPS

    ## method 1: use only enriched rbps
    ## method 2: Use all non-enriched rbps

    ## method 3: use all rbps: one at a time
    features['pseudo_feature'] = 0 # add fake feature
    sfs = [val for val in features.columns if val != 'label' and val != 'pseudo_feature']

    # keep only strong binding events
    features[sfs] = features[sfs].applymap(lambda val: 0 if val < 2 else val)

    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    auc_sf = {}
    for sf in sfs:
        sf_data = features[[sf,'pseudo_feature']]
        X, y = sf_data.values, features['label'].values

        # Initialize classifier and cross-validation
        clf = RandomForestClassifier(n_estimators=100, max_features= "sqrt", class_weight='balanced', n_jobs = -1, random_state=0)
        kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

        # Initialize arrays to store PR-AUC values
        pr_aucs = []

        # Loop over folds
        for i, (train, test) in enumerate(kf.split(X, y)):
            model = clf.fit(X[train], y[train])
            y_score = model.predict_proba(X[test])
            precision, recall, _ = precision_recall_curve(y[test], y_score[:, 1])
            pr_auc = auc(recall, precision)
            pr_aucs.append(pr_auc)

        # Calculate mean PR-AUC and confidence interval
        pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
        pr_auc_ci = float("%.2f" % (pr_auc_range[1] - pr_auc_range[0])) / 2
        mean_pr_auc = "%.2f" % (pr_auc_range[1] - pr_auc_ci)
        auc_sf[sf] = mean_pr_auc

    with open(f"0_Files/Post-processing/{hm}_aucs.json", "w") as outfile:
        json.dump(auc_sf, outfile, indent=4, sort_keys=False)

        
def visualise_dist(hm):

    # colordict
    hms = [  "H3K27ac","H3K27me3","H3K4me3","H3K9me3", "H3K36me3", "H3K4me1"]
    color_dict = dict(zip(hms,["#AD50D3", "#FA5557", "#FA55BA", "#FCB10C", "#91C820", "#33ABCC"]))

    # read aucs
    with open(f"0_Files/Post-processing/{hm}_aucs.json") as f:
        aucs = json.load(f)

    #plot
    dist = sns.displot(sorted(list(aucs.values())),kde=True, color=color_dict[hm], height=10, aspect=1.0)

    # Adjust y-axis ticks to be natural numbers and add grid
    for ax in dist.axes.flat:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        # Adjust tick label size
        ax.tick_params(axis='both', which='major', labelsize=8.5)

    # Set title and axis labels with increased font size
    dist.fig.suptitle(f'Distribution of AUCs - {hm}', fontsize=15)
    dist.set_axis_labels('AUC', 'Frequency', fontsize=12)


    # Adjust subplot parameters to ensure title is not cut off
    plt.subplots_adjust(top=0.9, bottom=0.07)
    plt.tight_layout()
    plt.savefig(f'0_Files/Post-processing/{hm}_aucs.png')


if __name__ == "__main__":
    
    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1']

for hm in hms:
        stratified_hms_classifier_indiv(hm=hm)
        visualise_dist(hm)
