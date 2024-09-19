import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
import shap
import os
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.metrics import average_precision_score



def stratified_hms_classifier( hm):

    output_dir = '0_Files/Post-processing/SHAP'
    # Create a directory specific to the hm value
    hm_dir = os.path.join(output_dir, hm)
    os.makedirs(hm_dir, exist_ok=True)  # Create the directory if it doesn't exist


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
    # sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data, features['label']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)


    # Initialize classifier and cross-validation
    clf = RandomForestClassifier(n_estimators=100, max_features= "sqrt", class_weight='balanced', n_jobs = -1, random_state=0)
    clf.fit(X_train,y_train)

    # Compute SHAP values
    explainer = shap.TreeExplainer(clf)
    shap_values = explainer.shap_values(X_test)

    shap_values_class_1 = np.array(shap_values[:,:,1])


    # # Plot summary plot for SHAP values of the positive class (class 1)
    # shap.summary_plot(shap_values_class_1, X_test, feature_names=sf_data.columns, show=False, max_display=20)

    # # Use matplotlib to add a title
    # plt.title(f'SHAP Summary Plot for {hm}', fontsize=12)

    # # Display the plot
    # plt.savefig(f'{hm_dir}/{hm}_SHAP.png')
    # plt.close()

    mean_abs_shap_values = np.abs(shap_values_class_1).mean(axis=0)
    top_5_feature_indices = np.argsort(mean_abs_shap_values)[-5:]
    top_5_features = sf_data.columns[top_5_feature_indices]

    # Plot SHAP dependence plots for the top 5 features
    i = 1
    for feature in ['CELF1', 'RBFOX3', 'RBFOX2', 'RBFOX1', "IGF2BP2", "MBNL1", 'HNRPLL']:
        print(feature)
        shap.dependence_plot(feature, shap_values_class_1, X_test, show=False, interaction_index='auto')
        plt.title(f'SHAP Dependence Plot for {feature}', fontsize=16)
        # plt.savefig(f'{hm_dir}/{i}_{feature}_{hm}_dependence_plot.png')
        # plt.close()
        plt.show()
        i +=1


def stratified_hms_classifier_cv( hm):

    output_dir = '0_Files/Post-processing/SHAP'
    # Create a directory specific to the hm value
    hm_dir = os.path.join(output_dir, hm)
    os.makedirs(hm_dir, exist_ok=True)  # Create the directory if it doesn't exist

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
    kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=3) 

    # Initialize list to store SHAP values
    shap_values_list = []
    test_set_list = []

    # Loop over folds
    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        
        # Calculate SHAP values
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X[test])
        shap_values_list.append(shap_values)

        # save test set
        test_set_list.append(X[test])

    # Combine SHAP values across all folds
    combined_shap_values = np.concatenate([shap_values[:,:,1] for shap_values in shap_values_list], axis=0)
    combined_X_test = np.concatenate([set for set in test_set_list], axis=0)

    # Plot summary plot for combined SHAP values
    shap.summary_plot(combined_shap_values, combined_X_test, feature_names=sf_data.columns, show=False, max_display=50)

    # Use matplotlib to add a title
    plt.title(f'SHAP Summary Plot for {hm}', fontsize=12)

    # Display the plot
    plt.savefig(f'{hm_dir}/{hm}_SHAP.png')
    plt.close()

    mean_abs_shap_values = np.abs(combined_shap_values).mean(axis=0)
    top_feature_indices = np.argsort(mean_abs_shap_values)[-10:]
    top_features = sf_data.columns[top_feature_indices]

    # write top 20 features into json file
    with open(f"{output_dir}/top_{hm}.txt", 'w') as f:
        for line in list(reversed(list(top_features))):
            f.write(f"{line}\n")

    # # validate_rf(top_features, hm)

    # # Plot SHAP dependence plots of all features
    for feature in features:
        if feature == 'label':
            continue
        shap.dependence_plot(feature, combined_shap_values, combined_X_test, feature_names=sf, show=False, interaction_index=None)
        plt.title(f'{feature} Binding in {hm} Flanks', fontsize=10)
        plt.xlabel(f'Binding Scores of {feature}')
        plt.tight_layout()
        plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.2)

        #Display the plot
        plt.savefig(f'{hm_dir}/{feature}.png')
        plt.close()


def validate_rf(impt_rbps, hm):

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
    sf_data = features[impt_rbps]
    # sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
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
        avg_precision = average_precision_score(y[test], y_score[:, 1])
        pr_aucs.append(avg_precision)

    # Calculate mean PR-AUC and confidence interval
    pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
    pr_auc_ci = float("%.2f" % (pr_auc_range[1] - pr_auc_range[0])) / 2
    mean_pr_auc = "%.2f" % (pr_auc_range[1] - pr_auc_ci)
    print(f'EPI vs NONEPI - {hm}:', mean_pr_auc)


if __name__ == "__main__":
    
    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1']

    for hm in hms:
        print(hm)
        stratified_hms_classifier_cv(hm=hm)