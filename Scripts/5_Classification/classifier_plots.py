from pathlib import Path
import pandas as pd
import numpy as np
import os, csv
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate
import shap
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def SHAP_imptRBPs_plot(hm, output_dir):

    op_dir =  f"{output_dir}/SHAP/impt_RBPs"
    Path(op_dir).mkdir(parents=True, exist_ok=True)

    #get SHAP epi and nonepiRBPs to plot
    features_to_plot = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPs/SHAP_epiRBPs/rbps_{hm}.txt", "r")
    features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPs/SHAP_nonepiRBPs/rbps_{hm}.txt", "r")
    features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])

    # get RBPs corr with SHAP RBPs to plot
    corr_features_to_plot = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPs/epiRBPs_{hm}.txt", "r")
    corr_features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPs/nonepiRBPs_{hm}.txt", "r")
    corr_features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    features_to_plot = list(set(features_to_plot + corr_features_to_plot))

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
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
    clf = RandomForestClassifier(class_weight='balanced', n_jobs = -1, random_state=0, n_estimators=200)
    kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10) 

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

    # Plot summary plot for combined SHAP values (only epi and nonepirbpS)

    # Get indices of the features to plot
    feature_indices = [sf_data.columns.get_loc(f) for f in features_to_plot if f in sf_data.columns]
    # Filter combined SHAP values and test set based on the selected features
    filtered_shap_values = combined_shap_values[:, feature_indices]
    filtered_X_test = combined_X_test[:, feature_indices]
    filtered_feature_names = [sf_data.columns[i] for i in feature_indices]

    # Plot summary plot for filtered features
    shap.summary_plot(filtered_shap_values, filtered_X_test, feature_names=filtered_feature_names, show=False, max_display=50)

    # Use matplotlib to add a title
    # plt.title(f'SHAP Summary Plot for {hm}', fontsize=12)

    # # Display the plot
    plt.savefig(f'{op_dir}/{hm}_SHAP_impt.png')
    plt.close()


def correlated_to_shap_RBPs(hm):

    print(hm)

    for mode in ['epi', 'nonepi']:

        # Get episplicing RBPs
        rbps_file = open(f"0_Files/Post-processing/{mode}RBPs/SHAP_{mode}RBPs/rbps_{hm}.txt", "r")
        shap_rbps = list(set([rbp for rbp in rbps_file.read().split('\n') if rbp]))

        features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
        features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
        features = pd.concat([features1,features2], axis=1)
        features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
        features.fillna(0,inplace=True)
        # features = shuffle(features, random_state=42)

        # extract features for hm
        features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
        # features = features[features.label =='epigene']
        features = features.drop(['type','label'], axis=1) # drop hm info
        features = features.applymap(lambda val: 0 if val < 2 else val)

        # Perform Pearson correlation
        rbps = [val for val in features.columns]
        corr_coeffs_pvalues = []
        for i in range(len(rbps)):
            for j in range(i + 1, len(rbps)):
                col1 = rbps[i]
                col2 = rbps[j]
                corr, p_value = pearsonr(features[col1], features[col2])
                corr_coeffs_pvalues.append((col1, col2, corr, p_value))

        corr_df = pd.DataFrame(corr_coeffs_pvalues, columns=["RBP1", "RBP2", "coeff", "pval"]) # convert to df
        corr_df["adj_pval"] = multipletests(corr_df["pval"], method="fdr_bh")[1] # adjust pvalues


        # filtered_corr = corr_df[(corr_df["coeff"] >= 0.5) & (corr_df["adj_pval"] <= 0.05)] #get correlated rbps
        filtered_corr = corr_df[(corr_df["coeff"] >= 0.7) & (corr_df["adj_pval"] <= 0.05)] #get highly correlated rbps

        # get rbps hghly correlated with episplicicng RBPs
        shap_rbp_corr = filtered_corr[filtered_corr.RBP1.isin(shap_rbps)]
        final_rbps = list(set(shap_rbp_corr.RBP2.values.tolist() + shap_rbps))
        final_rbps.sort()
        with open(f"0_Files/Post-processing/{mode}RBPs/{mode}RBPs_{hm}.txt", 'w') as f:
            for line in final_rbps:
                f.write("%s\n" % line)

        print(f'\n\nManually-selected {mode}RBPs:', len(shap_rbps), f'\nManually-selected + Correlated {mode}RBPs:' ,len(final_rbps), '\n\n')


def imptRBPs_overlap():

    for mode in ['epi', 'nonepi']:
        hm_rbps = {}
        hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
        for hm in hms:
            rbps_file = open(f"0_Files/Post-processing/{mode}RBPs/{mode}RBPs_{hm}.txt", "r")
            hm_rbps[hm] = [rbp for rbp in rbps_file.read().split('\n') if rbp]

        
        all_rbps = list(set([rbp for elem_list in hm_rbps.values() for rbp in elem_list]))
        all_rbps.sort()

        # Get occurence of epiRBP across all HMs

        # Open a TSV file to write the results
        with open(f"0_Files/Post-processing/{mode}RBPs/{mode}RBPs.tsv", 'w', newline='') as tsvfile:
            tsv_writer = csv.writer(tsvfile, delimiter='\t')

            # Write the header row
            tsv_writer.writerow([
                'RBP', 'H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3'
            ])

            tsv_rbps = {rbp:[] for rbp in all_rbps}
            for rbp in all_rbps:
                for hm in hms:
                    if rbp in hm_rbps[hm]:
                        tsv_rbps[rbp].append(1)
                    else:
                        tsv_rbps[rbp].append(0)

                tsv_writer.writerow([rbp] + tsv_rbps[rbp]) # write occurrence of each RBP


if __name__ == "__main__":

    op_dir = str(Path(os.getcwd())) + "/0_Files/Post-processing/Analyses"
    Path(op_dir).mkdir(parents=True, exist_ok=True)

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:
        ## STEP 0: Manually select imptRBPs (Epi and nonepi)
        correlated_to_shap_RBPs(hm=hm) # STEP 1: RBPs strongly correlated with manually-slected important RBPs
        SHAP_imptRBPs_plot(hm, output_dir=op_dir) # STEP 2: SHAP plot of Manually-selected + Correlated epi/nonepiRBPs

    imptRBPs_overlap() # Table3,4
        
    