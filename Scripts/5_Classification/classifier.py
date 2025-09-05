from pathlib import Path
import pandas as pd
import numpy as np
import sys, os, csv
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate
from sklearn.metrics import precision_recall_curve, average_precision_score, make_scorer, fbeta_score, confusion_matrix
from sklearn.utils import shuffle
import shap
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from sklearn.dummy import DummyClassifier
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def all_prauc_together():

    os.makedirs("0_Files/Post-processing/Analyses/PRAUC/", exist_ok=True)  # Create the directory if it doesn't exist
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    types = ['epi', 'nonepi', 'epi_nonspliced']
    label = ['epigene', 'non-epigene', 'epi_nonspliced_gene']

    for i in range(len(types)):
        for j in range(len(types)):
            if i >= j:
                continue

            # Initialize plot
            fig, axes = plt.subplots(1, 2, figsize=(14,6), sharey=True)

            for hm in hms:

                features1 = pd.read_csv(f'0_Files/Post-processing/features_all_{types[i]}_vs_{types[j]}_132.csv', delimiter='\t')
                features2 = pd.read_csv(f'0_Files/Post-processing/features_all_{types[i]}_vs_{types[j]}_47.csv', delimiter='\t')
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
                sf_data = features[sf]        
                sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val) # strong binding events
                
                features['label'] = features['label'].map({label[i]: 1, label[j]: 0}).astype(int)
                X, y = sf_data.values, features['label'].values

                # Classifiers
                clf = RandomForestClassifier(class_weight='balanced', n_jobs = -1, random_state=0, n_estimators=200)
                dummy_clf = DummyClassifier(strategy='stratified', random_state=42)

                # cross-validation
                kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10)


                def evaluate_model(clf_model, X, y, SHAP):
                    all_precisions = []
                    mean_pr_aucs = []
                    # loop over folds
                    for i, (train, test) in enumerate(kf.split(X, y)):
                        model = clf_model.fit(X[train], y[train])
                        y_score = model.predict_proba(X[test])
                        avg_precision = average_precision_score(y[test], y_score[:, 1])
                        mean_pr_aucs.append(avg_precision)      
                        precision, recall, _ = precision_recall_curve(y[test], y_score[:, 1])
                        all_precisions.append(np.interp(np.linspace(0, 1, 100), recall[::-1], precision[::-1])) ## Interpolate precision values for consistent recall points
                    # calculate the mean precision-recall curve
                    mean_precision = np.mean(all_precisions, axis=0)
                    mean_recall = np.linspace(0, 1, 100)
                    mean_pr_auc = np.mean(mean_pr_aucs)
                    pr_auc_ci = np.std(mean_pr_aucs) * 1.96  # 95% CI approximation

                    return mean_recall, mean_precision, mean_pr_auc, pr_auc_ci
                
                # classifiers
                rf_recall, rf_precision, rf_auc, rf_ci = evaluate_model(clf, X, y, SHAP=True)
                dummy_recall, dummy_precision, dummy_auc, dummy_ci = evaluate_model(dummy_clf, X, y, SHAP=False)

                # RF model
                axes[0].plot(rf_recall, rf_precision, color=color_dict[hm], label = f'{hm} Mean PR-AUC = {rf_auc:.2f} ± {rf_ci:.2f}')
                axes[0].set_xlabel("Recall", fontsize=14)
                axes[0].set_ylabel("Precision", fontsize=14)
                axes[0].legend(loc='lower left', fontsize=12)

                # Baseline
                axes[1].plot(dummy_recall, dummy_precision, color=color_dict[hm], label = f'Baseline {hm} Mean PR-AUC = {dummy_auc:.2f} ± {dummy_ci:.2f}', linestyle="--")
                axes[1].set_xlabel("Recall", fontsize=14)
                axes[1].legend(loc='upper right', fontsize=12)

            plt.tight_layout(pad = 0, h_pad=0, w_pad=0,rect=[0, 0, 1, 1])
            plt.savefig(f"0_Files/Post-processing/Analyses/PRAUC/PRAUC_{types[i]}_vs_{types[j]}.png", bbox_inches='tight', dpi=300)
            plt.close()


def plot_metrics():

    metrics_df = pd.read_csv(f"0_Files/Post-processing/Analyses/Metrics/metrics.tsv", delimiter='\t')
    metrics = ['Mean Recall', 'Mean Specificity', 'Mean Precision']

    values = []
    errors = []

    for metric in metrics:
        metric_values = metrics_df[metric].str.extract(r'([\d.]+) ± ([\d.]+)').astype(float)
        values.append(metric_values[0].values)
        errors.append(metric_values[1].values)

    # Create the bar plot
    x = np.arange(len(metrics_df))  # Group positions
    bar_width = 0.2
    group_spacing = 0.85
    patterns = ['**', 'xx', '..']
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    colors = ["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]



    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot bars for each metric
    for i, metric in enumerate(metrics):
        for j, group_color in enumerate(colors):
            ax.bar(
                x[j] + i * bar_width, values[i][j], bar_width, yerr=errors[i][j], capsize=5,
                color=group_color, edgecolor='white', hatch=patterns[i], ecolor='grey'
            )


    # Add labels and legend
    ax.set_xticks(x + bar_width * (len(metrics) - 1) / 2)
    ax.set_xticklabels(metrics_df["Model"])

    # Set x-tick labels iteratively for each group
    ax.set_ylabel("Mean Scores", fontsize=12)
    # ax.set_title("Model Metrics", fontsize=14)

    # Create custom legend
  
    pattern_patches = [mpatches.Patch(facecolor='grey', edgecolor='white', hatch=pattern, label=metric) for pattern, metric in zip(patterns, metrics)]
    ax.legend(handles=pattern_patches, title="Metrics", loc='upper left', bbox_to_anchor=(0.715, 1.0), prop={'size': 8}, handleheight=2, handlelength=3)


    plt.tight_layout()
    plt.savefig(f"0_Files/Post-processing/Analyses/Metrics/metrics.png", bbox_inches='tight', dpi=300)
    plt.close()


def metrics():

    os.makedirs("0_Files/Post-processing/Analyses/Metrics/", exist_ok=True)  # Create the directory if it doesn't exist

    # Custom scoring function for specificity
    def specificity_score(y_true, y_pred):
        cm = confusion_matrix(y_true, y_pred)
        tn, fp, fn, tp = cm.ravel()
        specificity = tn / (tn + fp)
        return specificity

    # Create a specificity scorer, f-beta scores
    specificity_scorer = make_scorer(specificity_score)
    fbeta_scorer = make_scorer(fbeta_score, beta=0.5)

    # Open a TSV file to write the results
    with open(f"0_Files/Post-processing/Analyses/Metrics/metrics.tsv", 'w', newline='') as tsvfile:
        tsv_writer = csv.writer(tsvfile, delimiter='\t')

        # Write the header row
        tsv_writer.writerow([
            'Model', 'Mean F1', 'Mean F-beta', 'Mean Recall', 
            'Mean Specificity', 'Mean Precision', 'Mean PRAUC'
        ])


        hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
        for hm in hms:
        
            features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
            features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
            features = pd.concat([features1,features2], axis=1)
            features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
            features.fillna(0,inplace=True)
            features = shuffle(features, random_state=42)

            # extract features for hm
            features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
            features = features.drop('type', axis=1) # drop hm info


            sf = [val for val in features.columns if val != 'label']

            ## STEP 0: keep only strong binding events
            sf_data = features[sf]
            sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

            features['label'] = features['label'].map({'epigene': 0, 'non-epigene': 1}).astype(int)
            X, y = sf_data.values, features['label'].values

            ## STEP 1: Initialize classifier and cross-validation
            clf = RandomForestClassifier(class_weight='balanced', n_jobs = -1, random_state=0, n_estimators=200)
            kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10) 

            # Evaluate  model
            scoring = {
                'f1': 'f1',
                'f-beta': fbeta_scorer,
                'recall': 'recall',
                'precision': 'precision',
                'average_precision': 'average_precision',
                'specificity': specificity_scorer
            }
            scores = cross_validate(clf, X, y, scoring=scoring, cv=kf)
            print(f'\n\n{hm}')

            # Calculate mean and confidence intervals
            def mean_ci(values):
                mean = np.mean(values)
                std_error = np.std(values, ddof=1) / np.sqrt(len(values))
                ci = 1.96 * std_error #95% CI
                return f'{mean:.3f} ± {ci:.3f}'
                
            mean_f1 = mean_ci(scores['test_f1'])
            mean_fbeta = mean_ci(scores['test_f-beta'])
            mean_recall = mean_ci(scores['test_recall'])
            mean_specificity = mean_ci(scores['test_specificity'])
            mean_precision = mean_ci(scores['test_precision'])
            mean_avg_precision = mean_ci(scores['test_average_precision'])
            
            # Write the model name and scores with CIs to the TSV file
            tsv_writer.writerow([
                hm, mean_f1, mean_fbeta, mean_recall,
                mean_specificity, mean_precision, mean_avg_precision
            ])

    plot_metrics()


def confusion_matrix_plot(hm):

    os.makedirs("0_Files/Post-processing/Analyses/Metrics/", exist_ok=True)  # Create the directory if it doesn't exist

    features1 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_epi_vs_nonepi_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)
    features = shuffle(features, random_state=42)

    # extract features for hm
    features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
    features = features.drop('type', axis=1) # drop hm info


    sf = [val for val in features.columns if val != 'label']

    ## STEP 0: keep only strong binding events
    sf_data = features[sf]
    sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    ## STEP 1: Initialize classifier and cross-validation
    clf = RandomForestClassifier(class_weight='balanced', n_jobs = -1, random_state=0, n_estimators=200)
    kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10) 

            
    y_true = []
    y_pred = []
    y_proba = []

    # Manually perform cross-validation to collect predictions
    for train_idx, test_idx in kf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        clf.fit(X_train, y_train)
        
        # Store the true labels and predictions for confusion matrix
        y_true.extend(y_test)
        y_pred.extend(clf.predict(X_test))

    # Convert lists to numpy arrays
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)


    ## STEP 3: Construct Confusion Matrix

    cm = confusion_matrix(y_true, y_pred)
    # Convert the confusion matrix to percentages
    cm_percentage = cm.astype('float') / cm.sum() * 100

    plt.figure(figsize=(8, 6))
    sns.heatmap(cm_percentage, annot=True, fmt=".2f", cmap="binary", cbar_kws={'label': 'Percentage'}, vmin=1, vmax=100)

    # Add labels
    plt.xlabel("Predicted Labels")
    plt.ylabel("True Labels")
    # plt.title(f'Confusion Matrix - {hm}')
    plt.savefig(f"0_Files/Post-processing/Analyses/Metrics/{hm}_CM.png", bbox_inches='tight', dpi=300)
    plt.close()


def SHAP( hm):

    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    hm_dir = f"0_Files/Post-processing/Analyses/SHAP/{hm}"
    os.makedirs(hm_dir, exist_ok=True)  # Create the directory if it doesn't exist

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

    # Plot summary plot for combined SHAP values
    shap.summary_plot(combined_shap_values, combined_X_test, feature_names=sf_data.columns, show=False, max_display=50)

    # Use matplotlib to add a title
    plt.title(f'SHAP Summary Plot for {hm}', fontsize=12)

    # # Display the plot
    plt.savefig(f'{hm_dir}/{hm}_SHAP.png')
    plt.close()

    mean_abs_shap_values = np.abs(combined_shap_values).mean(axis=0)


    # # Plot SHAP dependence plots of all features
    # for feature in sf:
    #     feature_index = sf.index(feature)
    #     shap.dependence_plot(feature_index, combined_shap_values, combined_X_test, feature_names=sf, show=False, interaction_index=None)

    #     # to use custom color 
    #     ax = plt.gca()
    #     for artist in ax.collections: # remove default scatter plots
    #         artist.remove()
    #     plt.scatter(combined_X_test[:,feature_index], combined_shap_values[:, feature_index], s=3, color=color_dict[hm])

    #     # plt.title(f'{feature} Binding in {hm}-marked Exon Flanks', fontsize=12)
    #     plt.xlabel(f'Predicted Binding Scores of {feature}', fontsize=12)
    #     plt.ylabel(f'SHAP Values for {feature}', fontsize=12)
    #     plt.xticks(fontsize=12)
    #     plt.yticks(fontsize=12)
    #     plt.tight_layout()
    #     plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.2)

    #     #Display the plot
    #     plt.savefig(f'{hm_dir}/{feature}_{hm}.png')
    #     plt.close()

    # Create a DataFrame for feature importance
    feature_importance_df = pd.DataFrame({
        'Feature': sf_data.columns, 
        'Mean Absolute SHAP Value': mean_abs_shap_values
    })

    # Sort the features by importance
    feature_importance_df = feature_importance_df.sort_values(by='Mean Absolute SHAP Value', ascending=False)
    feature_importance_df.to_csv(f'0_Files/Post-processing/Analyses/SHAP/{hm}_shap_feature_importance.csv', sep='\t', index=False)


if __name__ == "__main__":

    all_prauc_together()
    metrics()

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:
        confusion_matrix_plot(hm)
        SHAP(hm=hm)

        # empty files for manual selection of impt RBPs based on SHAP plots
        for mode in ['epi', 'nonepi']:
            output_dir = str(Path(os.getcwd())) + f"/0_Files/Post-processing/{mode}RBPs/SHAP_{mode}RBPs"
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            open(f'{output_dir}/rbps_{hm}.txt', 'a').close()
    