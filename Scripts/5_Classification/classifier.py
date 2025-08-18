import pandas as pd
import numpy as np
import sys, os, csv
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate
from sklearn.metrics import precision_recall_curve, average_precision_score, make_scorer, fbeta_score, confusion_matrix
from sklearn.utils import shuffle
import shap
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


def all_prauc_together():

    os.makedirs("0_Files/Post-processing/Analyses/PRAUC/", exist_ok=True)  # Create the directory if it doesn't exist
    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    # Initialize plot
    fig, axes = plt.subplots(1, 2, figsize=(14,6), sharey=True)

    for hm in hms:

        features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
        features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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
        
        features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
        X, y = sf_data.values, features['label'].values

        # Classifiers
        clf = RandomForestClassifier(class_weight='balanced', n_jobs = -1, random_state=0, n_estimators=200)
        dummy_clf = DummyClassifier(strategy='uniform', random_state=42)

        # cross-validation
        kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10)         

        def evaluate_model(clf_model, X, y):

            all_precisions = []
            mean_pr_aucs = []

            # Loop over folds
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
        rf_recall, rf_precision, rf_auc, rf_ci = evaluate_model(clf, X, y)
        dummy_recall, dummy_precision, dummy_auc, dummy_ci = evaluate_model(dummy_clf, X, y)

        # RF model
        axes[0].plot(rf_recall, rf_precision, color=color_dict[hm], label = f'{hm} Mean PR-AUC = {rf_auc:.2f} ± {rf_ci:.2f}')
        axes[0].set_xlabel("Recall")
        axes[0].set_ylabel("Precision")
        axes[0].legend(loc='lower left')

        # Baseline
        axes[1].plot(dummy_recall, dummy_precision, color=color_dict[hm], label = f'Baseline {hm} Mean PR-AUC = {dummy_auc:.2f} ± {dummy_ci:.2f}', linestyle="--")
        axes[1].set_xlabel("Recall")
        axes[1].legend(loc='lower left')

    plt.tight_layout(pad = 0, h_pad=0, w_pad=0,rect=[0, 0, 1, 1])
    plt.savefig(f"0_Files/Post-processing/Analyses/PRAUC/PRAUC_affinities.png", bbox_inches='tight', dpi=300)
    plt.close()
    

def stratified_hms_classifier(output_dir, hm):

    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    parameters = {
        "H3K27ac" : {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200},
        "H3K27me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K36me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K9me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K4me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 5, 'n_estimators': 200}
    }

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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
    clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
    kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=3) 

    # Initialize arrays to store PR-AUC values
    pr_aucs = []

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')

    # Loop over folds
    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        avg_precision = average_precision_score(y[test], y_score[:, 1])
        pr_aucs.append(avg_precision)
        precision, recall, _ = precision_recall_curve(y[test], y_score[:, 1])
        plt.plot(recall, precision, color_dict[hm])

    # Calculate mean PR-AUC and confidence interval
    pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
    pr_auc_ci = float("%.2f" % (pr_auc_range[1] - pr_auc_range[0])) / 2
    mean_pr_auc = "%.2f" % (pr_auc_range[1] - pr_auc_ci)
    print(f'EPI vs NONEPI - {hm}:', mean_pr_auc)

    # Plot single point for mean PR-AUC
    plt.plot(1, float(mean_pr_auc), 'o', color=color_dict[hm], label='Mean PR-AUC = {} ± {}'.format(mean_pr_auc, pr_auc_ci))

    # Set plot properties
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve - {hm}')
    plt.legend(loc='lower right')

    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.savefig(output_dir + f'{hm}_PR_curve.png')


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
        
            features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
            features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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
    plt.title(f'SHAP Summary Plot for {hm}', fontsize=10)

    # # Display the plot
    plt.savefig(f'{hm_dir}/{hm}_SHAP.png')
    plt.close()

    mean_abs_shap_values = np.abs(combined_shap_values).mean(axis=0)


    # Plot SHAP dependence plots of all features
    for feature in sf:
        feature_index = sf.index(feature)
        shap.dependence_plot(feature_index, combined_shap_values, combined_X_test, feature_names=sf, show=False, interaction_index=None)

        # to use custom color 
        ax = plt.gca()
        for artist in ax.collections: # remove default scatter plots
            artist.remove()
        plt.scatter(combined_X_test[:,feature_index], combined_shap_values[:, feature_index], s=3, color=color_dict[hm])

        # plt.title(f'{feature} Binding in {hm}-marked Exon Flanks', fontsize=10)
        plt.xlabel(f'Predicted Binding Scores of {feature}', fontsize=10)
        plt.ylabel(f'SHAP Values for {feature}', fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout()
        plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.2)

        #Display the plot
        plt.savefig(f'{hm_dir}/{feature}_{hm}.png')
        plt.close()

    # Create a DataFrame for feature importance
    feature_importance_df = pd.DataFrame({
        'Feature': sf_data.columns, 
        'Mean Absolute SHAP Value': mean_abs_shap_values
    })

    # Sort the features by importance
    feature_importance_df = feature_importance_df.sort_values(by='Mean Absolute SHAP Value', ascending=False)
    feature_importance_df.to_csv(f'0_Files/Post-processing/Analyses/SHAP/{hm}_shap_feature_importance.csv', sep='\t', index=False)


def SHAP_imptRBPs_plot(hm):

    #get SHAP epi and nonepiRBPs to plot
    features_to_plot = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPS/SHAP_epiRBPs/rbps_{hm}.txt", "r")
    features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPS/SHAP_nonepiRBPs/rbps_{hm}.txt", "r")
    features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])

    # get RBPs corr with SHAP RBPs to plot
    corr_features_to_plot = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPS/epiRBPs_{hm}.txt", "r")
    corr_features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPS/nonepiRBPs_{hm}.txt", "r")
    corr_features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    features_to_plot = list(set(features_to_plot + corr_features_to_plot))

    hm_dir = f"0_Files/Post-processing/imptRBPS"
    os.makedirs(hm_dir, exist_ok=True)  # Create the directory if it doesn't exist

    features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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
    plt.savefig(f'{hm_dir}/{hm}_SHAP_impt.png')
    plt.close()


def correlated_to_shap_RBPs(hm):

    print(hm)

    for mode in ['epi', 'nonepi']:

        # Get episplicing RBPs
        rbps_file = open(f"0_Files/Post-processing/{mode}RBPS/SHAP_{mode}RBPs/rbps_{hm}.txt", "r")
        shap_rbps = [rbp for rbp in rbps_file.read().split('\n') if rbp]

        features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
        features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
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
        with open(f"0_Files/Post-processing/{mode}RBPS/{mode}RBPs_{hm}.txt", 'w') as f:
            for line in final_rbps:
                f.write("%s\n" % line)

        print(len(shap_rbps), len(final_rbps))


def imptRBPS_overlap():

    for mode in ['epi', 'nonepi']:
        hm_rbps = {}
        hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
        for hm in hms:
            rbps_file = open(f"0_Files/Post-processing/{mode}RBPS/{mode}RBPs_{hm}.txt", "r")
            hm_rbps[hm] = [rbp for rbp in rbps_file.read().split('\n') if rbp]

        
        all_rbps = list(set([rbp for elem_list in hm_rbps.values() for rbp in elem_list]))
        all_rbps.sort()

        # Get occurence of epiRBP across all HMs

        # Open a TSV file to write the results
        with open(f"0_Files/Post-processing/{mode}RBPS/{mode}RBPS.tsv", 'w', newline='') as tsvfile:
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

    all_prauc_together() # Figure 2a
    metrics() # Figure 2b
    imptRBPS_overlap() # Table3,4

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']

    for hm in hms:
        stratified_hms_classifier(output_dir=sys.argv[1], hm=hm) # indiv PRAUC plots
        confusion_matrix_plot(hm) # Fig 13 (Suppl)
        SHAP(hm=hm) 
        ## Manually select imptRBPs (Epi and nonepi)
        SHAP_imptRBPs_plot(hm) # Fig 3,15
        correlated_to_shap_RBPs(hm=hm) # Fig 16 (Suppl)  
    