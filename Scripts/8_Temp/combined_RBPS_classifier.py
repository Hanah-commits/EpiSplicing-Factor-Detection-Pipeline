import pandas as pd
import numpy as np
import sys, os, csv
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate, GridSearchCV, train_test_split
from sklearn.metrics import precision_recall_curve, average_precision_score, make_scorer, fbeta_score
from sklearn.metrics import confusion_matrix
from sklearn.utils import shuffle
import shap
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests



def untuned_vs_tuned():

    
    hms= [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    eval_scores = {}
    os.makedirs("0_Files/Post-processing/Analyses/Metrics/", exist_ok=True)  # Create the directory if it doesn't exist

    for hm in hms:

        print('\n\n',hm, '\n')

        features1 = pd.read_csv('0_Files/Post-processing/features_all_132.csv', delimiter='\t')
        features2 = pd.read_csv('0_Files/Post-processing/features_all_47.csv', delimiter='\t')
        features = pd.concat([features1,features2], axis=1)
        features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
        features.fillna(0,inplace=True)
        features = shuffle(features, random_state=42)

        # extract features for hm
        features = features[features['type'].apply(lambda x: any(item in [hm] for item in x.split(',')))]
        features = features.drop('type', axis=1) # drop hm info

        # keep strong binding events
        sf = [val for val in features.columns if val != 'label']

        # keep only strong binding events
        sf_data = features[sf]
        sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

        features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
        X, y = sf_data, features['label']

        # Train-test split: 20% test set
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)
        
        ## STEP 1: Initialize untuned classifier
        
        # cross validation method
        kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10) 

        # Step 1: Train-test split (70%-30%)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)

        modes = ['untuned', 'tuned'] 

        for mode in modes:

            if mode == 'untuned':

                untuned_clf = RandomForestClassifier(class_weight='balanced', bootstrap=False, n_jobs = -1, random_state=12)

                # Train untuned model
                pr_aucs = []
                for i, (train, test) in enumerate(kf.split(X_train, y_train)):

                    X_train_fold = X_train.iloc[train]
                    y_train_fold = y_train.iloc[train]
                    X_test_fold = X_train.iloc[test]
                    y_test_fold = y_train.iloc[test]

                    model = untuned_clf.fit(X_train_fold, y_train_fold)
                    y_score = model.predict_proba(X_test_fold)
                    avg_precision = average_precision_score(y_test_fold, y_score[:, 1])
                    pr_aucs.append(avg_precision)


                # Training performance
                pr_auc_range = np.percentile(pr_aucs, (2.5, 97.5))
                pr_auc_ci = float("%.2f" % (pr_auc_range[1] - pr_auc_range[0])) / 2
                untuned_train_score = (pr_auc_range[1] - pr_auc_ci)
                print(f'Performance of untuned model (training): {untuned_train_score:.4f} \n')

                #  Test the untuned model
                y_pred_untuned = untuned_clf.predict_proba(X_test)[:, 1]
                untuned_test_score = average_precision_score(y_test, y_pred_untuned)
                print(f"Performance of untuned model (test): {untuned_test_score:.4f} \n")

            elif mode == 'tuned':

                ## Train tuned model (all features)
                param_grid = {
                "n_estimators": [100, 200],
                "max_features": ["sqrt", "log2", None],
                "max_depth": [3, 5, None],
                "min_samples_split" : [2, 5]
                }
                
                scorer = "average_precision"

                tuned_clf = GridSearchCV(
                estimator=RandomForestClassifier(class_weight='balanced', bootstrap=False, n_jobs = -1, random_state=12),
                param_grid=param_grid,
                scoring=scorer,
                cv=kf,
                n_jobs=-1,
                verbose=1,
                )

                # Fit GridSearchCV on the training set
                tuned_clf.fit(X_train, y_train)
                print(f"Best hyperparameters: {tuned_clf.best_params_} \n")
                tuned_train_score = tuned_clf.best_score_
                print(f"Performance of tuned model (training): {tuned_train_score:.4f}")

                # Test the tuned model
                best_model = tuned_clf.best_estimator_
                y_pred_tuned = best_model.predict_proba(X_test)[:, 1]
                tuned_test_score = average_precision_score(y_test, y_pred_tuned)
                print(f"Performance of tuned model: {tuned_test_score:.4f} \n")
                

        # save scores
        eval_scores[hm] = [untuned_train_score, untuned_test_score, tuned_train_score, tuned_test_score]

    # print(eval_scores)
    # eval_scores = {'H3K27ac': [0.8725345520790744, 0.8342091060071373, 0.8870289046935795, 0.9066179046212118], 'H3K27me3': [0.7936426767676767, 0.753075560187311, 0.8370049747882313, 0.8510441809190956], 'H3K36me3': [0.80523916702003, 0.8115620890083819, 0.8165537169654283, 0.8613831726744965], 'H3K9me3': [0.9045557190451009, 0.9338054364686202, 0.9360760624546037, 0.9366370403195617], 'H3K4me3': [0.8054078360910345, 0.8934364712567686, 0.8353645629288127, 0.9298684806477222]}
    plot_evaluation_scores(eval_scores=eval_scores)


def plot_evaluation_scores(eval_scores):
    ##Plot evaluation scores

    hms= [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))
    score_types = ['Untuned Training Score', 'Untuned Test Score', 'Tuned Training Score', 'Tuned Test Score']

    # positions for each score type on the x-axis
    score_positions = np.arange(len(score_types))

    jiggle_strength = 0.05 

    # store handles for the legend
    handles = []
    labels = []

    # Plot metric color based on group
    fig, ax = plt.subplots()
    for i, (key, val_list) in enumerate(eval_scores.items()):

        jiggle = np.random.uniform(-jiggle_strength, jiggle_strength, size=len(score_types))
        x_positions = score_positions + jiggle  # Add the jiggle to x positions
        for j, val in enumerate(val_list):
            scatter = ax.scatter(x_positions[j], val, color=color_dict[key], label=key if j == 0 else "", s=20)

        # add a handle for the group if it's the first appearance
        if key not in labels:
            handles.append(scatter)
            labels.append(key)
        
        ax.plot(x_positions, val_list, color=color_dict[key], linestyle='-', linewidth=2)   # Connect the dots with a line

    # customize plot
    ax.set_ylabel('Average Precision', fontsize=8)
    ax.set_title('Evaluation of Untuned vs Tuned Models', fontsize=8)
    ax.set_xticks(score_positions)
    ax.set_xticklabels(score_types, fontsize=6)
    ax.tick_params(axis='y', labelsize=6) 
    ax.legend(handles=handles, labels=labels, fontsize=6)

    # Display the plot
    plt.savefig(f"0_Files/Post-processing/Analyses/Metrics/Untuned_vs_Tuned.png", bbox_inches='tight', dpi=300)
    plt.close()


def all_prauc_together():

    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    parameters = {
        "H3K27ac" : {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200},
        "H3K27me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K36me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K9me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K4me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 5, 'n_estimators': 200}
    }

    # Initialize plot
    plt.figure(figsize=(7, 7))
    plt.axes().set_aspect('equal', 'datalim')

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

        # Classifier and cross-validation
        clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
        kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10) 

        all_precisions = []
        mean_pr_aucs = []

        # Loop over folds
        for i, (train, test) in enumerate(kf.split(X, y)):
            model = clf.fit(X[train], y[train])
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

        # plot the mean curve
        plt.plot(mean_recall, mean_precision, color=color_dict[hm], label=f'{hm} Mean PR-AUC = {mean_pr_auc:.2f} ± {pr_auc_ci:.2f}')

    # customize plot
    plt.xlabel('Recall',  fontsize=10)
    plt.ylabel('Precision',  fontsize=10)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.legend(loc='lower right')

    plt.title(f'Precision-Recall Curve', fontsize=10)
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
    ax.set_ylabel("Mean Scores")
    ax.set_title("Model Metrics", fontsize=10)

    # Create custom legend
  
    pattern_patches = [mpatches.Patch(facecolor='grey', edgecolor='white', hatch=pattern, label=metric) for pattern, metric in zip(patterns, metrics)]
    ax.legend(handles=pattern_patches, title="Metrics", loc='upper left', bbox_to_anchor=(0.7, 0.95), prop={'size': 8}, handleheight=2, handlelength=3)


    plt.tight_layout()
    plt.savefig(f"0_Files/Post-processing/Analyses/Metrics/metrics.png", bbox_inches='tight', dpi=300)


def metrics():

    # Custom scoring function for specificity
    def specificity_score(y_true, y_pred):
        cm = confusion_matrix(y_true, y_pred)
        tn, fp, fn, tp = cm.ravel()
        specificity = tn / (tn + fp)
        return specificity

    # Create a specificity scorer, f-beta scores
    specificity_scorer = make_scorer(specificity_score)
    fbeta_scorer = make_scorer(fbeta_score, beta=0.5)

    parameters = {
        "H3K27ac" : {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200},
        "H3K27me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K36me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K9me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K4me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 5, 'n_estimators': 200}
    }

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
            clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
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


    sf = [val for val in features.columns if val != 'label']

    ## STEP 0: keep only strong binding events
    sf_data = features[sf]
    sf_data = sf_data.applymap(lambda val: 0 if val < 2 else val)

    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    ## STEP 1: Initialize classifier and cross-validation
    clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
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
    sns.heatmap(cm_percentage, annot=True, fmt=".2f", cmap="viridis", cbar_kws={'label': 'Percentage'})

    # Add labels
    plt.xlabel("Predicted Labels")
    plt.ylabel("True Labels")
    plt.title(f'Confusion Matrix - {hm}')
    plt.savefig(f"0_Files/Post-processing/Analyses/Metrics/{hm}_CM.png", bbox_inches='tight', dpi=300)
    plt.close()


def SHAP( hm):


    parameters = {
        "H3K27ac" : {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200},
        "H3K27me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K36me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K9me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K4me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 5, 'n_estimators': 200}
    }

    hms = [  "H3K27ac","H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3"]
    color_dict = dict(zip(hms,["#9A71F8", "#69D4EC", "#B0D212", "#FF9900", "#ED588A"]))

    output_dir = sys.argv[1]
    # Create a directory specific to the hm value
    hm_dir = os.path.join(output_dir, hm)
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
    clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
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

        plt.title(f'{feature} Binding in {hm} Flanks', fontsize=10)
        plt.xlabel(f'Binding Scores of {feature}', fontsize=10)
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
    feature_importance_df.to_csv(output_dir + f'/{hm}_shap_feature_importance.csv', sep='\t', index=False)


def SHAP_imptRBPs_plot(hm):

    parameters = {
        "H3K27ac" : {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200},
        "H3K27me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K36me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K9me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K4me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 5, 'n_estimators': 200}
    }

    #get epi and nonepiRBPs to plot
    features_to_plot = []
    rbps_file = open(f"0_Files/Post-processing/epiRBPS/SHAP_epiRBPs/rbps_{hm}.txt", "r")
    features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])
    rbps_file = open(f"0_Files/Post-processing/nonepiRBPS/SHAP_nonepiRBPs/rbps_{hm}.txt", "r")
    features_to_plot.extend([rbp for rbp in rbps_file.read().split('\n') if rbp])


    output_dir = sys.argv[1]
    # Create a directory specific to the hm value
    hm_dir = os.path.join(output_dir, hm)
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
    clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
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
    shap.summary_plot(filtered_shap_values, filtered_X_test, feature_names=filtered_feature_names, show=False)

    # Use matplotlib to add a title
    plt.title(f'SHAP Summary Plot for {hm}', fontsize=12)

    # # Display the plot
    plt.savefig(f'{hm_dir}/{hm}_SHAP_impt.png')
    plt.close()


def correlated_to_epi_RBPs(hm):

    print(hm)

    # Get episplicing RBPs
    rbps_file = open(f"0_Files/Post-processing/epiRBPS/SHAP_epiRBPs/rbps_{hm}.txt", "r")
    epi_rbps = [rbp for rbp in rbps_file.read().split('\n') if rbp]

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
    epi_rbp_corr = filtered_corr[filtered_corr.RBP1.isin(epi_rbps)]
    final_epi_rbps = list(set(epi_rbp_corr.RBP2.values.tolist() + epi_rbps))
    final_epi_rbps.sort()
    with open(f"0_Files/Post-processing/epiRBPS/epiRBPs_{hm}.txt", 'w') as f:
        for line in final_epi_rbps:
            f.write("%s\n" % line)

    print(len(epi_rbps), len(final_epi_rbps))


def epiRBPS_overlap():

    hm_rbps = {}
    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:
        rbps_file = open(f"0_Files/Post-processing/epiRBPS/epiRBPs_{hm}.txt", "r")
        hm_rbps[hm] = [rbp for rbp in rbps_file.read().split('\n') if rbp]

    
    all_rbps = list(set([rbp for elem_list in hm_rbps.values() for rbp in elem_list]))
    all_rbps.sort()

    # Get occurence of epiRBP across all HMs

     # Open a TSV file to write the results
    with open(f"0_Files/Post-processing/epiRBPS/epiRBPs.tsv", 'w', newline='') as tsvfile:
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


def SHAP_weights(hm):
    parameters = {
        "H3K27ac" : {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200},
        "H3K27me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K36me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K9me3": {'max_depth': None, 'max_features': 'sqrt', 'min_samples_split': 2, 'n_estimators': 200}, 
        "H3K4me3": {'max_depth': None, 'max_features': 'log2', 'min_samples_split': 5, 'n_estimators': 200}
    }


    output_dir = sys.argv[1]
    # Create a directory specific to the hm value
    hm_dir = os.path.join(output_dir, hm)
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
    clf = RandomForestClassifier(n_estimators=parameters[hm]['n_estimators'], max_features= parameters[hm]['max_features'], max_depth=parameters[hm]['max_depth'], min_samples_split= parameters[hm]['min_samples_split'], class_weight='balanced', n_jobs = -1, random_state=0)
    kf = RepeatedStratifiedKFold(n_splits=5, random_state=42, n_repeats=10) 

   # Initialize variables to store cumulative SHAP values and counts
    cumulative_shap_values = np.zeros((X.shape[0], X.shape[1]))
    fold_counts = np.zeros(X.shape[0])

    # Loop over folds
    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        
        # Calculate SHAP values
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X[test])
        
        # Add SHAP values for this fold 
        cumulative_shap_values[test, :] += shap_values[:,:,1]

        # Count occurrences of each test index (for averaging)
        fold_counts[test] += 1

    # Compute mean SHAP values across folds
    mean_shap_values = cumulative_shap_values / fold_counts[:, None]

    shap_df = pd.DataFrame(mean_shap_values, columns=sf)
    

    # Ssve shap values
    shap_df.to_csv(output_dir + f'/{hm}_shap_weights.tsv', sep='\t', index=False)



if __name__ == "__main__":

    untuned_vs_tuned() # Tune parameters
    all_prauc_together() # Figure 2a
    metrics() # Figure 2b
    epiRBPS_overlap() # Table3,4

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']

    for hm in hms:
        stratified_hms_classifier(output_dir=sys.argv[1], hm=hm) # indiv PRAUC plots
        confusion_matrix_plot(hm) # Fig 13 (Suppl)
        SHAP(hm=hm) 
        ## Manually select imptRBPs (Epi and nonepi)
        SHAP_imptRBPs_plot(hm) # Fig 3,15
        correlated_to_epi_RBPs(hm=hm) # Fig 16 (Suppl)  
        # SHAP_weights(hm)