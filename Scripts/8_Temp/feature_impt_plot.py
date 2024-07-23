import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys


def viusalise(hm, op_dir):

    op_dir = os.path.join(op_dir, 'feature_impt')
    os.makedirs(op_dir, exist_ok=True)  # Create the directory if it doesn't exist


    # STEP 1: Get permutation importance scores
    sf_sorted = open(f'0_Files/impt_features/sfs_{hm}.txt','r').read().split('\n')
    mean_importance_sorted = open(f'0_Files/impt_features/scores_{hm}.txt','r').read().split('\n')

    perm_impt = dict(zip(sf_sorted, mean_importance_sorted))
    perm_impt = {key: value for key, value in perm_impt.items() if key and value is not ''}
    perm_impt = {key: perm_impt[key] for key in sorted(perm_impt)}
    perm_impt = {key: round(float(value), 2) for key, value in perm_impt.items()}
    
    # STEP 2: Get indiv RBP PRAUC scores
    with open(f'0_Files/Post-processing/indiv_rbp_aucs/{hm}_aucs.json') as f:
            indiv_auc= json.load(f)

    indiv_auc = {key: indiv_auc[key] for key in sorted(indiv_auc)}
    indiv_auc = {key: round(float(value), 2) for key, value in indiv_auc.items()}

    df = pd.DataFrame({'AUC': indiv_auc, 'LOSS': perm_impt}, index=perm_impt.keys())

    # STEP 3: Plot both scores as a heatmap
    heatmap(df,hm,op_dir )


def heatmap(df, hm, op_dir):
     # Plot settings
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot heatmap
    sns.heatmap(df.T, annot=False, cmap='coolwarm', ax=ax, linewidths=0.01, linecolor='white')

    # Labeling
    ax.set_xticks(range(132))
    ax.set_xticklabels(df.index, rotation=90, ha='right', fontsize=5)
    ax.set_xlabel('RBPs')
    ax.set_ylabel('Metrics')
    ax.set_title(f'Feature Importance and Dependence - {hm}')

    # Show plot
    plt.savefig(f'{op_dir}/{hm}.png')
    plt.close()


if __name__ == "__main__":

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3', 'H3K4me1']
    for hm in hms:
        viusalise(hm, op_dir=sys.argv[1])