import pandas as pd
import numpy as np
from pandas import Series
import sys
import json


def adjust_pvalue(df, col):

    # get indices of null values
    na_idx = df[df[col].isnull()].index.tolist()

    # adjust non-null p values
    pvals = df[col].values.tolist()
    pvals = [x for x in pvals if str(x) != 'nan']
    adj_pval = p_adjust_bh(pvals).tolist()

    # insert null at original indices
    for idx in na_idx:
        adj_pval.insert(idx, None)

    # adjusted p values as new df
    df['adj_pval'] = adj_pval
    return df


def p_adjust_bh(p):
    ## multiple hypothesis testing
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
with open('paths.json') as f:
    d = json.load(f)
name = d["tissue1"] + d["tissue2"]
file = sys.argv[1]+'MAJIQ/deltapsi/' + name + '.deltapsi.tsv'
voila = pd.read_csv(file, delimiter='\t')
col_list = ['gene_id', 'lsv_id', 'mean_dpsi_per_lsv_junction', 'probability_changing', 'junctions_coords', 'num_exons'] #, 'exons_coords']
voila = voila[col_list]

# FILTER 1: filter erroneous LSVs
voila["num_exons"] = voila["num_exons"].replace('na' ,'0')
voila["num_exons"] = pd.to_numeric(voila["num_exons"])
voila = voila[voila['num_exons'] > 0]

# split column values to multiple lines
voila = voila.assign(junctions_coords=voila['junctions_coords'].str.split(';'), mean_dpsi_per_lsv_junction=voila['mean_dpsi_per_lsv_junction'].str.split(';'), probability_changing=voila['probability_changing'].str.split(';'))

# explode the list in the columns to create individual rows
voila = voila.explode(['junctions_coords', 'mean_dpsi_per_lsv_junction', 'probability_changing']).reset_index(drop=True)

voila = voila[col_list]
# skipping nan -> 99464127-nan
voila = voila[~voila['junctions_coords'].str.contains("nan")]

# FILTER 2: Drop rows with changing probability < 0.05. (remove DJU events)
voila['pval'] = 1- pd.to_numeric(voila['probability_changing'])
voila = adjust_pvalue(voila, col='pval')
voila = voila[pd.to_numeric(voila['adj_pval']) > 0.05]

# STEP 3: Get the source and target indices of splice junctions
voila[['source', 'target']] = voila['junctions_coords'].str.split('-', n=1, expand=True)

del(voila['junctions_coords'])
voila_temp = voila.copy()
del(voila_temp['target'])
del(voila['source'])

voila.rename(columns={'target': 'junction0'}, inplace = True)
voila_temp.rename(columns={'source': 'junction0'}, inplace=True)

voila = pd.concat([voila_temp, voila]).sort_index(kind='merge')
voila.to_csv('0_Files/majiq_junctions_control.csv', index=False, sep='\t', header=True)