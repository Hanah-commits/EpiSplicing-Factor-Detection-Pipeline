import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde



# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
file = '/Users/hanah/EpiSplicing_RMATS/Output/neuro-H1/RMATS/A3SS.MATS.JCEC.txt'
rmats = pd.read_csv(file, delimiter='\t')

col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE' ]
rmats = rmats[col_list]

# STEP 2 : Split into multiple rows, keeping one exon coord in one row.

# Create two DataFrames, one for each row
row1 = rmats[['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'longExonStart_0base', 'longExonEnd']]
row2 = rmats[['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'shortES', 'shortEE']]

# Rename columns 
row1.columns = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
row2.columns = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']

# Concatenate the two DataFrames to get the output
rmats = pd.concat([row1, row2], ignore_index=True)

# FILTER 1: Keep single set of exon-coords for every A3SS event

# same event mentioned more than once
rmats = rmats.drop_duplicates()
# same exon mentioned more than once
rmats = rmats.drop_duplicates(subset=['exonStart_0base'])

# Reset the index
rmats = rmats.reset_index(drop=True)

# FILTER 1: Get AS ( |dPSI| > 0.2, FDR < 0.05) and CS exons ( |dPSI| < 0.2, FDR < 0.05)
rmats_AS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() >= 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]
rmats_CS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() < 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]


# STEP 2: Prepare bedtools input

dfs = [rmats_AS, rmats_CS]
type = ['AS', 'CS']
for i in range(0,2):

    df = dfs[i].copy()
    df_temp = df.copy()
    del(df_temp['exonStart_0base'])
    del(df['exonEnd'])

    df.rename(columns={'exonStart_0base': 'exon_coord0'}, inplace = True)
    df_temp.rename(columns={'exonEnd': 'exon_coord0'}, inplace=True)

    df = pd.concat([df_temp, df]).sort_index(kind='merge')

    keep_cols = ['chr', 'exon_coord0', 'strand']
    df_bed = df[keep_cols]
    df_bed = df_bed.drop_duplicates()
    # to fit bedtools input requirements
    df_bed['exon_coord1'] = pd.to_numeric(df_bed['exon_coord0']) + 1
    df_bed['feature'] = "flank"
    df_bed['score'] = "."
    

    df_bed = df_bed[['chr', "exon_coord0", "exon_coord1", "feature", "score", "strand"]]
    df_bed.to_csv(f'0_Files/A3SS_{type[i]}.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
    df.to_csv(f'0_Files/A3SS_exons_{type[i]}.csv', index=False, sep='\t', header=True)
