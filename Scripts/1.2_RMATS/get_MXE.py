import pandas as pd
import numpy as np
import sys



# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
file = sys.argv[1]+ 'RMATS/MXE.MATS.JC.txt'
rmats = pd.read_csv(file, delimiter='\t')

col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd']
rmats = rmats[col_list]

print('Processing RMATS output: Mutually Exclusive Exons \n')
print('# genes reported:                ', len(set(rmats.geneSymbol.values.tolist()))) # log

# use | dPSI | and only true values
rmats['IncLevelDifference'] = rmats['IncLevelDifference'].abs()
rmats = rmats[rmats['FDR'] <=0.05]

print('FDR-adj pvalue <= 0.05:          ', len(set(rmats.geneSymbol.values.tolist()))) # log

if len(rmats) == 0:
    print(' No mutually exclusive exons to process \n')
    sys.exit(0)


# STEP 2 : Split into multiple rows, keeping one exon coord in one row.

# Create two DataFrames, one for each row
row1 = rmats[['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', '1stExonStart_0base', '1stExonEnd']]
row2 = rmats[['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', '2ndExonStart_0base', '2ndExonEnd']]

# Rename columns 
row1.columns = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
row2.columns = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']

# mark exon order
row1['exon_order'] = 1
row2['exon_order'] = 2

# Concatenate the two DataFrames to get the output
rmats = pd.concat([row1, row2], ignore_index=True)

# # housekeeping
# os.system('rm 0_Files/RMATS/rmats_*.bed')

# STEP 3: Get dPSI scores based on inclusion exon
# NOTE: the inclusion isoform includes the exon that is “earlier” in the transcript.

rmats['dPSI'] = np.where(
    (rmats['strand'] == '+') & (rmats['exon_order'] == 1),
    rmats['IncLevelDifference'],
    np.where(
        (rmats['strand'] == '+') & (rmats['exon_order'] == 2),
        1 - rmats['IncLevelDifference'],
        np.where(
            (rmats['strand'] == '-') & (rmats['exon_order'] == 2),
            rmats['IncLevelDifference'],
            np.where(
                (rmats['strand'] == '-') & (rmats['exon_order'] == 1),
                1 - rmats['IncLevelDifference'],
                np.nan  # default value for other cases
            )
        )
    )
)

# FILTER 1: Get true MXE events
try:
    SE_exons = pd.read_csv("0_Files/RMATS/SE_exons.csv", delimiter='\t')
    SE_exons = list(set(SE_exons[SE_exons.dPSI > 0.2].exon_coord0.values.tolist()))
    rmats = rmats[(~rmats['exonStart_0base'].isin(SE_exons)) & (~rmats['exonEnd'].isin(SE_exons))] # covers A3SS,A5SS versions of skipped exons

    print('True MXE:                        ', len(set(rmats.geneSymbol.values.tolist()))) # log
except:
    print('No skipped exons available: All are considered True MXE')

# FILTER 2: Get AS ( |dPSI| > 0.2, FDR < 0.05) and CS exons ( |dPSI| < 0.2, FDR < 0.05)
rmats_AS = rmats[(pd.to_numeric(rmats['dPSI'] >= 0.2)) & (pd.to_numeric(rmats['FDR']) <= 0.05)]

# FILTER 3: Keep coords of single version of exon if A3SS/A5SS events exist (to prevent 2+ flanks per exon)

# #                   GeneID geneSymbol    chr strand  IncLevelDifference       FDR  exonStart_0base   exonEnd   dPSI
# # 4634  ENSG00000126456.15       IRF3  chr19      -               0.449  0.000202         49664442  49664673  0.449
# # 4636  ENSG00000126456.15       IRF3  chr19      -               0.584  0.000012         49664552  49664673  0.584
# # 4640  ENSG00000126456.15       IRF3  chr19      -               0.363  0.015466         49664586  49664673  0.363

def A3SS_A5SS_filter(group, subset_column):
    group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
    group.drop_duplicates(subset=['exonStart_0base'], keep='first', inplace=True)
    group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
    group.drop_duplicates(subset=['exonEnd'], keep='first', inplace=True)
    return group

rmats_AS = rmats_AS.groupby('GeneID').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
rmats_AS.reset_index(drop=True, inplace=True)

print('| IncLevelDifference | > 0.2:    ', len(set(rmats_AS.geneSymbol.values.tolist()))) # log

df = rmats_AS.copy()
# temp output fiilee
df['feature'] = "Exon"
df['score'] = "."
df['exonStart_0base'] = pd.to_numeric(df['exonStart_0base']) + 1
df[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "GeneID", "dPSI"]].to_csv(f'0_Files/RMATS/MXE_exons.tsv', index=False, sep='\t', header=True)
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
df_bed.to_csv(f'0_Files/RMATS/MXE.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
df.to_csv(f'0_Files/RMATS/MXE_exons.csv', index=False, sep='\t', header=True)

