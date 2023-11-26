import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde



# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
file = '/Users/hanah/EpiSplicing_RMATS/Output/neuro-H1/RMATS/SE.MATS.JC.txt'
rmats = pd.read_csv(file, delimiter='\t')
col_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'IncLevelDifference', 'FDR', 'exonStart_0base', 'exonEnd']
rmats = rmats[col_list]

# use | dPSI | and only true values
rmats['IncLevelDifference'] = rmats['IncLevelDifference'].abs()
rmats = rmats[rmats['FDR'] <=0.05]

# FILTER 1: Get AS ( |dPSI| > 0.2, FDR < 0.05) and CS exons ( |dPSI| < 0.2, FDR < 0.05)
rmats_AS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() >= 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]
rmats_CS = rmats[(pd.to_numeric(rmats['IncLevelDifference']).abs() < 0.2) & (pd.to_numeric(rmats['FDR']) <= 0.05)]

# FILTER 2: If skipped exon is reported many times,  pick single dPSI score (can happen if down/upstream exons vary)

## get the largest dPSI value for AS exons (most differentially used score)
for i, df in enumerate([rmats_AS, rmats_CS]):
    # Create 'dPSI' and 'dPSI' columns
    df['dPSI'] = df.groupby('exonStart_0base')['IncLevelDifference'].transform(lambda x: ','.join(x.astype(str)))
    df['dPSI'] = df['dPSI'].str.split(',').apply(lambda x: max(map(float, x)) if x[0] else None)

    # Keep only rows where 'IncLevelDifference' is equal to 'dPSI'
    df = df[df['IncLevelDifference'] == df['dPSI']]

    # FILTER 3: Drop duplicate exon entries
    df = df.drop_duplicates(subset=["geneSymbol", "strand", "exonStart_0base", "exonEnd"], keep='first')

    # Assign the modified DataFrame back to the original variable
    if i == 0:
        rmats_AS = df
    else:
        rmats_CS = df

# FILTER 4: Get only the exons from rmats_CS that are unique to it (not in rmats_AS)
merged_df = pd.merge(rmats_CS, rmats_AS[["geneSymbol", "strand", "exonStart_0base", "exonEnd"]], on=["geneSymbol", "strand", "exonStart_0base", "exonEnd"], how='left', indicator=True)
rmats_CS = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])

# FILTER 5: Keep coords of single version of exon if A3SS/A5SS events exist (to prevent 2+ flanks per exon)

# A5SS
rmats_AS.sort_values(by=['exonStart_0base', 'dPSI'], ascending=[True, False], inplace=True)
rmats_AS.drop_duplicates(subset=['exonStart_0base'], keep='first', inplace=True)
rmats_CS.sort_values(by=['exonStart_0base', 'dPSI'], ascending=[True, False], inplace=True)
rmats_CS.drop_duplicates(subset=['exonStart_0base'], keep='first', inplace=True)

# A3SS
rmats_AS.sort_values(by=['exonEnd', 'dPSI'], ascending=[True, False], inplace=True)
rmats_AS.drop_duplicates(subset=['exonEnd'], keep='first', inplace=True)
rmats_CS.sort_values(by=['exonEnd', 'dPSI'], ascending=[True, False], inplace=True)
rmats_CS.drop_duplicates(subset=['exonEnd'], keep='first', inplace=True)


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
    df_bed.to_csv(f'0_Files/SE_{type[i]}.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
    df.to_csv(f'0_Files/SE_exons_{type[i]}.csv', index=False, sep='\t', header=True)

    
    










def view_dpsi(type):

    if type == "sub":
        data1 = rmats_AS['IncLevelDifference'].abs().values.tolist()
        data2 = rmats_CS['IncLevelDifference'].abs().values.tolist()

        kde1 = gaussian_kde(data1)
        kde2 = gaussian_kde(data2)

        # Generate points on the x-axis for the KDE plots
        x1 = np.linspace(min(data1), max(data1), 1000)
        x2 = np.linspace(min(data2), max(data2), 1000)

        # Calculate the KDE values for both data lists
        kde_values1 = kde1(x1)
        kde_values2 = kde2(x2)

        # Create two subplots side by side
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        # Plot the KDE for data1 on the first subplot
        axs[0].plot(x1, kde_values1)
        axs[0].set_xlabel('| dPSI | values of DJU events (Alternative)')
        axs[0].set_ylabel('Density')


        # Plot the KDE for data2 on the second subplot
        axs[1].plot(x2, kde_values2)
        axs[1].set_xlabel('| dPSI | values of non-DJU events (constitutive)')
        axs[1].set_ylabel('Density')


        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()

    elif type == "main":
        data = rmats['IncLevelDifference'].abs().values.tolist()
        kde = gaussian_kde(data)
       

        # Generate points on the x-axis for the KDE plots
        x1 = np.linspace(min(data), max(data), 1000)

        # Calculate the KDE values for both data lists
        kde_values1 = kde(x1)

        # Create two subplots side by side
        fig, axs = plt.subplots(1, figsize=(12, 5))

        # Plot the KDE for data1 on the first subplot
        axs.plot(x1, kde_values1)
        axs.set_xlabel('| dPSI | values of SE events')
        axs.set_ylabel('Density')

        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()
