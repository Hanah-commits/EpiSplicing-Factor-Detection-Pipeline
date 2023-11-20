import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import json


def kde_MAJIQ():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    for hm in hms:

        df1 = pd.read_csv('0_Files/Filtered_MValues.csv', delimiter='\t')
        df2 = pd.read_csv('0_Files/Filtered_MValues_control.csv', delimiter='\t')
        
        data1 = df1[hm].abs().values.tolist()
        data2 = df2[hm].abs().values.tolist()

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
        axs[0].set_xlabel(f'{hm} peak scores near DJU events (Alternative)')
        axs[0].set_ylabel('Density')
        

        # Plot the KDE for data2 on the second subplot
        axs[1].plot(x2, kde_values2)
        axs[1].set_xlabel(f'{hm} peaks scores near non-DJU events (constitutive)')
        axs[1].set_ylabel('Density')
        

        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()


def kde_MAJIQ_epigenes():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    for hm in hms:

        dPSI = pd.read_csv('0_Files/Filtered_dPSI_control.csv', delimiter='\t')
        peaks = pd.read_csv('0_Files/Filtered_MValues_control.csv', delimiter='\t')

        dPSI.drop_duplicates(inplace=True)
        peaks.drop_duplicates(inplace=True)
        control_flanks = pd.merge(dPSI, peaks, how="outer")

        dju_file = '0_Files/dPSI_Mval_epi_' + hm + '.csv'
        dju_hm_flanks = pd.read_csv(dju_file, delimiter='\t')
        dju_genes = list(set(dju_hm_flanks['gene_id'].values.tolist()))
        nondju_hm_flanks = control_flanks[control_flanks['gene_id'].isin(dju_genes)]

        #add label
        dju_hm_flanks['type'] = "dju"
        nondju_hm_flanks['type'] = "non-dju"

        # impute missing data points
        dju_hm_flanks.fillna(0,inplace=True)
        nondju_hm_flanks.fillna(0, inplace=True)
            
        data1 = dju_hm_flanks[hm].abs().values.tolist()
        data2 = nondju_hm_flanks[hm].abs().values.tolist()

        print(data1)

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
        axs[0].set_xlabel(f'{hm} peak scores near DJU events of epigenes (Alternative)')
        axs[0].set_ylabel('Density')
        

        # Plot the KDE for data2 on the second subplot
        axs[1].plot(x2, kde_values2)
        axs[1].set_xlabel(f'{hm} peaks scores near non-DJU events of epigenes(constitutive)')
        axs[1].set_ylabel('Density')
        

        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()


def kde_rmats():

    with open('paths.json') as f:
        d = json.load(f)

    hms = d["Histone modifications"]

    for hm in hms:

        dPSI = pd.read_csv('0_Files/Filtered_dPSI_AS.csv', delimiter='\t')
        peaks = pd.read_csv('0_Files/Filtered_MValues_AS.csv', delimiter='\t')

        dPSI.drop_duplicates(inplace=True)
        peaks.drop_duplicates(inplace=True)
        control_flanks = pd.merge(dPSI, peaks, how="outer")

        dju_file = '0_Files/dPSI_Mval_epi_' + hm + '.csv'
        dju_hm_flanks = pd.read_csv(dju_file, delimiter='\t')
        dju_genes = list(set(dju_hm_flanks['gene_id'].values.tolist()))
        nondju_hm_flanks = control_flanks[control_flanks['gene_id'].isin(dju_genes)]

        #add label
        dju_hm_flanks['type'] = "dju"
        nondju_hm_flanks['type'] = "non-dju"

        # impute missing data points
        dju_hm_flanks.fillna(0,inplace=True)
        nondju_hm_flanks.fillna(0, inplace=True)
            
        data1 = dju_hm_flanks[hm].abs().values.tolist()
        data2 = nondju_hm_flanks[hm].abs().values.tolist()

        print(data1)

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
        axs[0].set_xlabel(f'{hm} peak scores near DJU events of epigenes (Alternative)')
        axs[0].set_ylabel('Density')
        

        # Plot the KDE for data2 on the second subplot
        axs[1].plot(x2, kde_values2)
        axs[1].set_xlabel(f'{hm} peaks scores near non-DJU events of epigenes(constitutive)')
        axs[1].set_ylabel('Density')
        

        # Adjust spacing between subplots
        plt.tight_layout()

        # Show the plots
        plt.show()


kde_rmats()