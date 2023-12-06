import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def kde_MAJIQ():
    col = 'mean_dpsi_per_lsv_junction'

    df1 = pd.read_csv('0_Files/Filtered_dPSI.csv', delimiter='\t')
    df2 = pd.read_csv('0_Files/Filtered_dPSI_control.csv', delimiter='\t')

    print(len(df1), len(df2))
    
    data1 = df1[col].abs().values.tolist()
    data2 = df2[col].abs().values.tolist()

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

def kde_RMATS():
    col = 'dPSI'

    df1 = pd.read_csv('0_Files/rmats_flanks200.bed', delimiter='\t')
    df1.columns = ['chr', "flank_start", "flank_stop", "feature", "score", "strand", "geneSymbol", "dPSI"]
    df2 = df1[df1.dPSI != 0.0]

    print(len(df1), len(df2))
    
    data1 = df1[col].abs().values.tolist()
    data2 = df2[col].abs().values.tolist()

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
    axs[0].set_xlabel('| dPSI | values of DEU and non-DEU events')
    axs[0].set_ylabel('Density')
    

    # Plot the KDE for data2 on the second subplot
    axs[1].plot(x2, kde_values2)
    axs[1].set_xlabel('| dPSI | values of only DEU events')
    axs[1].set_ylabel('Density')
    

    # Adjust spacing between subplots
    plt.tight_layout()

    # Show the plots
    plt.show()


kde_RMATS()