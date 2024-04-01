import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd

def kde_exonlength():
    col = 'exon_bp'

    df1 = pd.read_csv('0_Files/Post-processing/epi_exons.bed', delimiter='\t', header=None)
    df2 = pd.read_csv('0_Files/Post-processing/nonepi_exons.bed', delimiter='\t', header=None)
    
    for df in [df1, df2]:
        df.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']
        df[col] = df['exon_end'] - df['exon_start']
    
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
    axs[0].set_xlabel('Lengths of alternative exons - Epigenes')
    axs[0].set_ylabel('Density')

    # Plot the KDE for data2 on the second subplot
    axs[1].plot(x2, kde_values2)
    axs[1].set_xlabel('Lengths of alternative exons - Non-epigenes')
    axs[1].set_ylabel('Density')

    # Adjust y-axis limits to include only non-zero densities
    axs[0].set_ylim(0, max(kde_values1) * 1.1)
    axs[1].set_ylim(0, max(kde_values2) * 1.1)

    # Adjust x-axis limits to include only non-zero densities
    axs[0].set_xlim(0, max(x1) * 1.1)
    axs[1].set_xlim(0, max(x2) * 1.1)

    # Add a vertical line at the maximum x-axis value for each plot
    axs[0].axvline(x=max(data1), color='red', linestyle='--', label=max(data1))
    axs[1].axvline(x=max(data2), color='red', linestyle='--', label=max(data2))

    # Adjust spacing between subplots
    plt.tight_layout()

    # Show the plots
    # plt.savefig('exon_lengths.png')
    plt.show()

# kde_exonlength()

def percentage_histogram(data, bins):
    hist, _ = np.histogram(data, bins=bins)
    return hist / len(data) * 100

def histogram_exonlength():
    col = 'exon_bp'

    df1 = pd.read_csv('0_Files/Post-processing/epi_exons.bed', delimiter='\t', header=None)
    df2 = pd.read_csv('0_Files/Post-processing/nonepi_exons.bed', delimiter='\t', header=None)
    
    for df in [df1, df2]:
        df.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']
        df[col] = df['exon_end'] - df['exon_start']
    
    data1 = df1[col].abs().values
    data2 = df2[col].abs().values

    bins = np.linspace(0, max(max(data1), max(data2)), 50)

    hist1 = percentage_histogram(data1, bins)
    hist2 = percentage_histogram(data2, bins)

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    axs[0].bar(bins[:-1], hist1, width=bins[1] - bins[0], align='edge')
    axs[0].set_xlabel('Lengths of Alternative Exons - Epigenes (bp) ')
    axs[0].set_ylabel('Frequency of Occurrence (%)')

    axs[1].bar(bins[:-1], hist2, width=bins[1] - bins[0], align='edge')
    axs[1].set_xlabel('Lengths of Alternative Exons - Non-epigenes (bp)')
    axs[1].set_ylabel('Frequency of Occurrence (%)')

    plt.tight_layout()
    plt.show()

histogram_exonlength()