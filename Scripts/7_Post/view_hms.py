import pandas as pd
import matplotlib.pyplot as plt
from venn import venn


def view_hm_venn():
    # Read data from files
    df1 = pd.read_csv('0_Files/Post-processing/epi_exons.bed', delimiter='\t', header=None)
    df2 = pd.read_csv('0_Files/Post-processing/nonepi_exons.bed', delimiter='\t', header=None)
    
    i = 0
    for df in [df1, df2]:
        # Determine dataset name based on index
        name = 'Epi' if i == 0 else 'Nonepi'
        
        # Assign column names
        df.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'type', 'strand', 'gene_name']

        # Split column values into multiple lines
        df = df.assign(type=df['type'].str.split(','))

        # Explode the list in the columns to create individual rows
        df = df.explode('type').reset_index(drop=True)

        # Get unique HM types and their associated genes
        hms = df['type'].unique()
        hm_genes = {hm: set(df[df['type'] == hm]['gene_name']) for hm in hms}

        # Create Venn diagram data
        data = {str(hm): genes for hm, genes in hm_genes.items()}

        # Plot Venn diagram
        venn_data = {key: value for key, value in data.items() if len(value) > 0}
        venn(venn_data, cmap='plasma')
        plt.title(f"{len(df['gene_name'].unique())} {name}genes")
        plt.savefig(f'0_Files/Post-processing/hm_overlap_{name}.png')
        plt.close()
        
        i += 1

view_hm_venn()