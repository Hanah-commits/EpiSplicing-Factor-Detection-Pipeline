import pandas as pd
import sys
import requests
import json
import numpy as np
from pathlib import Path
import logomaker
import matplotlib.pyplot as plt


def map_gene_to_uniprot(gene_names):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    gene_to_uniprot = {}

    for gene in gene_names:
        # Query UniProt API for the gene name
        params = {
            "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true", # include only UniProtKB reviewed entries (Swiss-Prot)
            "fields": "accession,gene_names",
            "format": "json",
            "size": 1  # Retrieve only one result per gene
        }
        response = requests.get(base_url, params=params)
        
        if response.status_code == 200:
            data = response.json()
            if "results" in data and data["results"]:
                uniprot_id = data["results"][0]["primaryAccession"]
                gene_to_uniprot[uniprot_id] = gene
            else:
                gene_to_uniprot[gene] = None  # No match found
        else:
            print(f"Failed to fetch data for gene: {gene}")
            gene_to_uniprot[gene] = None

    return gene_to_uniprot



def conveert_to_uniprot(dir):

    # Get 160 rbps
    rbps_file = open(f'{dir}/rbps/rbps.txt', 'r')
    rbps = list(set([rbp.strip() for rbp in rbps_file.read().split('\n')]))

    # Convert to uniprot IDs
    rbps = map_gene_to_uniprot(rbps)# convert to uniprot ids
    with open(f'{dir}/rbps/rbps_uniprot.json', 'w', encoding='utf-8') as f:
        json.dump(rbps, f, ensure_ascii=False, indent=4)


    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:

        print('\n\n',hm)

        # Get HM-assocated proteins
        readers_file = open(f'{dir}/readers/{hm}_readers.txt', 'r')
        readers = list(set([reader.strip() for reader in readers_file.read().split('\n')]))
        
        readers = map_gene_to_uniprot(readers)# convert to uniprot ids
        with open(f'{dir}/readers/{hm}_readers_uniprot.json', 'w', encoding='utf-8') as f:
            json.dump(readers, f, ensure_ascii=False, indent=4)



def read_PPI_intact(dir):

    with open(f'{dir}/rbps/rbps_uniprot.json') as f: # read uniprot ids - RBPs
        rbps = json.load(f) 
    

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in ['H3K27ac', 'H3K36me3']:

        print('\n\n',hm)

        # Get HM-assocated proteins
        with open(f'{dir}/readers/{hm}_readers_uniprot.json') as f: # read uniprot ids - RBPs
            readers = json.load(f) 
       

        # Read PPI of current HM readers with RBPs
        ppi = pd.read_csv(f'{dir}/ppi/IntAct/{hm}_PPI.tsv', delimiter='\t')
        relevant_cols = ['# ID(s) interactor A', 'ID(s) interactor B', 'Interaction detection method(s)', 'Publication Identifier(s)', 'Confidence value(s)',  'Host organism(s)']
        ppi = ppi[relevant_cols]
        ppi.columns = ['RBP', 'Reader', 'Experiment_type', 'Experiment_ID', 'MI_Score', 'Host_organism']
    
        ## obtain relevant substrings
        ppi['RBP'] = ppi['RBP'].str.extract(r'uniprotkb:([A-Za-z0-9]+)')
        ppi['Reader'] = ppi['Reader'].str.extract(r'uniprotkb:([A-Za-z0-9]+)')
        ppi['Experiment_type'] = ppi['Experiment_type'].str.split('(').str[1].str.split(')').str[0]
        ppi['Experiment_ID'] = ppi['Experiment_ID'].str.split('pubmed:').str[1].str.split('|').str[0]
        ppi['MI_Score'] = ppi['MI_Score'].str.extract(r'intact-miscore:([0-9]*\.[0-9]+)')
        ppi['Host_organism'] = ppi['Host_organism'].str.split('taxid:').str[1].str.split('(').str[0]

        # Keep H. sapiens and S. cerevisea as host organisms
        ppi = ppi[ppi['Host_organism'].isin(['4932', '9606'])]

        rbp_ids = list(rbps.keys())
        reader_ids = list(readers.keys())

        # Remove PPIs between RBPs, between Readers
        ppi= ppi[(ppi['RBP'].isin(rbp_ids) & ~(ppi['Reader'].isin(rbp_ids)))].drop_duplicates()
        ppi= ppi[(~ppi['RBP'].isin(reader_ids) & (ppi['Reader'].isin(reader_ids)))].drop_duplicates()
        print(len(ppi))
        print(ppi)

        print(len(list(set(ppi['RBP'].values.tolist()))))
       


def read_PPI_string(dir):

    # Get 160 rbps
    rbps_file = open(f'{dir}/rbps/rbps.txt', 'r')
    rbps = list(set([rbp.strip() for rbp in rbps_file.read().split('\n')]))

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:

        print('\n\n',hm)

        # Get HM-assocated proteins
        readers_file = open(f'{dir}/readers/{hm}_readers.txt', 'r')
        readers = list(set([reader.strip() for reader in readers_file.read().split('\n')]))

        # Read PPI network
        ppi = pd.read_csv(f'{dir}/ppi/STRING/{hm}_PPI.tsv', delimiter='\t')[['#node1', 'node2', 'combined_score']]
        
        # Remove PPIs between RBPs, between Readers 
        ppi= ppi[(ppi['#node1'].isin(rbps) & ~(ppi['node2'].isin(rbps)))].drop_duplicates()
        ppi= ppi[(~ppi['#node1'].isin(readers) & (ppi['node2'].isin(readers)))].drop_duplicates()

        ppi.to_csv(f'{dir}/ppi/STRING/{hm}_filtered_PPI.tsv', sep='\t', index=False)

        

def expand_iupac(sequence):

    iupac_codes = {
        'R':['A', 'G'],
        'Y':['C','U'],
        'S':['G','C'],
        'W':['A','U'],
        'K':['G','U'],
        'M':['A','C'],
        'H': ['A','C','U'],
        'V': ['G','C','A'],
        'B': ['G','U','C'],
        'D': ['G','A','U']
    }

    """Expands IUPAC codes in a sequence to all possible nucleotides."""
    expanded = []
    for nucleotide in sequence.upper():
        if nucleotide in iupac_codes:
            expanded.append(iupac_codes[nucleotide])
        else:
            expanded.append([nucleotide])

            
    return expanded



def generate_pwm(sequences):

    """Generates a PWM from a list of sequences."""
    
    # Initialize PWM with zeros
    pwm = np.zeros((len(sequences), 4))  # 4 columns for A, C, G, U
    
    # Nucleotide indices for convenience
    nucleotide_indices = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
    
    # Fill in the PWM
    i = 0
    for seq in sequences:
        if len(seq) == 1:
            nucleotide = seq[0]
            if nucleotide in nucleotide_indices:
                pwm[i][nucleotide_indices[nucleotide]] = 1
        else:
            norm_factor = 1/len(seq)
            for nucleotide in seq:
                if nucleotide in nucleotide_indices:
                    pwm[i][nucleotide_indices[nucleotide]] = norm_factor

        i+=1

    return pwm



def generate_motifs_for_rbps_ppi(dir):

    # Get 160 rbps
    rbps_file = open(f'{dir}/rbps/rbps.txt', 'r')
    rbps = list(set([rbp.strip() for rbp in rbps_file.read().split('\n')]))

    features1 = pd.read_csv('0_Files/Post-processing/features_motifs_all_132.csv', delimiter='\t')
    features2 = pd.read_csv('0_Files/Post-processing/features_motifs_all_47.csv', delimiter='\t')
    features = pd.concat([features1,features2], axis=1)
    features = features.loc[:,~features.columns.duplicated()].copy() # drop duplicate columns
    features.fillna(0,inplace=True)

    hms = ['H3K27ac', 'H3K27me3','H3K36me3', 'H3K9me3', 'H3K4me3']
    for hm in hms:

        print('\n\n',hm)

        op_dir = f'{dir}/ppi/STRING/{hm}'
        Path(op_dir).mkdir(parents=True, exist_ok=True)

        # Get feature scores for current HM
        features_hm = features[features['type'] == hm]


        # Get HM-assocated proteins
        readers_file = open(f'{dir}/readers/{hm}_readers.txt', 'r')
        readers = list(set([reader.strip() for reader in readers_file.read().split('\n')]))


        # Read rbps in PPI network with chromatin readers of current HM
        ppi = pd.read_csv(f'{dir}/ppi/STRING/{hm}_filtered_PPI.tsv', delimiter='\t')


        # Remove PPIs between RBPs, between Readers 
        ppi= ppi[(ppi['#node1'].isin(rbps) & ~(ppi['node2'].isin(rbps)))].drop_duplicates()
        ppi= ppi[(~ppi['#node1'].isin(readers) & (ppi['node2'].isin(readers)))].drop_duplicates()


        # Get episplicing RBPs for current HM
        rbps_file = open(f"0_Files/Post-processing/epiRBPS/epiRBPS_{hm}.txt", "r")
        epi_rbps = [rbp for rbp in rbps_file.read().split('\n') if rbp]
        
        # drop PPIs involving epiRBPs
        ppi =  ppi[~(ppi['#node1'].isin(epi_rbps) & ~(ppi['node2'].isin(epi_rbps)))].drop_duplicates()
        
        # Get other RBPs
        ppi_rbps = list(set(ppi['#node1'].values.tolist() + ppi['node2'].values.tolist()))
        ppi_rbps = [rbp for rbp in ppi_rbps if rbp in rbps] # exclude reader proteins


        for sf in ppi_rbps:
            print(f'\n\n{sf}')          
            
            sf_motifs =  [seq for seq in list(set(features_hm[sf].values)) if seq != '0.0']
            
            # Generate PSSM
            for sequence in sf_motifs:
                expanded_sequence = expand_iupac(sequence)
                psssm = generate_pwm(expanded_sequence) ## Generate PSSM from expanded sequences
                pssm_df = pd.DataFrame(psssm, columns=['A', 'C', 'G', 'U'])

                # Step 6: Generate sequence logo
                plt.figure(figsize=(len(pssm_df), 2))  # Adjust width to match sequence length
                logo = logomaker.Logo(pssm_df)

                # Style adjustments for aesthetics
                logo.style_spines(visible=False)
                logo.style_xticks(visible=False)
                logo.ax.set_axis_off()

                # Save sequence logo
                plt.savefig(f'{op_dir}/{sf}_{sequence.upper()}.png', dpi=100, bbox_inches='tight', pad_inches=0)  # Adjust resolution and cropping
                plt.close()


if __name__ == "__main__":

    # conveert_to_uniprot(sys.argv[1])
    # read_PPI_string(dir=sys.argv[1])
    # generate_motifs_for_rbps_ppi(dir=sys.argv[1])
    read_PPI_intact(sys.argv[1])
