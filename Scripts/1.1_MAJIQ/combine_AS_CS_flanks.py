import pandas as pd
import os


AS_flanks = pd.read_csv(f'0_Files/Filtered_dPSI.csv', delimiter='\t')
all_flanks = pd.read_csv('0_Files/flanks200.bed', delimiter='\t', names=['seqid', "start", "stop", "feature", "score", "strand", "gene_id"], skiprows=1)

## STEP 1: fetch CS exons
# add label
AS_flanks['feature'] = 'AS'
AS_flanks['score'] = '.'
all_flanks['feature'] = 'CS'
all_flanks['mean_dpsi_per_lsv_junction'] = 0.0

# list of SE genes
AS_genes = list(set(AS_flanks.gene_id.values.tolist()))

combined_AS_CS = []

for gene in AS_genes:
    AS_flanks_gene = AS_flanks[AS_flanks.gene_id == gene]
    flanks_gene = all_flanks[all_flanks.gene_id == gene]
    flanks_gene= pd.concat([AS_flanks_gene, flanks_gene], ignore_index=True)
    flanks_gene = flanks_gene.drop_duplicates(subset=['seqid', 'start', 'stop', 'gene_id'])
    combined_AS_CS.append(flanks_gene)

# Concatenate the dataframes in the list
all_genes = pd.concat(combined_AS_CS, ignore_index=True)

#@ STEP 2: Get exon flanks

all_genes[['seqid', "start", "stop", "feature", "score", "strand", "gene_id", "mean_dpsi_per_lsv_junction"]].to_csv('0_Files/majiq_filtered_flanks.bed', index=False, sep='\t', header=False)

## FILTER 1: Drop flanked AS exns overlapping with TSS regions. CS exons are TSS-free since exon_coords.bed alreeady has TSS-filtered exons

os.system('bedtools intersect -wa -a 0_Files/majiq_filtered_flanks.bed -b 0_Files/TSS.bed -s -v > 0_Files/majiq_filtered_flanks_temp.bed && mv 0_Files/majiq_filtered_flanks_temp.bed 0_Files/majiq_filtered_flanks.bed')

#%# Note: Exons have varying lengths. Flanks can overlap.