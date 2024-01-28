import pandas as pd
import os
import json


with open('paths.json') as f:
    d = json.load(f)

ref = d['Reference genome']
fasta = d['Reference fasta']
ref_genome= fasta+".fai"

SE_AS = pd.read_csv(f'0_Files/RMATS/SE_exons_AS.tsv', delimiter='\t', names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI" ], skiprows=1)
MXE_AS = pd.read_csv(f'0_Files/RMATS/MXE_exons_AS.tsv', delimiter='\t',  names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI"], skiprows=1)
all_exons = pd.read_csv('0_Files/exon_coords.bed', delimiter='\t', names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol"], skiprows=1)


def A3SS_A5SS_filter(group, subset_column):
    group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
    # Drop duplicates based on 'exonStart_0base' and keep the one with the largest 'dPSI'
    group = group.groupby('exonStart_0base').apply(lambda x: x.loc[x[subset_column].idxmax()])
    

    group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
    # Drop duplicates based on 'exonEnd' and keep the one with the largest 'dPSI'
    group = group.groupby('exonEnd').apply(lambda x: x.loc[x[subset_column].idxmax()])
    
    return group


## STEP 1: get combined dPSI scores of AS exons
SE_MXE_exons =  pd.concat([SE_AS, MXE_AS], ignore_index=True)

## FILTER 1: A3SS/A5SS -  Avoid many A3SS/A5SS versions of AS exons
SE_MXE_exons = SE_MXE_exons.groupby('geneSymbol').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
SE_MXE_exons = SE_MXE_exons.reset_index(drop=True)

## STEP 2: fetch CS exons
# add label
SE_MXE_exons['label'] = 'AS'
all_exons['label'] = 'CS'
all_exons['dPSI'] = 0.0

# list of SE genes
AS_genes = list(set(SE_MXE_exons.geneSymbol.values.tolist()))

combined_AS_CS = []

for gene in AS_genes:
    AS_exons = SE_MXE_exons[SE_MXE_exons.geneSymbol == gene]
    exons_gene = all_exons[all_exons.geneSymbol == gene]
    exons_gene= pd.concat([AS_exons, exons_gene], ignore_index=True)
    exons_gene = exons_gene.drop_duplicates(subset=['chr', 'exonStart_0base', 'exonEnd', 'geneSymbol'])
    combined_AS_CS.append(exons_gene)

# Concatenate the dataframes in the list
all_genes = pd.concat(combined_AS_CS, ignore_index=True)

## FILTER 2: A3SS/A5SS -  Avoid marking A3SS/A5SS versions of AS exons as CS
all_genes = all_genes.groupby('geneSymbol').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
all_genes = all_genes.reset_index(drop=True)

#@ STEP 3: Get exon flanks

all_genes[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI"]].to_csv('0_Files/RMATS/rmats_exons_coords.bed', index=False, sep='\t', header=False)

 # exon boundary external flanks
os.system("bedtools flank -i 0_Files/RMATS/rmats_exons_coords.bed -g " + ref_genome + " -b 200 > 0_Files/flanks.bed" )

# separate start,stop flank coords
os.system("sed -n 'n;p' 0_Files/flanks.bed > 0_Files/stop.bed")
os.system("sed -n 'p;n' 0_Files/flanks.bed > 0_Files/start.bed")

# exon boundary internal flanks
os.system("bedtools slop -i 0_Files/start.bed -g " + ref_genome + " -l 0 -r 200 > 0_Files/start_flanks.bed")
os.system("bedtools slop -i 0_Files/stop.bed -g " + ref_genome +" -l 200 -r 0 > 0_Files/stop_flanks.bed")

# combine start,stop flank coords
os.system("paste -d'\n' 0_Files/start_flanks.bed 0_Files/stop_flanks.bed | sort -k1,1 -k2,2n > 0_Files/RMATS/rmats_flanks200.bed")

# remove intermediate files
os.system("rm 0_Files/start*.bed")
os.system("rm  0_Files/stop*.bed")
os.system("rm 0_Files/flanks.bed")


## FILTER 3: Drop flanked AS exns overlapping with TSS regions. CS exons are TSS-free since exon_coords.bed alreeady has TSS-filtered exons

os.system('bedtools intersect -wa -a 0_Files/RMATS/rmats_flanks200.bed -b 0_Files/TSS.bed -s -v > 0_Files/rmats_flanks200_temp.bed && mv 0_Files/rmats_flanks200_temp.bed 0_Files/RMATS/rmats_flanks200.bed')

#%# Note: Exons have varying lengths. Flanks can overlap.