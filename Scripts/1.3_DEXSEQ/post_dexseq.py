import pandas as pd
import os
import json
from pathlib import Path
import sys


with open('paths.json') as f:
    d = json.load(f)

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
fasta = d['Reference fasta']
ref_genome= fasta+".fai"

# STEP 0: Create directories to store DEXSEQ files
output_dir = str(Path(os.getcwd())) + "/0_Files/DEXSEQ/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

# read dexseq output
file = sys.argv[1] + f'DEXSEQ/DEXSEQ_{tissue1}_{tissue2}.tsv'
log2FC = f"log2fold_{tissue2}_{tissue1}"
dexseq = pd.read_csv(file, delimiter='\t')[['groupID', 'featureID', 'dispersion', 'stat', 'pvalue', 'padj', log2FC, 'genomicData.seqnames', 'genomicData.start', 'genomicData.end', 'genomicData.strand']]

print('Processing DEXSEQ output: \n')

### STEP 1: Get DEU exons

# remove exons with nan
dexseq = dexseq[~dexseq.padj.isna()]

temp = dexseq.assign(groupID=dexseq['groupID'].str.split('+')).explode('groupID')
print('# genes reported:                ', len(set(temp.groupID.values.tolist()))) # log

# get significant exons
dexseq = dexseq[dexseq.padj < 0.05]

temp = dexseq.assign(groupID=dexseq['groupID'].str.split('+')).explode('groupID')
print('FDR-adj pvalue <= 0.05:          ', len(set(temp.groupID.values.tolist()))) # log

## keep exons with | fold change | >= 2
dexseq = dexseq[dexseq[log2FC].abs() >= 1]

temp = dexseq.assign(groupID=dexseq['groupID'].str.split('+')).explode('groupID')
print('| log2FC | > 0.5:                ', len(set(temp.groupID.values.tolist()))) # log

# genes with 2+ exons
dexseq = dexseq.assign(groupID=dexseq['groupID'].str.split('+')).explode('groupID')
dexseq = dexseq[dexseq.groupby('groupID').groupID.transform(len) > 2]

temp = dexseq.assign(groupID=dexseq['groupID'].str.split('+')).explode('groupID')
print('Filtered out genes < 3 exons:    ', len(set(temp.groupID.values.tolist()))) # log

### STEP 2:  Get exons flanks
dexseq['feature'] = 'dexseq_exon'
dexseq['score'] = '.'

dexseq[['genomicData.seqnames', 'genomicData.start', 'genomicData.end', 'feature', 'score', 'genomicData.strand', 'groupID', log2FC]].to_csv('0_Files/DEXSEQ/dexseq_exons_coords.bed', index=False, sep='\t', header=False)

# exon boundary external flanks
os.system("bedtools flank -i 0_Files/DEXSEQ/dexseq_exons_coords.bed -g " + ref_genome + " -b 200 > 0_Files/flanks.bed" )

# separate start,stop flank coords
os.system("sed -n 'n;p' 0_Files/flanks.bed > 0_Files/stop.bed")
os.system("sed -n 'p;n' 0_Files/flanks.bed > 0_Files/start.bed")

# exon boundary internal flanks
os.system("bedtools slop -i 0_Files/start.bed -g " + ref_genome + " -l 0 -r 200 > 0_Files/start_flanks.bed")
os.system("bedtools slop -i 0_Files/stop.bed -g " + ref_genome +" -l 200 -r 0 > 0_Files/stop_flanks.bed")

# combine start,stop flank coords
os.system("paste -d'\n' 0_Files/start_flanks.bed 0_Files/stop_flanks.bed | sort -k1,1 -k2,2n > 0_Files/DEXSEQ/dexseq_flanks200.bed")

# remove intermediate files
os.system("rm 0_Files/start*.bed")
os.system("rm  0_Files/stop*.bed")
os.system("rm 0_Files/flanks.bed")


## FILTER 3: Drop flanked AS exns overlapping with TSS regions. CS exons are TSS-free since exon_coords.bed alreeady has TSS-filtered exons

os.system('bedtools intersect -wa -a 0_Files/DEXSEQ/dexseq_flanks200.bed -b 0_Files/TSS.bed -s -v > 0_Files/dexseq_flanks200_temp.bed && mv 0_Files/dexseq_flanks200_temp.bed 0_Files/DEXSEQ/dexseq_flanks200.bed')
