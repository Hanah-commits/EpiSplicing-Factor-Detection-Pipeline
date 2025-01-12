import pandas as pd
import os
import json
from pathlib import Path
import sys
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("output_dir")
p.add_argument("--process", "-p",
    help="The name of the process")
args = p.parse_args()
proc = args.process

tmp_out_dir = proc + '_0_Files'

with open('paths.json') as f:
    data = json.load(f)
d = data[proc]

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
fasta = d['Reference fasta']
ref_genome= fasta+".fai"

# STEP 0: Create directories to store DEXSEQ files
output_dir = str(Path(os.getcwd())) + f"/{tmp_out_dir}/DEXSEQ/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

# read dexseq output
file = args.output_dir + f'DEXSEQ/DEXSEQ_{tissue1}_{tissue2}.tsv'
try:
    log2FC = f"log2fold_{tissue2}_{tissue1}"
    dexseq = pd.read_csv(file, delimiter='\t')[['groupID', 'featureID', 'dispersion', 'stat', 'pvalue', 'padj', log2FC, 'genomicData.seqnames', 'genomicData.start', 'genomicData.end', 'genomicData.strand']]
except:
    log2FC = f"log2fold_{tissue1}_{tissue2}"
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

dexseq[['genomicData.seqnames', 'genomicData.start', 'genomicData.end', 'feature', 'score', 'genomicData.strand', 'groupID', log2FC]].to_csv(f'{tmp_out_dir}/DEXSEQ/dexseq_exons_coords.bed', index=False, sep='\t', header=False)

# exon boundary external flanks
os.system(f"bedtools flank -i {tmp_out_dir}/DEXSEQ/dexseq_exons_coords.bed -g {ref_genome} -b 200 > {tmp_out_dir}/flanks.bed" )

# separate start,stop flank coords
os.system(f"sed -n 'n;p' {tmp_out_dir}/flanks.bed > {tmp_out_dir}/stop.bed")
os.system(f"sed -n 'p;n' {tmp_out_dir}/flanks.bed > {tmp_out_dir}/start.bed")

# exon boundary internal flanks
os.system(f"bedtools slop -i {tmp_out_dir}/start.bed -g {ref_genome} -l 0 -r 200 > {tmp_out_dir}/start_flanks.bed")
os.system(f"bedtools slop -i {tmp_out_dir}/stop.bed -g {ref_genome} -l 200 -r 0 > {tmp_out_dir}/stop_flanks.bed")

# combine start,stop flank coords
os.system(f"paste -d'\n' {tmp_out_dir}/start_flanks.bed {tmp_out_dir}/stop_flanks.bed | sort -k1,1 -k2,2n > {tmp_out_dir}/DEXSEQ/dexseq_flanks200.bed")

# remove intermediate files
os.system(f"rm {tmp_out_dir}/start*.bed")
os.system(f"rm  {tmp_out_dir}/stop*.bed")
os.system(f"rm {tmp_out_dir}/flanks.bed")


## FILTER 3: Drop flanked AS exns overlapping with TSS regions. CS exons are TSS-free since exon_coords.bed alreeady has TSS-filtered exons

os.system(f'bedtools intersect -wa -a {tmp_out_dir}/DEXSEQ/dexseq_flanks200.bed -b {tmp_out_dir}/TSS.bed -s -v > {tmp_out_dir}/dexseq_flanks200_temp.bed && mv {tmp_out_dir}/dexseq_flanks200_temp.bed {tmp_out_dir}/DEXSEQ/dexseq_flanks200.bed')