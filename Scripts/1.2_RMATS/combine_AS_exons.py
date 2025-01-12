import pandas as pd
import os
import json
import sys
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("--process", "-p",
    help="The name of the process")
args = p.parse_args()
proc = args.process

tmp_out_dir = proc + '_0_Files'

with open('paths.json') as f:
    data = json.load(f)
d = data[proc]

fasta = d['Reference fasta']
ref_genome= fasta+".fai"

e = 0
try:
    SE = pd.read_csv(f'{tmp_out_dir}/RMATS/SE_exons.tsv', delimiter='\t', names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI" ], skiprows=1)
except:
    SE = pd.DataFrame()
    e +=1

try:
    MXE = pd.read_csv(f'{tmp_out_dir}/RMATS/MXE_exons.tsv', delimiter='\t',  names=['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI"], skiprows=1)
except:
    MXE = pd.DataFrame()
    e +=1

if e == 2:
    print('No RMATS exons available \n')
    sys.exit()


def A3SS_A5SS_filter(group, subset_column):
    group.sort_values(by=['exonStart_0base', subset_column], ascending=[True, False], inplace=True)
    # Drop duplicates based on 'exonStart_0base' and keep the one with the largest 'dPSI'
    group = group.groupby('exonStart_0base').apply(lambda x: x.loc[x[subset_column].idxmax()])
    

    group.sort_values(by=['exonEnd', subset_column], ascending=[True, False], inplace=True)
    # Drop duplicates based on 'exonEnd' and keep the one with the largest 'dPSI'
    group = group.groupby('exonEnd').apply(lambda x: x.loc[x[subset_column].idxmax()])
    
    return group


## STEP 1: get combined dPSI scores of AS exons

if len(SE) > 0 and len(MXE) >0 : 
    SE_MXE_exons =  pd.concat([SE, MXE], ignore_index=True)
elif len(SE) > 0:
    SE_MXE_exons = SE
else:
    SE_MXE_exons = MXE

## FILTER 1: A3SS/A5SS -  Avoid many A3SS/A5SS versions of AS exons
SE_MXE_exons = SE_MXE_exons.groupby('geneSymbol').apply(lambda x: A3SS_A5SS_filter(x, 'dPSI'))
SE_MXE_exons = SE_MXE_exons.reset_index(drop=True)

## FILTER 2: Keep only AS exons
SE_MXE_exons = SE_MXE_exons[SE_MXE_exons.dPSI > 0.2]

print('# genes with SE and/or MXE exons:   ', len(set(SE_MXE_exons.geneSymbol.values.tolist()))) # log

print(SE_MXE_exons[SE_MXE_exons.geneSymbol == 'ENSG00000105219.8'])

## STEP 2: Get exon flanks

SE_MXE_exons[['chr', "exonStart_0base", "exonEnd", "feature", "score", "strand", "geneSymbol", "dPSI"]].to_csv(f'{tmp_out_dir}/RMATS/rmats_exons_coords.bed', index=False, sep='\t', header=False)

## FILTER 3: Drop flanked AS exns overlapping with TSS regions. CS exons are TSS-free since exon_coords.bed alreeady has TSS-filtered exons

 # extend exon body by 200bp
os.system(f"bedtools slop -i {tmp_out_dir}/RMATS/rmats_exons_coords.bed -g {ref_genome} -b 200 > {tmp_out_dir}/RMATS/rmats_exons_flanked.bed" )

# #Drop flanked exons overlapping with TSS
os.system(f"bedtools intersect -wa -a {tmp_out_dir}/RMATS/rmats_exons_flanked.bed -b {tmp_out_dir}/TSS.bed -s -v > {tmp_out_dir}/RMATS/rmats_exons_filtered.bed")

# # get exon coordinates (remove flanking regions)
os.system(f"bedtools slop -i {tmp_out_dir}/RMATS/rmats_exons_filtered.bed -g {ref_genome} -l -200 -r -200 -s > {tmp_out_dir}/RMATS/rmats_exons_coords.bed")

 # exon boundary external flanks
os.system(f"bedtools flank -i {tmp_out_dir}/RMATS/rmats_exons_coords.bed -g {ref_genome}  -b 200 > {tmp_out_dir}/flanks.bed" )

# separate start,stop flank coords
os.system(f"sed -n 'n;p' {tmp_out_dir}/flanks.bed > {tmp_out_dir}/stop.bed")
os.system(f"sed -n 'p;n' {tmp_out_dir}/flanks.bed > {tmp_out_dir}/start.bed")

# exon boundary internal flanks
os.system(f"bedtools slop -i {tmp_out_dir}/start.bed -g {ref_genome} -l 0 -r 200 > {tmp_out_dir}/start_flanks.bed")
os.system(f"bedtools slop -i {tmp_out_dir}/stop.bed -g {ref_genome} -l 200 -r 0 > {tmp_out_dir}/stop_flanks.bed")

# combine start,stop flank coords
os.system(f"paste -d'\n' {tmp_out_dir}/start_flanks.bed {tmp_out_dir}/stop_flanks.bed | sort -k1,1 -k2,2n > {tmp_out_dir}/RMATS/rmats_flanks200.bed")

# remove intermediate files
os.system(f"rm {tmp_out_dir}/start*.bed")
os.system(f"rm  {tmp_out_dir}/stop*.bed")
os.system(f"rm {tmp_out_dir}/flanks.bed")
os.system(f"rm {tmp_out_dir}/RMATS/*flanked*.bed")
os.system(f"rm {tmp_out_dir}/RMATS/rmats_exons_filtered.bed")

# #%# Note: Exons have varying lengths. Flanks can overlap.