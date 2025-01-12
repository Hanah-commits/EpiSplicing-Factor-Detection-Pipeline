import os
import json
import sys
from pathlib import Path
from argparse import ArgumentParser

p = ArgumentParser()
p.add_argument("output_dir")
p.add_argument("--process", "-p",
help="The name of the process")
args = p.parse_args()
proc = args.process

with open('paths.json') as f:
        data = json.load(f)
d = data[proc]

bam_files = d['RNASeq files']
ref = d['Reference GTF']
tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
dexseq_scripts = d['DEXSEQ directory']

# create output directory
odir = args.output_dir + 'DEXSEQ/'
Path(odir).mkdir(parents=True, exist_ok=True)

# prepare flattened annotation file
os.system(f"python {dexseq_scripts}python_scripts/dexseq_prepare_annotation_subread.py {ref} -f {odir}DEXSEQ_reference.gtf {ref} {odir}DEXSEQ_reference.gff")

# count exons
os.system(f'bash HelperFunctions/count_all.sh {tissue1} {tissue2} {bam_files}')

# rundexseq
os.system(f'Rscript 1.3_DEXSEQ/run_dexseq.R {bam_files} {tissue1} {tissue2}')