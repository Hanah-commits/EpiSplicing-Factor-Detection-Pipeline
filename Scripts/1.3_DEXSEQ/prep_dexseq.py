import os
import json
import sys
from pathlib import Path


with open('paths.json') as f:
    d = json.load(f)

bam_files = d['RNASeq files']
ref = d['Reference genome']
tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
dexseq_scripts = d['DEXSEQ directory']

# create output directory
odir = sys.argv[1] + 'DEXSEQ/'
Path(odir).mkdir(parents=True, exist_ok=True)

# prepare flattened annotation file
os.system(f"python {dexseq_scripts}/python_scripts/dexseq_prepare_annotation.py {ref} {odir}/gencode.v24.DEXSEQ.gff")

# count exons
os.system(f'bash HelperFunction/count_all.sh {bam_files}{tissue1}_*.bam ')
os.system(f'bash HelperFunction/count_all.sh {bam_files}{tissue2}_*.bam ')