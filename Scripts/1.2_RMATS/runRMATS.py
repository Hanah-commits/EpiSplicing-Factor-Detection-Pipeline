import os
import sys
import json
import glob
from pathlib import Path

if __name__ == "__main__":

    
    with open('paths.json') as f:
        d = json.load(f)

    rna_files_dir = d['RNASeq files']
    ref = d['Reference GTF']
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]
    rmats_dir = d['RMATS directory']
    read_len = d['read_length']


    currdir = os.getcwd()
    
    output_dir = sys.argv[1] + 'RMATS/'
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # prepare input
    for tissue in [tissue1, tissue2]:
        
        # find input bam files
        tissue_files = glob.glob(rna_files_dir + f'{tissue}*.bam')

        # write to input file
        with open(f'{output_dir}/{tissue}.txt', 'a') as file:
            file.write(','.join(map(str, tissue_files)))

    # run rmats
    os.system(f'python {rmats_dir}rmats.py --b1 {output_dir}{tissue1}.txt --b2 {output_dir}{tissue2}.txt --gtf {ref} -t paired --readLength {read_len} --variable-read-length --libType fr-unstranded --allow-clipping --novelSS --nthread 4 --od {output_dir} --tmp {output_dir}tmp/ --task both')