import os
import json
import glob
from pathlib import Path
from argparse import ArgumentParser

# Get the process name, use it in the output directory
def get_argument_parser():
    p = ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("--process", "-p",
        help="The name of the process")
    return p

if __name__ == "__main__":

    
    p = get_argument_parser()
    args = p.parse_args()

    proc = args.process
    
    with open('paths.json') as f:
        data = json.load(f)
    d = data[proc]

    rna_files_dir = d['RNASeq files']
    ref = d['Reference GTF']
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]
    rmats_dir = d['RMATS directory']
    read_len = d['read_length']


    currdir = os.getcwd()
    
    output_dir = args.output_dir + 'RMATS/'
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # prepare input
    for tissue in [tissue1, tissue2]:
        
        # find input bam files
        tissue_files = glob.glob(rna_files_dir + f'{tissue}*.bam')

        # write to input file
        with open(f'{output_dir}/{tissue}.txt', 'w') as file:
            file.write(','.join(map(str, tissue_files)))

    # run rmats
    ## -t single (to include paired reads and those without a proper pair), --libType fr-firststrand (reverse-stranded), --allow-clipping
    os.system(f'python {rmats_dir}rmats.py --b1 {output_dir}{tissue1}.txt --b2 {output_dir}{tissue2}.txt --gtf {ref} -t single --variable-read-length --libType fr-firststrand --allow-clipping --nthread 8 --od {output_dir} --tmp {output_dir}tmp/ --task both')