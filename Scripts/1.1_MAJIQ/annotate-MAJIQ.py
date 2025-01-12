import os
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("--process", "-p",
    help="The name of the process")
args = p.parse_args()
proc = args.process

tmp_out_dir = proc + '_0_Files'

flanks = [50,100,200]

print("STEP 4: MAJIQ annotation is in progress")

for flank in flanks:

    os.system(f'bedtools intersect -loj -s -a {tmp_out_dir}/flanks {flank}.bed -b {tmp_out_dir}/MAJIQ/majiq.bed | sort | uniq > {tmp_out_dir}/MAJIQ/majiq_flanks{flank}.bed')