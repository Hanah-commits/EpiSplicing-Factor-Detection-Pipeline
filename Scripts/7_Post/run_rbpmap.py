import os
import json


# STEP 1: Run RBPmap

with open('paths.json') as f:
            d = json.load(f)

rbp_path = d["RBPmap directory"]

currdir = os.getcwd()
file = os.getcwd() + '/../RBPmap/input.txt'
os.chdir(rbp_path)

os.system("parallel --joblog parallel_log -a " + file + " perl RBPmap_EpiSplicing.pl -input {1} -genome 'human' -db 'hg38' -db_motifs all_human")

os.chdir(currdir)