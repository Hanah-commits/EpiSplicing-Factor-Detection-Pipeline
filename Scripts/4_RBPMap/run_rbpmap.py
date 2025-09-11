import os
import json


# STEP 1: Run RBPmap

with open('paths.json') as f:
            d = json.load(f)

rbp_path = list(set([sub_dict["RBPmap directory"] for sub_dict in d.values() if "RBPmap directory" in sub_dict]))[0]


currdir = os.getcwd()
file = os.getcwd() + '/../RBPmap/input.txt'
os.chdir(rbp_path)

os.system(f"parallel --joblog parallel_log -a {file} perl RBPmap_EpiSplicing.pl -input {{}} -genome 'human' -db 'hg38' -db_motifs all_human -stringency high -conservation on")

os.chdir(currdir)