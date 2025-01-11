import pandas as pd
import os
import sys


def pool_flanks(tool):

    # STEP 1: Get control flanks (nonepispliced exon flanks)
    os.system(f'python 8_Temp/get_epi_nonepi_flanks.py {tool} nonepi')

    # STEP 2: Get epi flanks
    os.system(f'python 8_Temp/get_epi_nonepi_flanks.py {tool} epi')

def preprocess_run_rbpmap():

    os.system('mkdir ../RBPmap')
    os.system('python 8_Temp/pre_rbp.py')
    os.system('python 8_Temp/run_rbpmap.py') # 132 internal RBPs


if __name__ == "__main__":

    tool = sys.argv[1]

    # STEP 1: Get control flanks (nonepispliced exon flanks)
    pool_flanks(tool)

    # STEP 2: Run RBPmap
    preprocess_run_rbpmap()

    # STEP 3: Get RBPmap predictions
    ## 132 RBPs from internal db, 47 RBPs (user input) using webserver 
    ## rename output directories accordingly
    os.system('python 8_Temp/post_rbp.py')