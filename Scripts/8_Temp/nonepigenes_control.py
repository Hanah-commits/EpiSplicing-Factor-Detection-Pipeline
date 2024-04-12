import pandas as pd
import os
import sys
import json
from get_epi_nonepi_flanks import get_nonepigenes

# STEP 1: Process nonepigenes

## STEP 1a: Get common nonepigene set
os.system('python 8_Temp/get_epi_nonepi_flanks.py RMATS epi')

# pre-rbp
os.system('mkdir ../RBPmap')
os.system('python 8_Temp/pre_rbp.py 1')

## STEP 1c: Run RBPmap
os.system('python 7_Post/run_rbpmap.py')

## STEP 1d: Process results
os.system('python 8_Temp/post_rbp.py 1')

## STEP 1e: Move dirs

os.system('mv ../RBPmap ./Post-processing/')