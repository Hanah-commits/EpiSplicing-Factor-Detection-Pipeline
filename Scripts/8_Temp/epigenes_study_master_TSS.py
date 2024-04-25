import pandas as pd
import os
import sys

tool = sys.argv[1]
mode = sys.argv[2]

# Rscript 8_Temp/geneid_to_genesymbol_TSS.R
# python 8_Temp/get_TSS_flanks.py


# get flanks of TSS exons of epigenes for current analysis
os.system(f'python 8_Temp/get_epi_nonepi_flanks_TSS.py {tool} {mode}')

# pre-rbp
os.system('mkdir ../RBPmap')
os.system('python 8_Temp/pre_rbp.py 0')

# run rbpmap
os.system('python 7_Post/run_rbpmap.py')

# post-rbpmap
os.system('python 8_Temp/post_rbp.py 0')

# get nonepi flanks
os.system(f'cp ./Post-processing/nonepi_flanks.bed ./0_Files/Post-processing/')
os.system(f'cp ./Post-processing/FilteredZscores_nonepi.csv ./0_Files/Post-processing/')

# features
os.system('python 7_Post/features.py flanks')
os.system('python 8_Temp/classifier_features.py')

# # classifier
os.system('python 8_Temp/classifier.py ./')

# #DEA
os.system('python 8_Temp/enrichment_hm.py')
