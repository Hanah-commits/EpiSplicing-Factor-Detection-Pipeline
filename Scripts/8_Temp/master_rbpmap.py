import sys
import os

mode = sys.argv[1]

# STEP 1: Get flanks of AS exons (from epi and nonepi genes)

os.system('python 8_Temp/get_epi_nonepi_flanks.py')

# STEP 2: Prepare RBPmap input

os.system(f'python 7_Post/pre_rbp.py {mode}')

# STEP 7: Run RBPmap

os.system('python 7_Post/run_rbpmap.py')

# STEP 8: Process RBPmap results

os.system(f'python 7_Post/post_rbp.py')

# STEP 9: Prepare feature matrix

os.system('python 7_Post/features.py')
os.system('python 7_Post/classifier_features.py')

# STEP 10: Run classifier
os.system('python 7_Post/classifier.py 0_Files/Post-processing/')

# STEP 11: Differential Enrichment Analysis
os.system('python 7_Post/enrichment.py')
