import os
import sys
from HelperFunctions.check_args import check_args_post_processing


def master_function():

    log_file_name = 'analysis_output.log'
    sys.stdout = open(log_file_name, 'a')

    # STEP 1: Pool epi and nonepiflanks
    try:
        print('\n\n Pooling Epispliced and Non-epispliced exon flanks \n\n', flush=True)
        os.system(f"python 3_Episplicing/get_epi_nonepi_flanks.py nonepi >> {log_file_name} 2>&1")
        os.system(f"python 3_Episplicing/get_epi_nonepi_flanks.py epi >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1)


    # STEP 2: Prepare RBPmap input
    try:
        print('\n\n Preparing Input for RBPmap \n\n', flush=True)
        os.system('mkdir ../RBPmap')
        os.system(f"python 4_RBPMap/pre-rbp.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1)

    # STEP 3: Execute RBPmap
    try:
        print('\n\n Executing RBPmap \n\n', flush=True)
        os.system(f"python 4_RBPMap/run_rbpmap.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1)

    # STEP 4: Process RBPMap output & Prep Feature Matrix
    ##### NOTE: Manuscript- 132 RBPs from internal db, 47 RBPs (user input pssms) using webserver #####
    ## name output directories accordingly
    ## rbps list in HelperFunctions/
    try:
        print('\n\n RBPmap: Analysing RBP-binding Predictions \n\n', flush=True)
        os.system(f"python 4_RBPMap/post-rbp.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1)

    # STEP 5: Binary Classification
    try:
        print('\n\n Binary Classifier: Performing Classification \n\n', flush=True)
        os.system(f"python 5_Classification/classifier.py  >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1)

    # STEP 6: Visualization: Additional plots, sequence logos
    try:
        print('\n\n Performing Differential RBP-Enrichment Analysis \n\n', flush=True)
        os.system(f"python 6_Visualization/visualize.py >> {log_file_name} 2>&1")
        os.system(f"python 6_Visualization/visualize_logo.py >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1)

    # STEP 7: eCLIP Annotation
    try:
        print('\n\n Performing Differential RBP-Enrichment Analysis \n\n', flush=True)
        os.system(f"python 7_Post/eclip_annotate.py >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        sys.exit(1) 



if __name__ == "__main__":
    # combine epi and non-epiflanks from individual analysis for downstream analyses
    print(f'Starting post-processing..')
    check_args_post_processing()
    master_function()
