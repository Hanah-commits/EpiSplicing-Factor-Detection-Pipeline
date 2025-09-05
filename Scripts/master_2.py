import os
import sys
from HelperFunctions.check_args import check_args_post_processing


def master_function(part, log_file_name):

    match part:
        case "pool":
            print(f'Starting post-processing..')
            check_args_post_processing()
            sys.stdout = open(log_file_name, 'a')

            # STEP 1: Pool epispliced, non-epispliced and epi non-spliced flanks
            try:
                print('\n\n Pooling Epispliced and Non-epispliced exon flanks \n\n', flush=True)
                os.system(f"python 3_Episplicing/get_epi_nonepi_flanks.py >> {log_file_name} 2>&1")
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

            print('\n\n NEXT: \n\
            1. cp ../RBPmap ../RBPmap_47 \n\
            2. Run RBPmap on the webserver: \n\
                - Query sequences: ../RBPmap_47/*.csv \n\
                - Query motifs: ./HelperFunctions/pssms.txt \n\
                - Ouptut: All_predictions.txt/All_predictions.csv \n\
            3. Save output files locally: \n\
                Eg: ../RBPmap_47/results<rbp_input_epi1.csv>/All_predictions.txt \n\n')

        case "rbpmap":

            # STEP 3: Execute RBPmap
            try:
                print('\n\n Executing RBPmap \n\n', flush=True)
                os.system(f"python 4_RBPMap/run_rbpmap.py   >> {log_file_name} 2>&1")
                print('\n\n----------- DONE -----------', flush=True)
            except Exception as ex:
                print(ex)
                sys.exit(1)

            print('\n\n NEXT: \n\
            1. mv ../RBPmap ../RBPmap_132 \n\n')

        case "classify":

            # STEP 4: Process RBPMap output & Prep Feature Matrix
            ##### NOTE: Manuscript- 132 RBPs from internal db, 47 RBPs (user input pssms) using webserver #####
            ## name output directories accordingly
            ## rbps list in HelperFunctions/
            try:
                print('\n\n Analysing RBP-binding Predictions \n\n', flush=True)
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

            print('\n\n NEXT: \n\
            For each <histone-mark> \n\
            1. View SHAP plot "./0_Files/Post-processing/Analyses/SHAP/<histone-mark>/<histone-mark>.png" \n\
            2. Manually select RBPs that distinguish the two exon classes: \n\
            3. Save RBPs that contribute to the prediction of the positive class as a \\n-separated list in "./0_Files/Post-processing/epiRBPs/SHAP_epiRBPs/rbps_<histone-mark>.txt" \n\
            3. Save RBPs that contribute to the prediction of the negative class as a \\n-separated list in "./0_Files/Post-processing/nonepiRBPs/SHAP_nonepiRBPs/rbps_<histone-mark>.txt" \n\n')

        case "visualize":

            # STEP 6: Visualization: SHAP plots of important features (epi/nonepisplicing RBPs)
             try:
                print('\n\n Collecting Important Features \n\n', flush=True)
                os.system(f"python 5_Classification/classifier.py  >> {log_file_name} 2>&1")
                print('\n\n----------- DONE -----------', flush=True)
            except Exception as ex:
                print(ex)
                sys.exit(1)

            # STEP 7: Visualization: Additional plots, sequence logos
            try:
                print('\n\n Creating Plots \n\n', flush=True)
                os.system(f"python 6_Visualization/visualize_results.py >> {log_file_name} 2>&1")
                os.system(f"python 6_Visualization/visualize_logo.py >> {log_file_name} 2>&1")
                print('\n\n----------- DONE -----------', flush=True)
            except Exception as ex:
                print(ex)
                sys.exit(1)


if __name__ == "__main__":
    # combine epi and non-epiflanks from individual analysis for downstream analyses
    arg = sys.argv[1]
    log_file_name = 'analysis_output.log'
    master_function(arg, log_file_name)
