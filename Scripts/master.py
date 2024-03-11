import os
import sys
from HelperFunctions.check_args import check_args, move_dirs


def master_function():
    # Check input arguments from paths.json
    output_dir = check_args()
    print(f'Starting execution of pipeline')
    log_file_name = 'output.log'
    sys.stdout = open(log_file_name, 'w')

    #check command line arguments
    weights = False
    if len(sys.argv) > 1 and sys.argv[1] == "-w":
        weights = True

    # STEP 0: Preprocessing

    # Differential Expression Analysis
    if weights:
        try:
            print('\n\n Peforming DGE Analysis \n\n')
            os.system(f"python PreProcessing/featureCounts.py   >> {log_file_name} 2>&1")
            os.system(f"Rscript PreProcessing/Limma.R   >> {log_file_name}")
            print('\n\n----------- DONE -----------')
        except Exception as ex:
            print(ex)
            move_dirs(output_dir)
            sys.exit(1)

    # Prepare flank reference : 50, 100, 200 bp
    try:
        print('\n\n Extracting Candidate Exons \n\n')
        os.system(f"python PreProcessing/prepare_FlanksRef.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## STEP 1: Differential Exon Usage

    ## 1.1 Execute MAJIQ 
    try:
        print('\n\n Executing MAJIQ \n\n')
        os.system(f"python 1.1_MAJIQ/runMAJIQ.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## 1.2 Execute RMATS
    try:
        print('\n\n Executing RMATS \n\n')
        os.system(f"python 1.2_RMATS/runRMATS.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)    

    ## 1.3 Execute DXESEQ
    try:
        print('\n\n Executing DEXSEQ \n\n')
        os.system(f"python 1.3_DEXSEQ/runDEXSEQ.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 2: Execute MANorm -  Differential Histone Modifications
    try:
        print('\n\n Executing MANorm \n\n')
        os.system(f"python 2_MANorm/manorm_all.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 3: Annotate all candidate exon flanks with MANorm peaks
        
    ## 3.1 Annotate exon flanks with peaks of available HMs    
    try:
        print('\n\n MANorm: Annotating Exons \n\n')
        os.system(f'python 2_MANorm/annotate_MANorm_all_exons.py {output_dir}   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## 3.2 Combine annotations into single file
    try:
        print('\n\n MANorm: Combining HM-exon Annotations \n\n')
        os.system(f'python 2_MANorm/combine_all_HMpeaks.py   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## 3.3 Convert ENSEMBL geneid to gene name
    try:
        print('\n\n MANorm: Converting from ENSEMBL GeneID to GeneSymbol \n\n')
        os.system(f'Rscript PreProcessing/gene_id_to_gene_symbol.R   >> {log_file_name}')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## STEP 4: Process MAJIQ output and annotate its exon flanks with peaks

    ## 4.1 Process majiq output
    try:
        print('\n\n MAJIQ: Analysing MAJIQ Output \n\n')
        os.system(f"python 1.1_MAJIQ/post-MAJIQ.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## 4.2 BEDTools - Annotate exon flanks with MAJIQ junctions
    try:
        print('\n\n MAJIQ: Annotating Candidate Flanks with DEU Score\n\n')
        os.system(f"python 1.1_MAJIQ/annotate-MAJIQ.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## 4.3 Process BEDTools output
    try:
        print('\n\n MAJIQ: Analysing BEDtools Output \n\n')
        os.system(f"python 1.1_MAJIQ/post-bedtools.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    ## 4.4 Process annotated flanks (TSS-filtering)
    try:
        print('\n\n MAJIQ: Including Constiutively Spliced Exon-flanks \n\n')
        os.system(f"python 1.1_MAJIQ/combine_AS_CS_flanks.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)   

    ## 4.5 Annotate TSS-filtered flanks with peaks
    try:
        print('\n\n MAJIQ: Annotating Candidate Flanks with DHM Score \n\n')
        os.system(f'python 2_MANorm/DHM_flanks_MAJIQ.py   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)


    # STEP 5: Process DEXSEQ output and annotate its exon flanks with peaks

    ## 5.1 Get deu exons, TSS-filter flanks of these exons    
    try:
        print('\n\n DEXSEQ: Analysing DEXSEQ Output & Annotating Exons with DEU Score \n\n')
        os.system(f'python 1.3_DEXSEQ/post_dexseq.py {output_dir}   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)   

    ## 5.2 Annotate TSS-filtered flanks with peaks
    try:
        print('\n\n DEXSEQ: Annotating Candidate Exons with DHM Score \n\n')
        os.system(f'python 2_MANorm/DHM_flanks_DEXSEQ.py   >> {log_file_name} 2>&1' )
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1) 

    ## STEP 6: Process RMATS output and annotate its exon flanks with peaks

    ## 6.1 Extract skipped exons   
    try:
        print('\n\n RMATS: Analysing RMATS Output - Skipped Exons \n\n')
        os.system(f'python 1.2_RMATS/get_SE.py {output_dir}   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1) 

    ## 6.2 Extract mutually exclusive exons
    try:
        print('\n\n RMATS: Analysing RMATS Output - Mutually Exclusive Exons \n\n')
        os.system(f'python 1.2_RMATS/get_MXE.py {output_dir}   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1) 

    ## 6.3 Get flanks of alternative exons
    try:
        print('\n\n RMATS: Combining DEU Scores of Skipped and Mutally Exclusive Exons \n\n')
        os.system(f'python 1.2_RMATS/combine_AS_exons.py   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1) 


    ## 6.4 annotate TSS-filtered flanks with peaks
    try:
        print('\n\n RMATS: Annotating Candidate Exons with DHM Score \n\n')
        os.system(f'python 2_MANorm/DHM_flanks_RMATS.py   >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1) 

    ## STEP 7: DEU - DHM Correlation
    try:
        print('\n\n Computing DEU-DHM Correlation \n\n')
        os.system(f"python 3_Episplicing/correlation_plot.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 9: Prepare RBPmap input
    try:
        print('\n\n Preparing Input for RBPmap \n\n')
        os.system(f"python 4_RBPMap/pre-rbp.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 10: Execute RBPmap
    try:
        print('\n\n Executing RBPmap \n\n')
        os.system(f"python 4_RBPMap/run_rbpmap.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 11: Process RBPMap output
    try:
        print('\n\n RBPmap: Analysing RBP-binding Predictions \n\n')
        os.system(f"python 4_RBPMap/post-rbp.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    try:    
        print('\n\n Preparing Feature Matrix \n\n')
        os.system(f"python 5_Classification/rbp_pvals.py   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 12: Add logFC weights to binding scores from RBPMap
    if weights:
        try:
            print('\n\n Adding Feature Weights \n\n')
            os.system(f"python PreProcessing/get_weights.py   >> {log_file_name} 2>&1")
            os.system(f"python 4_RBPMap/rbp-weights.py   >> {log_file_name} 2>&1")
            print('\n\n----------- DONE -----------')
        except Exception as ex:
            print(ex)
            move_dirs(output_dir)
            sys.exit(1)

    # STEP 13: Prep Feature Matrix
    try:
        print('\n\n Binary Classifier: Preparing Input \n\n')
        os.system(f"python 5_Classification/features.py {output_dir} {str(weights)}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    try:
        print('\n\n Binary Classifier: Preparing Input \n\n') 
        os.system(f"python 5_Classification/classifier_features.py {output_dir} {str(weights)}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 14: Binary Classification
    try:
        print('\n\n Binary Classifier: Performing Classification \n\n')
        os.system(f"python 5_Classification/classifier.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 15: Enrichment
    try:
        print('\n\n Performing Differential RBP-Enrichment Analysis \n\n')
        os.system(f"python 6_Enrichment/enrichment.py {output_dir}   >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------')
    except Exception as ex:
        print(ex)
        move_dirs(output_dir)
        sys.exit(1)

    # STEP 15: Move files generated from current pipeline run to
    move_dirs(output_dir)
    print(f'Execution of   pipeline finished successfully!')

    sys.stdout.close()


if __name__ == "__main__":
    master_function()