import os
import sys
from HelperFunctions.check_args import check_args, move_dirs
from multiprocessing import Process


# Check input arguments from paths.json
list_of_processes, output_dirs = check_args()

# The function that will be executed by all processes 
def master_function(proc, output_dir):
    print(f'Starting execution of {proc} pipeline')
    log_file_name = f'{proc}_output.log'
    sys.stdout = open(log_file_name, 'a')

    # STEP 0: Preprocessing

    # Prepare flank reference : 50, 100, 200 bp
    try:
        print('\n\n Extracting Candidate Exons \n\n', flush=True)
        os.system(f"python PreProcessing/prepare_FlanksRef.py -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## STEP 1: Differential Exon Usage

    ## 1.1 Execute MAJIQ 
    try:
        print('\n\n Executing MAJIQ \n\n', flush=True)
        os.system(f"python 1.1_MAJIQ/runMAJIQ.py {output_dir} -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## 1.2 Execute RMATS
    try:
        print('\n\n Executing RMATS \n\n', flush=True)
        os.system(f"python 1.2_RMATS/runRMATS.py {output_dir} -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)    

    ## 1.3 Execute DXESEQ
    try:
        print('\n\n Executing DEXSEQ \n\n', flush=True)
        os.system(f"python 1.3_DEXSEQ/runDEXSEQ.py {output_dir} -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    # STEP 2: Execute MANorm -  Differential Histone Modifications
    try:
        print('\n\n Executing MANorm \n\n', flush=True)
        os.system(f"python 2_MANorm/manorm_all.py {output_dir} -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    # STEP 3: Annotate all candidate exon flanks with MANorm peaks
        
    ## 3.1 Annotate exon flanks with peaks of available HMs    
    try:
        print('\n\n MANorm: Annotating Exons \n\n', flush=True)
        os.system(f'python 2_MANorm/annotate_MANorm_all_exons.py {output_dir} -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## 3.2 Combine annotations into single file
    try:
        print('\n\n MANorm: Combining HM-exon Annotations \n\n', flush=True)
        os.system(f'python 2_MANorm/combine_all_HMpeaks.py -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## 3.3 Convert ENSEMBL geneid to gene name
    try:
        print('\n\n MANorm: Converting from ENSEMBL GeneID to GeneSymbol \n\n', flush=True)
        os.system(f'Rscript PreProcessing/gene_id_to_gene_symbol.R -p {proc} >> {log_file_name}')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## STEP 4: Process MAJIQ output and annotate its exon flanks with peaks

    ## 4.1 Process majiq output
    try:
        print('\n\n MAJIQ: Analysing MAJIQ Output \n\n', flush=True)
        os.system(f"python 1.1_MAJIQ/post-MAJIQ.py {output_dir} -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## 4.2 BEDTools - Annotate exon flanks with MAJIQ junctions
    try:
        print('\n\n MAJIQ: Annotating Candidate Flanks with DEU Score\n\n', flush=True)
        os.system(f"python 1.1_MAJIQ/annotate-MAJIQ.py -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## 4.3 Process BEDTools output
    try:
        print('\n\n MAJIQ: Analysing BEDtools Output \n\n', flush=True)
        os.system(f"python 1.1_MAJIQ/post-bedtools.py -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    ## 4.4 Process annotated flanks (TSS-filtering)
    try:
        print('\n\n MAJIQ: Including Constiutively Spliced Exon-flanks \n\n', flush=True)
        os.system(f"python 1.1_MAJIQ/combine_AS_CS_flanks.py -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)   

    ## 4.5 Annotate TSS-filtered flanks with peaks
    try:
        print('\n\n MAJIQ: Annotating Candidate Flanks with DHM Score \n\n', flush=True)
        os.system(f'python 2_MANorm/DHM_flanks_MAJIQ.py -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)


    # STEP 5: Process DEXSEQ output and annotate its exon flanks with peaks

    ## 5.1 Get deu exons, TSS-filter flanks of these exons    
    try:
        print('\n\n DEXSEQ: Analysing DEXSEQ Output & Annotating Exons with DEU Score \n\n', flush=True)
        os.system(f'python 1.3_DEXSEQ/post_dexseq.py {output_dir} -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)   

    ## 5.2 Annotate TSS-filtered flanks with peaks
    try:
        print('\n\n DEXSEQ: Annotating Candidate Exons with DHM Score \n\n', flush=True)
        os.system(f'python 2_MANorm/DHM_flanks_DEXSEQ.py -p {proc} >> {log_file_name} 2>&1' )
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1) 

    ## STEP 6: Process RMATS output and annotate its exon flanks with peaks

    ## 6.1 Extract skipped exons   
    try:
        print('\n\n RMATS: Analysing RMATS Output - Skipped Exons \n\n', flush=True)
        os.system(f'python 1.2_RMATS/get_SE.py {output_dir} -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1) 

    ## 6.2 Extract mutually exclusive exons
    try:
        print('\n\n RMATS: Analysing RMATS Output - Mutually Exclusive Exons \n\n', flush=True)
        os.system(f'python 1.2_RMATS/get_MXE.py {output_dir} -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1) 

    ## 6.3 Get flanks of alternative exons
    try:
        print('\n\n RMATS: Combining DEU Scores of Skipped and Mutally Exclusive Exons \n\n', flush=True)
        os.system(f'python 1.2_RMATS/combine_AS_exons.py -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1) 


    ## 6.4 annotate TSS-filtered flanks with peaks
    try:
        print('\n\n RMATS: Annotating Candidate Exons with DHM Score \n\n', flush=True)
        os.system(f'python 2_MANorm/DHM_flanks_RMATS.py -p {proc} >> {log_file_name} 2>&1')
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1) 

    ## STEP 7: DEU - DHM Correlation
    try:
        print('\n\n Computing DEU-DHM Correlation \n\n', flush=True)
        os.system(f"python 3_Episplicing/correlation_plot.py -p {proc} >> {log_file_name} 2>&1")
        print('\n\n----------- DONE -----------', flush=True)
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    # STEP 8: Move files generated from current pipeline run to
    move_dirs(output_dir, proc)
    print(f'Execution of   pipeline finished successfully!')


active_processes = []

# Start the processes
for pr, out_dir in zip(list_of_processes, output_dirs):
    proc = Process(target=master_function, args=(pr, out_dir,), daemon=True)
    active_processes.append(proc)
    proc.start()

# Wait for the processes to complete and free resources
for proc in active_processes:
        proc.join()