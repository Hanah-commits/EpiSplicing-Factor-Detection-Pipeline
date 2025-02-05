import os, sys
import subprocess
import pandas as pd


def annotate_eclip_peaks(eclip_directory, flanks_dir):

    # Get list of files containing 'HepG2.bed' or 'K562.bed'
    bed_files = sorted([f for f in os.listdir(eclip_directory) if f.endswith('.bed') and ('HepG2' in f or 'K562' in f)])
    flanks_file = f'{flanks_dir}/Post-processing/epi_flanks.bed'
    annotated_flanks_file = f'{flanks_dir}/Post-processing/epi_flanks_annotated.bed'
    subprocess.run(f"cp {flanks_file} {annotated_flanks_file}", shell=True, check=True)

    # Iterate over each file
    for bed_file in bed_files:
        bed_file_path = os.path.join(eclip_directory, bed_file)
        print(bed_file)
        
        # Command to run bedtools intersect
        intersect_command = f"bedtools intersect -loj -s -a {annotated_flanks_file} -b {bed_file_path} | sort | uniq > bedtools_output.bed"
        subprocess.run(intersect_command, shell=True, check=True)
        
        # Load the temp_output.bed file into a DataFrame
        temp_df = pd.read_csv('bedtools_output.bed', sep='\t', header=None)
        temp_df[10] = temp_df[10].replace(-1, None)

        # Modify DataFrame based on column 11 and 12
        temp_df[8] = temp_df.apply(
        lambda row: f"{row[8]},{row[12].replace('_IDR', '')}" if row[8] != 'eCLIP' and row[10] is not None 
        else row[12].replace('_IDR', '') if row[8] == 'eCLIP' and row[10] is not None
        else row[8],  # No change if column 10 is None
        axis=1)

        # Drop all columns after the 8th column
        temp_df = temp_df.iloc[:, :9]

        # Drop duplicates
        temp_df.groupby([0,1,2,6])[8].apply(','.join).reset_index()
        temp_df[8] = temp_df[8].apply(lambda x: ','.join(sorted(set(x.split(','))))) # remove duplicate eCLIP annotation

        # Save modified DataFrame back to temp_output.bed
        temp_df.drop_duplicates().to_csv('bedtools_output.bed', sep='\t', header=False, index=False)
        
        # Replace file1.bed with the updated temp_output.bed
        subprocess.run(f"mv bedtools_output.bed {annotated_flanks_file}", shell=True, check=True)


def get_deu_dhm_info(op_dir):

    epi_flanks = pd.read_csv(f'{op_dir}/Post-processing/epi_flanks_annotated.bed', delimiter='\t', header=None, names=['chr', 'flank_start', 'flank_stop', 'feature', 'score', 'strand', 'gene', 'DHM', 'eCLIP'])
    exon_coords = pd.read_csv(f'{op_dir}/RMATS/rmats_exons_coords.bed', delimiter='\t', names=['chr', 'exon_start', 'exon_stop', 'feature', 'score', 'strand', 'gene', 'dPSI'])

    epi_flanks['exon_coord0'] = epi_flanks['flank_start'] + 200


    # STEP 1: Get start/stop coords of exons of flanks
    epi_flanks['coord_tuple'] = epi_flanks['exon_coord0'].apply(
        lambda coord: (
            exon_coords[exon_coords['exon_start'] == coord][['exon_start', 'exon_stop']].iloc[0].tolist() 
            if not exon_coords[exon_coords['exon_start'] == coord].empty else 
            exon_coords[exon_coords['exon_stop'] == coord][['exon_start', 'exon_stop']].iloc[0].tolist()
            if not exon_coords[exon_coords['exon_stop'] == coord].empty else None
        )
    )

    epi_flanks[['exon_start', 'exon_stop']] = pd.DataFrame(epi_flanks['coord_tuple'].tolist(), index=epi_flanks.index)
    epi_flanks = epi_flanks.drop(columns=['coord_tuple', 'exon_coord0']) # delete unnecc cols

    # STEP 2: Get dPSI and inclusion status of exons
    SE = pd.read_csv(f'{op_dir}/RMATS/SE_exons.csv', delimiter='\t')
    MXE = pd.read_csv(f'{op_dir}/RMATS/MXE_exons.csv', delimiter='\t')
    SE_MXE_exons =  pd.concat([SE, MXE], ignore_index=True)


    epi_flanks['dPSI'] = epi_flanks.apply(
    lambda row: SE_MXE_exons[SE_MXE_exons['exon_coord0'].isin([row['exon_start'], row['exon_stop']])]['dPSI'].iloc[0] 
    if not SE_MXE_exons[SE_MXE_exons['exon_coord0'].isin([row['exon_start'], row['exon_stop']])].empty 
    else None,
    axis=1
    )

    epi_flanks['DEU'] = epi_flanks.apply(
    lambda row: SE_MXE_exons[SE_MXE_exons['exon_coord0'].isin([row['exon_start'], row['exon_stop']])]['InclusionStatus'].iloc[0] 
    if not SE_MXE_exons[SE_MXE_exons['exon_coord0'].isin([row['exon_start'], row['exon_stop']])].empty 
    else None,
    axis=1
    )

    # STEP 3: Get Mvalue and DHM status of flanks
    hms = ['H3K27ac', 'H3K4me3', 'H3K36me3', 'H3K9me3']
    for hm in hms:
        epi_flanks[f'M_value_{hm}'] = None

        dhm_file = pd.read_csv(f'{op_dir}/MANorm/{hm}_all_exons.bed', delimiter='\t', header=None)
        dhm_file.drop([12, 13, 15, 16], axis=1, inplace=True)
        dhm_file.columns = ['chr', "exon_start", "exon_stop", "feature", "score", "strand", "geneSymbol", 'chr_2', 'peak_start', 'peak_end', 'summit', 'M_value', 'signal_status']
        dhm_file = dhm_file[['chr', "exon_start", "exon_stop",'M_value','signal_status']]

        epi_flanks[f'M_value_{hm}'] = epi_flanks.apply(
        lambda row: dhm_file[(dhm_file['exon_start'].isin([row['exon_start'], row['exon_stop']])) | (dhm_file['exon_stop'].isin([row['exon_start'], row['exon_stop']]))]['M_value'].iloc[0] 
        if not dhm_file[(dhm_file['exon_start'].isin([row['exon_start'], row['exon_stop']])) | (dhm_file['exon_stop'].isin([row['exon_start'], row['exon_stop']]))].empty 
        else None,
        axis=1
        ) 

        epi_flanks[f'Signal_{hm}'] = epi_flanks.apply(
        lambda row: dhm_file[(dhm_file['exon_start'].isin([row['exon_start'], row['exon_stop']])) | (dhm_file['exon_stop'].isin([row['exon_start'], row['exon_stop']]))]['signal_status'].iloc[0].split('_peak_unique')[0] 
        if not dhm_file[(dhm_file['exon_start'].isin([row['exon_start'], row['exon_stop']])) | (dhm_file['exon_stop'].isin([row['exon_start'], row['exon_stop']]))].empty 
        else None,
        axis=1
        )

    epi_flanks.drop_duplicates().to_csv(f'{op_dir}/Post-processing/epi_flanks_annotated.bed', sep='\t',index=False)


def process_annotated_eclip(op_dir):
    
    df = pd.read_csv(f'{op_dir}/Post-processing/epi_flanks_annotated.bed', delimiter='\t')
    
    os.makedirs(f"{op_dir}/Post-processing/eclip/", exist_ok=True)  # Create the directory if it doesn't exist

    # split column values to multiple lines
    df['eCLIP'] = df['eCLIP'].str.replace(".,", "")
    df = df.assign(eCLIP=df['eCLIP'].str.split(','))
    df = df.explode(['eCLIP']).reset_index(drop=True)
    df = df.assign(DHM=df['DHM'].str.split(','))
    df = df.explode(['DHM']).reset_index(drop=True)

    hms = [ "H3K27ac",
        "H3K4me3",
        "H3K9me3",
        "H3K36me3"]

    exact_df = []
    exact_opp = []
    partial_dhm = []
    partial_rbp = []

    for hm in hms:
        
        eclip = df[df['DHM'].str.contains(hm)]

        cells = ['HepG2', 'K562']
        for cell in cells:

            # EXACT MATCHES
            filtered = eclip[eclip.apply(lambda row: cell in row['eCLIP'] and cell in row[f'Signal_{hm}'] and cell == row['DEU'], axis=1)]
            exact_df.append(filtered)

            # EXACT OPP MATCHES
            filtered = eclip[eclip.apply(lambda row: cell in row['eCLIP'] and cell in row[f'Signal_{hm}'] and cell != row['DEU'], axis=1)]
            exact_opp.append(filtered)

            # PARTIAL MATCHES - DEU == DHM != eCLIP
            opp_cell = [c for c in cells if c != cell][0]
            filtered = eclip[eclip.apply(lambda row: opp_cell in row['eCLIP'] and cell in row[f'Signal_{hm}'] and cell == row['DEU'], axis=1)]
            partial_dhm.append(filtered)

            # PARTIAL MATCHES - DEU != DHM == eCLIP
            filtered = eclip[eclip.apply(lambda row: cell in row['eCLIP'] and opp_cell in row[f'Signal_{hm}'] and cell == row['DEU'], axis=1)]
            partial_rbp.append(filtered)

        
        
    exact_match = pd.concat(exact_df, ignore_index=True)#
    exact_opp_match = pd.concat(exact_opp, ignore_index=True)
    partial_dhm_match = pd.concat(partial_dhm, ignore_index=True)
    partial_rbp_match = pd.concat(partial_rbp, ignore_index=True)


    # Write to output file
    exact_match.to_csv(f'{op_dir}/Post-processing/eclip/exact_match.tsv', sep='\t', index=False)
    exact_opp_match.to_csv(f'{op_dir}/Post-processing/eclip/exact_opp_match.tsv', sep='\t', index=False)
    partial_dhm_match.to_csv(f'{op_dir}/Post-processing/eclip/partial_dhm_match.tsv', sep='\t', index=False)
    partial_rbp_match.to_csv(f'{op_dir}/Post-processing/eclip/partial_rbp_match.tsv', sep='\t', index=False)


if __name__ == "__main__":
    
    eclip_dir = sys.argv[1]
    op_dir = sys.argv[2]

    # STEP 1: Annotate eclip peaks
    ## Run epigenes_study_master.py with a confg file for HepG2-K562 (SE, MXE: min 10 read support) -> epi_flanks.bed
    # annotate_eclip_peaks(eclip_directory = eclip_dir, flanks_dir = op_dir)
    
    # STEP 2: Get epi exon coords, deu and dhm status
    # get_deu_dhm_info(op_dir)


    # # STEP 2: Split eclip annot by deu, dhm, eclip patterns
    process_annotated_eclip(op_dir)