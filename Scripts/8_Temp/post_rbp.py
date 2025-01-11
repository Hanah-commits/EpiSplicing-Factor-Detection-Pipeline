import pandas as pd
import re
from statsmodels.stats.multitest import multipletests # type: ignore


def post_rbp(rbp_num):
    '''
    process RBPmap output -> obtain motif used in each flank by each RBP
    '''
    

    for type in ['epi', 'nonepi']:
        opdir = '0_Files/Post-processing'
                    
        rbps = open(f'../RBPmap_{str(rbp_num)}/rbps.txt', 'r')
        proteins = list(set([rbp for rbp in rbps.read().split('\n') if len(rbp)]))

        # # proteins = ['RBM15'] #if single rbp

        filename=f'../RBPmap_{str(rbp_num)}/results_rbp_input_{type}1.csv/All_Predictions.txt'

        parsed_data = []

        # open and read the input file
        with open(filename, 'r') as file:
            current_region = None
            current_strand = None
            current_protein = None
            region_data = {protein: [] for protein in proteins}  # (motif, p_value)

            for line in file:
                line = line.strip()

                # get genomic regions and strands
                if re.match(r'^chr', line):
                    if current_region:
                        # process region data before appending
                        row_data = []
                        for protein in proteins:
                            z_p_tuples = region_data[protein]

                            if z_p_tuples:
                                # separate Z-scores and p-values
                                motifs, p_values = zip(*z_p_tuples)

                                # Benjamini-Hochberg correction
                                _, adjusted_pvals, _, _ = multipletests(p_values, method='fdr_bh')

                                # select Z-score corresponding to the smallest adjusted p-value
                                min_p_index = adjusted_pvals.argmin()
                                selected_motif = motifs[min_p_index]

                                # filter by adjusted p-value <= 0.05
                                if adjusted_pvals[min_p_index] <= 0.05:
                                    row_data.append(selected_motif)
                                else:
                                    row_data.append(0.0)
                            else:
                                row_data.append(0.0)

                        # append the region data to parsed_data
                        parsed_data.append(current_region + [current_strand] + row_data)

                    # split line to extract chromosome, start, end, and strand information
                    region_parts = re.split(r'[:-]', line)
                    current_region = [region_parts[0], int(region_parts[1]), int(region_parts[2])]
                    current_strand = region_parts[-1]

                    # Reset region data for new region
                    region_data = {protein: [] for protein in proteins}

                # get protein blocks
                elif line.startswith('Protein:'):
                    if rbp_num == 47:
                        current_protein = line.split('_')[1]
                    else:
                        current_protein = line.split('Protein: ')[1].split('(Hs/Mm)')[0]

                # skip header lines
                elif line.startswith('Sequence Position') or line.startswith('Z-score'):
                    continue

                # extract Z-scores and p-values
                elif line and current_protein in proteins:
                    line_parts = line.split()
                    if len(line_parts) >= 6:  # check there are enough elements
                        try:
                            motif = str(line_parts[-3]) # save most impt binding score
                            p_value = float(line_parts[-1])
                            region_data[current_protein].append((motif, p_value))
                        except ValueError:
                            continue

            # append the last region data
            if current_region:
                row_data = []
                for protein in proteins:
                    z_p_tuples = region_data[protein]

                    if z_p_tuples:
                        motifs, p_values = zip(*z_p_tuples)

                        # apply Benjamini-Hochberg correction
                        _, adjusted_pvals, _, _ = multipletests(p_values, method='fdr_bh')

                        # select Z-score corresponding to the smallest adjusted p-value
                        min_p_index = adjusted_pvals.argmin()
                        selected_motif = motifs[min_p_index]

                        # filter by adjusted p-value <= 0.05
                        if adjusted_pvals[min_p_index] <= 0.05:
                            row_data.append(selected_motif)
                        else:
                            row_data.append(0.0)
                    else:
                        row_data.append(0.0)

                parsed_data.append(current_region + [current_strand] + row_data)


        # # Create DataFrame from parsed data
        df = pd.DataFrame(parsed_data, columns=['chr', 'flank_start', 'flank_end', 'strand'] + proteins)

        # Display the resulting DataFrame
        # df[proteins].to_csv(f"{opdir}/FilteredZscores_{type}.csv", sep=',', index=None)
        df.to_csv(f"{opdir}/FilteredZscores_{type}.csv", sep=',', index=None)


def feature_matrix_1():

    exons_files = ['0_Files/Post-processing/epi_flanks.bed', '0_Files/Post-processing/nonepi_flanks.bed']
    Zscore_files = ['0_Files/Post-processing/FilteredZscores_epi.csv', '0_Files/Post-processing/FilteredZscores_nonepi.csv']

    for i in range(2):
        exons = pd.read_csv(exons_files[i], delimiter='\t', header=None)
        exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name', 'type']

        rbp = pd.read_csv(Zscore_files[i], delimiter=',')

        features = pd.concat([exons, rbp], axis=1)

        # double-check rbpmap results parsing
        if len(features[features.exon_start != features.flank_start]) == 0 and len(features[features.exon_end != features.flank_end]) == 0: # check 
            features = features.drop(columns=['flank_start', 'flank_end'])
            features = features.loc[:, ~features.columns.duplicated(keep='first')] # keep first occurence of columns: chr, strand
        else:
            raise ValueError
        
        name = ''
        if i == 0:
            name = 'epi'
        else:
            name = 'nonepi'

        features = features.loc[:, ~features.columns.duplicated()]
        features.to_csv('0_Files/Post-processing/features_' + name + '.csv', sep='\t', index=False)


def feature_matix_2(rbp_num):
    epi_features = pd.read_csv('0_Files/Post-processing/features_epi.csv', delimiter='\t')
    nonepi_features = pd.read_csv('0_Files/Post-processing/features_nonepi.csv', delimiter='\t')

    # keep nonepi flanks for hms available in current study
    nonepi_features = nonepi_features[nonepi_features['type'].apply(lambda x: any(item in ['H3K27ac', 'H3K27me3', 'H3K9me3', 'H3K4me3', 'H3K36me3'] for item in x.split(',')))]
    # nonepi_features.to_csv('0_Files/Post-processing/features_nonepi.csv', sep='\t', index=False)

    epi_features['label'] = 'epigene'
    nonepi_features['label'] = 'non-epigene'

    all_features = pd.concat([epi_features, nonepi_features], axis=0)

    # remove genes with both labels
    common_genes = list(set(epi_features.gene_name.values.tolist()) & set(nonepi_features.gene_name.values.tolist()))
    all_features = all_features[~((all_features.gene_name.isin(common_genes)) & (all_features.label == 'non-epigene'))]
    all_features = all_features[~((all_features.gene_name.isin(common_genes)) & (all_features.label == 'epigene'))]

    for hm in ['H3K27ac', 'H3K27me3', 'H3K9me3', 'H3K4me3', 'H3K36me3']:
        print('\n',hm)
        temp_features = all_features[all_features['type'].apply(lambda x: hm in x.split(','))]
        print('Epigenes:', len(set(temp_features[temp_features.label=='epigene'].gene_name.values)), 'Flanks:', len(temp_features[temp_features.label=='epigene']))
        print('Nonepigenes:', len(set(temp_features[temp_features.label=='non-epigene'].gene_name.values)), 'Flanks:', len(temp_features[temp_features.label=='non-epigene']))

    all_features.drop(['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name'], axis=1, inplace=True)


    col = all_features.pop("label")
    all_features.insert(0, col.name, col)

    all_features.to_csv(f'0_Files/Post-processing/features_all_{rbp_num}.csv', sep='\t', index=False)



if __name__ == "__main__":
    
    # # prep feature motif matrix -132 RBPS
    post_rbp(132)
    feature_matrix_1()
    feature_matix_2(132)
    
    # prep feature motif matrix - 47 RBPS
    post_rbp(47)
    feature_matrix_1()
    feature_matix_2(47)