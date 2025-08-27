import pandas as pd, numpy as np
import re
import os
import glob


def post_rbp(rbp_num):

    # get list of RBPs
    os.system(f'cp HelperFunctions/rbps_{str(rbp_num)}.txt ../RBPmap_{str(rbp_num)}/rbps.txt')

    for type in ['epi', 'nonepi', 'epi_nonspliced']:
        opdir = '0_Files/Post-processing'
                    
        rbps = open(f'../RBPmap_{str(rbp_num)}/rbps.txt', 'r') #'HelperFunctions/rbps_{rbp_num}.txt'
        proteins = list(set([rbp for rbp in rbps.read().split('\n') if len(rbp)]))

        # # proteins = ['RBM15'] #if single rbp
        try:
            filename=f'../RBPmap_{str(rbp_num)}/resultsrbp_input_{type}1.csv*/All_Predictions.txt' # local RBPmap/old versions
            filename = [file for file in glob.glob(filename)][0]
        except:
            filename=f'../RBPmap_{str(rbp_num)}/resultsrbp_input_{type}1.csv*/All_Predictions.csv' # webserver/ new version
            filename = [file for file in glob.glob(filename)][0]
        
        if 'Predictions.csv' in filename:
            parse_csv_file(filename, proteins, type, opdir)
        else:
            parse_txt_file(filename, proteins, type, rbp_num, opdir)
        

def parse_csv_file(filename, proteins, type, opdir):
    coor_scores = pd.DataFrame(columns=proteins)
    coordinate = ""
    score_dict = {}

    with open(filename, 'r') as file:

        for line in file:
            line = line.strip()

            # skip empty lines and metadata
            if not line or line.startswith(("Predictions for job", "Calculation parameters", "Genome:", "Selected motifs:", "Stringency level:", "Conservation filter:", "Protein")):
                continue

            if line.startswith("chr") and ":" in line and "-" in line:
                coordinate = line
            elif line == 'No motifs found.':
                score_dict[coordinate] = {}
                score_dict[coordinate][proteins[0]] = 0.0 #pseudo score

            fields = [f.strip() for f in line.split(",")]
            if len(fields) != 7:
                continue

            raw_protein = fields[0]
            z_score_str = fields[5]
            protein_match = re.search(r"(?:User_)?([a-zA-Z0-9]+)(?:\()?", raw_protein) #"U2AF2(Hs/Mm)", "User_U2AF2"
            # print(protein_match)
            if not protein_match or protein_match.group(1) not in proteins:
                continue
            protein = protein_match.group(1)
            

            try:
                z_score = float(z_score_str)
            except ValueError:
                continue

            # keep max binding score of RBP in current flanks
            if coordinate and protein:
                if coordinate not in score_dict:
                    score_dict[coordinate] = {}
                current = score_dict[coordinate].get(protein, 0.0)
                score_dict[coordinate][protein] = max(current, z_score)

        for coord, protein_scores in score_dict.items():
            for prot, val in protein_scores.items():
                coor_scores.loc[coord, prot] = val

        coor_scores.infer_objects(copy=False).fillna(0.0, inplace=True)
        coor_scores.index = coor_scores.index.astype(str)

        # # Split coords into columns
        coor_scores['coord'] = coor_scores.index
        coor_scores[['chr', 'start_end', 'strand']] = coor_scores['coord'].str.split(':', expand=True)
        coor_scores[['flank_start', 'flank_end']] = coor_scores['start_end'].str.split('-', expand=True)
        #drop index and cols, rearrange cols
        coor_scores = coor_scores.drop(columns=['coord','start_end'])
        coor_scores = coor_scores[['chr', 'flank_start', 'flank_end', 'strand'] + proteins]
        coor_scores.reset_index(drop=True, inplace=True)
        coor_scores.to_csv(f"{opdir}/FilteredZscores_{type}.csv", sep=',', index=None)


def parse_txt_file(filename, proteins, type, rbp_num, opdir):
    parsed_data = []

    # open and read the input file
    with open(filename, 'r') as file:
        current_region = None
        current_strand = None
        current_protein = None
        region_data = {protein: [] for protein in proteins}  # (zscore, p_value)

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
                            zscores, p_values = zip(*z_p_tuples)
                            # select max Z-score 
                            selected_zscore = np.max(np.array(zscores, dtype=float)) 
                            row_data.append(selected_zscore)

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
                        zscore = str(line_parts[-2]) # save most impt binding score
                        p_value = float(line_parts[-1])
                        region_data[current_protein].append((zscore, p_value))
                    except ValueError:
                        continue

        # append the last region data
        if current_region:
            row_data = []
            for protein in proteins:
                z_p_tuples = region_data[protein]

                if z_p_tuples:
                    zscores, p_values = zip(*z_p_tuples)
                    # select max Z-score 
                    selected_zscore = np.max(np.array(zscores, dtype=float)) 
                else:
                    row_data.append(0.0)

            parsed_data.append(current_region + [current_strand] + row_data)


    # # Create DataFrame from parsed data
    df = pd.DataFrame(parsed_data, columns=['chr', 'flank_start', 'flank_end', 'strand'] + proteins)
    df.to_csv(f"{opdir}/FilteredZscores_{type}.csv", sep=',', index=None)


def feature_matrix_1():

    exons_files = []
    Zscore_files = []
    types = ['epi', 'nonepi', 'epi_nonspliced']
    for type in types:
        exons_files.append(f'0_Files/Post-processing/{type}_flanks.bed')
        Zscore_files.append(f'0_Files/Post-processing/FilteredZscores_{type}.csv')

    for i in range(len(exons_files)):
        exons = pd.read_csv(exons_files[i], delimiter='\t', header=None)
        exons.columns = ['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name', 'type']
        rbp = pd.read_csv(Zscore_files[i], delimiter=',').fillna(0.0)
        features = pd.concat([exons, rbp], axis=1)

        # double-check rbpmap results parsing
        if len(features[features.exon_start != features.flank_start]) == 0 and len(features[features.exon_end != features.flank_end]) == 0: # check 
            features = features.drop(columns=['flank_start', 'flank_end'])
            features = features.loc[:, ~features.columns.duplicated(keep='first')] # keep first occurence of columns: chr, strand
        else:
            raise ValueError
        
        name = types[i]

        features = features.loc[:, ~features.columns.duplicated()]
        features.to_csv('0_Files/Post-processing/features_' + name + '.csv', sep='\t', index=False)

    # remove intermediate files
    os.system(f'rm 0_Files/Post-processing/FilteredZscores*.csv')


def feature_matix_2(rbp_num):

    types = ['epi', 'nonepi', 'epi_nonspliced']
    labels = ['epigene', 'non-epigene', 'epi_nonspliced_gene']

    # read feature matrices
    feature_dfs = []
    for type in types:
        feature_df = pd.read_csv(f'0_Files/Post-processing/features_{type}.csv', delimiter='\t')
        if type == 'nonepi': # keep nonepi flanks for hms available in current study
            feature_df = feature_df[feature_df['type'].apply(lambda x: any(item in ['H3K27ac', 'H3K27me3', 'H3K9me3', 'H3K4me3', 'H3K36me3'] for item in x.split(',')))]
        feature_dfs.append(feature_df)

    # construct feature matrics of each comparison pair
    for i in range(len(types)):
        for j in range(len(types)):
            if i >= j:
                continue
            
            features1 = feature_dfs[i]
            features2 = feature_dfs[j]

            # label exon classes
            features1['labels'] = labels[i]
            features2['labels'] = labels[j]

            # merge into single matrix
            all_features = pd.concat([features1, features2], axis=0)

            filtered_features = []
            for hm in ['H3K27ac', 'H3K27me3', 'H3K9me3', 'H3K4me3', 'H3K36me3']:
                print('\n',hm)
                temp_features = all_features[all_features['type'].apply(lambda x: hm in x.split(','))]
                features1 = temp_features[temp_features.labels == labels[i]]
                features2 = temp_features[temp_features.labels != labels[i]]
                
                # remove genes with both labels
                common_genes = list(set(features1.gene_name.values.tolist()) & set(features2.gene_name.values.tolist()))
                if len(common_genes):
                    temp_features = temp_features[~((temp_features.gene_name.isin(common_genes)) & (temp_features.labels == labels[j]))]
                    temp_features = temp_features[~((temp_features.gene_name.isin(common_genes)) & (temp_features.labels == labels[i]))]

                print(f'{labels[i].title()}s:', len(set(temp_features[temp_features.labels==labels[i]].gene_name.values)), 'Flanks:', len(temp_features[temp_features.labels==labels[i]]))
                print(f'{labels[j].title()}s:', len(set(temp_features[temp_features.labels==labels[j]].gene_name.values)), 'Flanks:', len(temp_features[temp_features.labels==labels[j]]))

                filtered_features.append(temp_features)

            all_features = pd.concat(filtered_features, axis=0).drop_duplicates()
            all_features.drop(['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name'], axis=1, inplace=True)


            col = all_features.pop("labels")
            all_features.insert(0, col.name, col)

            all_features.to_csv(f'0_Files/Post-processing/features_all_{types[i]}_vs_{types[j]}_{rbp_num}.csv', sep='\t', index=False)
            print( '\n-----------------------\n')

    # remove intermediate files
    for type in types:
        os.system(f'rm 0_Files/Post-processing/features_{type}.csv')


if __name__ == "__main__":
    
    # prep feature zscore matrix -132 RBPS
    post_rbp(132)
    feature_matrix_1()
    feature_matix_2(132)
    
    # prep feature zscore matrix - 47 RBPS
    post_rbp(47)
    feature_matrix_1()
    feature_matix_2(47)