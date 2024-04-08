import pandas as pd
import os
import sys

def exon_flanks(i):


    input_files = []

    name = 'epi' if i == '0' else 'nonepi'
    op_dir = "0_Files/Post-processing" if i =='0' else "Post-processing"

    flanks = pd.read_csv(f'{op_dir}/{name}_flanks.bed', delimiter='\t', header=None)
    flanks.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']
    flanks['flanks'] = flanks[['flank_start', 'flank_end']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    
    # Get input sequences for RBPmap
    n = 0
    i = 1
    while n < len(flanks):
            file_name = f'../RBPmap/rbp_input_{name}{i}.csv'
            input_files.append(file_name)
            if n+5000 <= len(flanks):
                    flanks[['chr', 'flanks', 'strand']].iloc[n:n+5000].to_csv(file_name, index=False, sep=':', header=False)
                    n += 5000
            else:
                    flanks[['chr', 'flanks', 'strand']].iloc[n:len(flanks)+1].to_csv(file_name, index=False, sep=':', header=False)
                    break
            i += 1


    # # get flanks for feature matrix preparation step
    flanks[['gene_name', 'flanks']].to_csv(f'{op_dir}/query_flanks_{name}.csv', sep='\t', index=False)

    input_files = [os.getcwd() + '/' + file for file in input_files]
    with open('../RBPmap/input.txt', 'w') as f:
        for item in input_files:
            f.write("%s\n" % item)



if __name__ == "__main__":
     
    exon_flanks(sys.argv[1])