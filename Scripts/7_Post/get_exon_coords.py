import pandas as pd
import os

op_dir = "0_Files/Post-processing"
for i in range(2):
    if i == 0:
        output = f"{op_dir}/epi_flanks_meta.bed"
        name = "epi"
        os.system(f"bedtools intersect -wo -a 0_Files/Post-processing/epi_flanks.bed -b 0_Files/flanks_meta.bed -s | awk -F'\t' '$NF == 400' | sort | uniq > {output}")
        merged_df = pd.read_csv(output, delimiter='\t', header=None)[[0,1,2,3,4,5,6,12,7]]

    else:
        output = f"{op_dir}/nonepi_flanks_meta.bed"
        name = "nonepi"
        os.system(f"bedtools intersect -wo -a 0_Files/Post-processing/nonepi_flanks.bed -b 0_Files/flanks_meta.bed -s | awk -F'\t' '$NF == 400' | sort | uniq > {output}")
        merged_df = pd.read_csv(output, delimiter='\t', header=None)[[0,1,2,3,4,5,6,12,7]]
    
    merged_df.columns = ['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'exon_coords', 'type']

    # drop duplicates
    merged_df = merged_df[merged_df['exon_coords'].notna()]
    merged_df = merged_df.drop_duplicates()

    # Get indiv exon coords
    merged_df[['exon_start', 'exon_end']] = merged_df['exon_coords'].str.split('-', expand=True)
    merged_df['exon_start'] = merged_df['exon_start'].astype(int)
    merged_df['exon_end'] = merged_df['exon_end'].astype(int)

    # Get exon length
    merged_df['exon_bp'] = merged_df['exon_end'] - merged_df['exon_start']

    # set default value
    merged_df['drop'] = False
    duplicates = merged_df.duplicated(subset=['chr', 'flank_start', 'flank_end', 'gene_name'],  keep=False)

    ## get longest version of a3ss/a3SS exons AS for RBPmap
    merged_df.loc[duplicates, 'drop'] = merged_df[duplicates].groupby(['exon_start', 'gene_name'])['exon_bp'].rank(method='dense', ascending=False) != 1
    merged_df.loc[duplicates, 'drop'] = merged_df[duplicates].groupby(['exon_end', 'gene_name'])['exon_bp'].rank(method='dense', ascending=False) != 1

    print('# merged flanks ', len(merged_df))

    merged_df = merged_df[~merged_df['drop']]
    flanks_df = merged_df[['chr', 'flank_start', 'flank_end', 'feature', 'score', 'strand', 'gene_name', 'type']].drop_duplicates()
    flanks_df.to_csv(output,sep='\t', header=False, index=False)
    merged_df[['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']].drop_duplicates()[['chr', 'exon_start', 'exon_end', 'feature', 'score', 'strand', 'gene_name']].to_csv(f'{op_dir}/{name}_exons.bed', sep='\t', header=False, index=False)

    print('# flanks with exon coords', len(flanks_df))

