import os
import pandas as pd


def get_seqid_strand():

    dju = pd.read_csv('0_Files/majiq_junctions.csv', delimiter='\t')
    non_dju = pd.read_csv('0_Files/majiq_junctions_control.csv', delimiter='\t')
    non_dju_genes = list(set(non_dju.gene_id.values.tolist()))
    dju_genes = list(set(dju.gene_id.values.tolist()))
    common_genes = list(set(non_dju_genes) & set(dju_genes))


    ## remove genes with no non-DJU events to use as control from further analysis
    if len(dju_genes) > len(common_genes):
        dju = dju[dju['gene_id'].isin(common_genes)]
        # prepare bedtools input again
        keep_cols = ['seqid', 'junction0', 'strand']
        majiq_bed = dju[keep_cols]
        majiq_bed = majiq_bed.drop_duplicates()
        # to fit bedtools input requirements
        majiq_bed['junction1'] = pd.to_numeric(majiq_bed['junction0']) + 1
        majiq_bed['feature'] = "flank"
        majiq_bed['score'] = "."
        # rearrange
        majiq_bed = majiq_bed[['seqid', "junction0", "junction1", "feature", "score", "strand"]]
        majiq_bed.to_csv('0_Files/majiq.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
        dju.to_csv('0_Files/majiq_junctions.csv', index=False, sep='\t', header=True)


    ## only non-dju events for which dju events are found
    non_dju = non_dju[non_dju['gene_id'].isin(common_genes)]

    # create empty columns
    non_dju['seqid'] = [0] * len(non_dju)
    non_dju['strand'] = ['-'] * len(non_dju)

    for gene in common_genes:
        # get chr and strand details from dju event
        seqid = dju[dju['gene_id']==gene]['seqid'].values.tolist()[0]
        strand = dju[dju['gene_id']==gene]['strand'].values.tolist()[0]
        # add chr and strand details to non-dju event
        non_dju.loc[non_dju['gene_id'] == gene, 'seqid'] = seqid
        non_dju.loc[non_dju['gene_id'] == gene, 'strand'] = strand

    keep_cols = ['seqid', 'junction0', 'strand']
    majiq_bed = non_dju[keep_cols]
    majiq_bed = majiq_bed.drop_duplicates()
    # to fit bedtools input requirements
    majiq_bed['junction1'] = pd.to_numeric(majiq_bed['junction0']) + 1
    majiq_bed['feature'] = "flank"
    majiq_bed['score'] = "."
    # rearrange
    majiq_bed = majiq_bed[['seqid', "junction0", "junction1", "feature", "score", "strand"]]
    majiq_bed.to_csv('0_Files/majiq_control.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
    non_dju.to_csv('0_Files/majiq_junctions_control.csv', index=False, sep='\t', header=True)


def annotate_control():
    
    flanks = [50,100,200]
    for flank in flanks:
        os.system('bedtools intersect -loj -s -a ' + '0_Files/flanks' + str(flank) +'.bed -b 0_Files/majiq_control.bed | sort | uniq > 0_Files/majiq_flanks' + str(flank) + '_control.bed')


if __name__ == "__main__":
    get_seqid_strand()
    annotate_control()