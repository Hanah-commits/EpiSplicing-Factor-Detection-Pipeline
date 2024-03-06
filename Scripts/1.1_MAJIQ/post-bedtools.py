import pandas as pd
import os
import json


flank_lens = [50, 100, 200]
junctions = pd.read_csv('0_Files/MAJIQ/majiq_junctions.csv', delimiter='\t')
with open('paths.json') as f:
    d = json.load(f)

fasta = d['Reference fasta']
ref_genome= fasta+".fai"

annot = []
for length in flank_lens:

        # adjust flank length to 200bp +- without exceeding chromosome bounds
        if length < 200:
                file = "majiq_flanks" + str(length) + ".bed"
                adjust_size = str(200 -length)

                # separate start,stop flank coords
                os.system("cut -f 1-3 0_Files/MAJIQ/"+  file +" > 0_Files/coords.bed")

                # adjust flank boundaries
                os.system("bedtools slop -i 0_Files/coords.bed" + " -g " + ref_genome + " -b " + adjust_size + " > 0_Files/coords_adjusted.bed")

                # replace flank coords with adjusted coords
                os.system("awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' 0_Files/coords_adjusted.bed 0_Files/MAJIQ/" + file + " > 0_Files/adjusted_flanks.bed")
                os.system("awk 'FNR==NR{a[NR]=$3;next}{$3=a[FNR]}1' 0_Files/coords_adjusted.bed 0_Files/adjusted_flanks.bed > 0_Files/adjusted.bed")
                os.system("sed 's/ /\t/g' 0_Files/adjusted.bed > 0_Files/MAJIQ/"+file)

                # remove intermediate files
                os.system("rm 0_Files/coords*.bed")
                os.system("rm  0_Files/adjusted*.bed")


        flanks = pd.read_csv('0_Files/MAJIQ/majiq_flanks' + str(length) + '.bed', delimiter='\t', header=None)

        # drop flanks that have no junction
        ##chr1   10324545   10325245  Exon  .  +  ENSG00000054523.17  chr1   10324895   10324896  flank  .  +  ENSG00000054523.17 
        ##chr10   100041843       100042543       Exon    .       -     ENSG00000274847.1      .       -1      -1      .       -1      .       ENSG00000274847.1 

        flanks = flanks[flanks[8] != -1]

        # merge the flanks df with the jns df
        flanks[14] = flanks[[1, 2]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        flanks.drop([3, 4, 5, 7, 9, 10, 11, 12], axis=1, inplace=True)
        flanks.columns = ['seqid', 'start', 'stop', 'gene_id', 'junction0', 'gene_id2', 'flanks']

        ## drop flanks with wrong junctions mapped to them (multi-geneic exons)
        flanks = flanks[~(flanks['gene_id'] != flanks['gene_id2'])]
        del (flanks['gene_id2'])


        junctions['index'] = junctions.index
        flank_jns = pd.merge(junctions, flanks, on=['junction0', 'seqid', 'gene_id'])

        flank_jns.sort_values(['index'], inplace=True)
        flank_jns = flank_jns.reset_index(drop=True)
        del (flank_jns['index'])
        flank_jns['dpsi_'+str(length)] = pd.to_numeric(flank_jns['mean_dpsi_per_lsv_junction'])

        # # get all junctions that belong to each flank
        flank_jns_group = flank_jns.groupby(['flanks', 'gene_id'])['dpsi_'+str(length)] \
                .apply(lambda val: ','.join(str(v) for v in val)).reset_index()
        
        # unqiue index
        flank_jns_group['idx'] = flank_jns_group['gene_id'] + flank_jns_group['flanks']

        annot.append(flank_jns_group)


for a in annot:
        a.set_index('idx',inplace=True)

df = pd.concat(annot,axis=1,sort=False).reset_index()
df.dpsi_50.fillna(df.dpsi_100, inplace=True)
df.dpsi_50.fillna(df.dpsi_200, inplace=True)
flank_jns_group = df.drop(['dpsi_100', 'dpsi_200'], axis=1)
flank_jns_group = flank_jns_group.iloc[:, [-1, -2, 3]]
flank_jns_group.columns = ['gene_id', 'flanks', 'mean_dpsi_per_lsv_junction']

# FILTER 3: if flank has 1+ junctions, keep junction with highest dPSI value
flank_jns_group['max_dPSI'] = flank_jns_group['mean_dpsi_per_lsv_junction']\
        .apply(lambda x: min(map(float, x.split(','))) if isinstance(x, str) else x) # string -> list of strings -> list of floats -> max float

# # get the corresponding junction for each flank's max dPSI value
flank_jns_group = pd.merge(flank_jns_group[['flanks', 'max_dPSI', 'gene_id']], flank_jns, on=['flanks', 'gene_id'], how='inner')

# Group by 'flanks' and 'gene_id', then apply filtering using a lambda function
flank_jns_group = (flank_jns_group.groupby(['flanks', 'gene_id'])
                              .apply(lambda x: x[x['mean_dpsi_per_lsv_junction'] == x['max_dPSI']])
                              .reset_index(drop=True))

# FILTER 4: If flank has 1+ junctions with same max dPSI value, keep one
flank_jns_group.drop_duplicates(subset=['flanks', 'gene_id'], keep='first', inplace=True)

# # bookkeeping
del(flank_jns_group['max_dPSI'])
flank_jns_group = flank_jns_group[['gene_id', 'lsv_id', 'seqid', 'junction0', 'mean_dpsi_per_lsv_junction',
        'probability_changing', 'flanks', 'start', 'stop', 'strand']]

# Get all filtered flanks
flank_jns_group.drop_duplicates().to_csv('0_Files/MAJIQ/all_flanks.csv', sep='\t', index=False)


# Get dPSI values of filtered flanks
flank_jns_group[['seqid', 'strand', 'start', 'stop', 'gene_id', 'mean_dpsi_per_lsv_junction']].drop_duplicates().to_csv('0_Files/MAJIQ/Filtered_dPSI.csv', index=False, sep='\t')