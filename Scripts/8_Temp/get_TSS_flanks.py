
import json
import os

with open('paths.json') as f:
    d = json.load(f)

ref = d['Reference GFF3']
fasta = d['Reference fasta']

# genome file
ref_genome= fasta+".fai"


# exon boundary external flanks
os.system("bedtools flank -i 0_Files/TSS.bed -g " + ref_genome + " -b "+ "200" + " > 0_Files/flanks.bed" )

# separate start,stop flank coords
os.system("sed -n 'n;p' 0_Files/flanks.bed > 0_Files/stop.bed")
os.system("sed -n 'p;n' 0_Files/flanks.bed > 0_Files/start.bed")

# exon boundary internal flanks
os.system("bedtools slop -i 0_Files/start.bed -g " + ref_genome + " -l 0 -r " + "200" + " > 0_Files/start_flanks.bed")
os.system("bedtools slop -i 0_Files/stop.bed -g " + ref_genome +" -l " + "200" + " -r 0 > 0_Files/stop_flanks.bed")

# combine start,stop flank coords
# ( 2 * num exons != num flanks: exons can be A3SS/A5SS versions, can have duplicate flanks)
os.system("paste -d'\n' 0_Files/start_flanks.bed 0_Files/stop_flanks.bed | sort | uniq > 0_Files/Post-processing/TSS_flanks" + "200" + ".bed")

# remove intermediate files
os.system("rm 0_Files/start*.bed")
os.system("rm  0_Files/stop*.bed")
os.system("rm 0_Files/flanks.bed")