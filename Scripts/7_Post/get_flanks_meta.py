import json
import os

with open('paths.json') as f:
    d = json.load(f)


ref_genome = d['Reference fasta'] + ".fai"

# Prepare flanks keeping exon coord info
os.system('awk \'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $2"-"$3, $6,$7}\' 0_Files/exon_coords.bed > 0_Files/exons_meta.bed')

# exon boundary external flanks
os.system("bedtools flank -i 0_Files/exons_meta.bed -g " + ref_genome + " -b 200 > 0_Files/flanks.bed" )

# separate start,stop flank coords
os.system("sed -n 'n;p' 0_Files/flanks.bed > 0_Files/stop.bed")
os.system("sed -n 'p;n' 0_Files/flanks.bed > 0_Files/start.bed")

# exon boundary internal flanks
os.system("bedtools slop -i 0_Files/start.bed -g " + ref_genome + " -l 0 -r 200 > 0_Files/start_flanks.bed")
os.system("bedtools slop -i 0_Files/stop.bed -g " + ref_genome +" -l 200 -r 0 > 0_Files/stop_flanks.bed")

# combine start,stop flank coords
# ( 2 * num exons != num flanks: exons can be A3SS/A5SS versions, can have duplicate flanks)
os.system("paste -d'\n' 0_Files/start_flanks.bed 0_Files/stop_flanks.bed | sort | uniq > 0_Files/flanks_meta.bed")

# remove intermediate files
os.system("rm 0_Files/start*.bed")
os.system("rm  0_Files/stop*.bed")
os.system("rm 0_Files/flanks.bed")
os.system("rm 0_Files/exons_meta.bed")