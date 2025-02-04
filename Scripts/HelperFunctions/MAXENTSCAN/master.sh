awk -F'\t' '{mid = int(($2 + $3) / 2); print $0 "\t" mid}' epi_flanks.bed > temp.bed && mv temp.bed flanks.bed

bash 1_get_exon_coords.sh flanks.bed

bash 2_merge.sh flanks.bed 11 > temp.bed && mv temp.bed flanks.bed 

bash 3_get_indiv_exons.sh flanks.bed | sort | uniq > temp.bed && mv temp.bed flanks.bed 

bash 4_get_splice_sites.sh

#### Run MAXENTSCAN ###

bash 5_process_scores.sh
