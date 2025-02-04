# Generate 3' splice site output
awk '{ if ($6 == "+") { start = $2 - 20; end = $2 + 3; } else { start = $3 - 3; end = $3 + 20; } print $1"\t"start"\t"end"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3 }' flanks.bed > 3prime.bed

# Generate 5' splice site output
awk '{ if ($6 == "+") { start = $3 - 3; end = $3 + 6; } else { start = $2 - 6; end = $2 + 3; } print $1"\t"start"\t"end"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3 }' flanks.bed > 5prime.bed

# Adjust using bedtools slop
bedtools slop -i 3prime.bed -g ~/data/ref_genome/v24/GRCh38.PRI.fa.fai -b 0 | sort | uniq > temp.bed && mv temp.bed 3prime.bed
bedtools slop -i 5prime.bed -g ~/data/ref_genome/v24/GRCh38.PRI.fa.fai -b 0 | sort | uniq >  temp.bed && mv temp.bed 5prime.bed

# Extract FASTA sequences
bedtools getfasta -fi ~/data/ref_genome/v24/GRCh38.PRI.fa -bed 3prime.bed -fo 3prime_sequences.fa
bedtools getfasta -fi ~/data/ref_genome/v24/GRCh38.PRI.fa -bed 5prime.bed -fo 5prime_sequences.fa
