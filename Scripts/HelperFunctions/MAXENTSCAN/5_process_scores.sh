# process results 
awk 'BEGIN { FS="\t"; OFS="\t" }
    /^>/ {
        key = substr($1, 2);
        getline;
        split($0, score_parts, "MAXENT:");
        score = score_parts[2];
        gsub(/[ \t]+/, "", score);
        scores[key] = score;
        next;
    }
    {
        key = $1 ":" $2 "-" $3;
        if (key in scores) {
            print $0, scores[key];
        }
    }' 3prime_scored.fa 3prime.bed > 3prime_scored.bed

# repeat for 5'scores
awk 'BEGIN { FS="\t"; OFS="\t" }
    /^>/ {
        key = substr($1, 2);
        getline;
        split($0, score_parts, "MAXENT:");
        score = score_parts[2];
        gsub(/[ \t]+/, "", score);
        scores[key] = score;
        next;
    }
    {
        key = $1 ":" $2 "-" $3;
        if (key in scores) {
            print $0, scores[key];
        }
    }' 5prime_scored.fa 5prime.bed > 5prime_scored.bed

# Combine 3' and 5' splice site scores into single file
awk 'BEGIN { FS="\t"; OFS="\t" }
    # Store information from 5prime_scored.bed using a key of chrom:start:end
    NR==FNR { key = $1 ":" $7 "-" $8; scores[key] = $9; next }

    # For 3prime_scored.bed, construct the same key, and print the matching values from both files
    { key = $1 ":" $7 "-" $8;
      if (key in scores) {
          print $1, $7, $8, $4, $5, $6, scores[key], $9
      }
    }' 5prime_scored.bed 3prime_scored.bed > epi_exons_splicesite_scores.bed
