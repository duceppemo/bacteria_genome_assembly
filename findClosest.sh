#!/bin/bash

#Find closest genome in refseq using Minhash (mash)


# I/O
output_folder="$1"
assembly="$2"
cpu="$3"

sample=$(basename "${assembly%.*}")

#resources
[ -z "$cpu" ] && cpu=$(nproc)

#Functions
function fixTable()
{
    table="$1"
    #Add a column with number of mash matching
    cat "$table" \
        | awk -F $'\t' 'BEGIN{OFS=FS} {split($5,array,"/"); $6=array[1]; print}' \
        > "${table}".tmp

    mv "${table}".tmp "$table"

    #sort table with the column added
    cat "$table" \
        | sort -rnk6,6 \
        > "${table}".tmp

    mv "${table}".tmp "$table"

    # #add header to table
    # echo -e "Reference-ID\tQuery-ID\tMash-distance\tP-value\tMatching-hashes\tmatch" > "${table}".tmp
    # cat "$table" >> "${table}".tmp

    # mv "${table}".tmp "$table"
}


############################################
#                                          #
#   Make a sketch from reference genomes   #
#                                          #
############################################


#create a list of genomes to compare
[ -e "${output_folder}"/"${sample}"_genomes.list ] && rm "${output_folder}"/"${sample}"_genomes.list #remove if exists because were are using append
for g in $(find "${output_folder}" -type f -name "*.fna.gz"); do
    echo "$g" >> "${output_folder}"/"${sample}"_genomes.list
done

#sketch reference genomes
mash sketch \
    -p "$cpu" \
    -l "${output_folder}"/"${sample}"_genomes.list \
    -o "${output_folder}"/"${sample}"_genomes.msh

#compute distance
# Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes:
mash dist \
    -p "$cpu" \
    "${output_folder}"/"${sample}"_genomes.msh \
    "$assembly" \
    > "${output_folder}"/"${sample}".distance.tsv

fixTable "${output_folder}"/"${sample}".distance.tsv

