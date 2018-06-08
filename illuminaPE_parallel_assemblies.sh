#!/bin/bash

version="0.2.3"


######################
#                    #
#    User Defined    #
#                    #
######################


#Analysis folder
export baseDir=""${HOME}"/analyses/salmonella_pigeon"

#reads
export reads="/media/6tb_raid10/data/salmonella_pigeon/merged"

#program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

# Centrifuge DB to use
export centrifuge_db="/media/6tb_raid10/db/centrifuge/2017-10-12_bact_vir_h"

#Kraken DB to use
export kraken_db="/media/6tb_raid10/db/kraken/refseq_BV_old"

#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=8

#Annotation
export kingdom="Bacteria"
export genus="Salmonella"
export species="enterica"
export gram="neg"
export locus_tag="TOCHANGE"
export centre="OLF"

# For assembly trimming
export smallest_contig=1000
export size=4900000


#################
#               #
#   Resources   #
#               #
#################


#computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"
export memCdHit=$((mem*1000))


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


#Folder structure
export fastq=""${baseDir}"/fastq"
export logs=""${baseDir}"/logs"
export qc=""${baseDir}"/QC"
export fastqc=""${qc}"/fastqc"
export kat=""${qc}"/kat"
export genomescope=""${qc}"/genomescope"
export blob=""${qc}"/blobtools"
export trimmed=""${baseDir}"/trimmed"
export merged=""${baseDir}"/merged"
export corrected=""${baseDir}"/corrected"
export assembly=""${baseDir}"/assembly"
export ordered=""${baseDir}"/ordered"
export annotation=""${baseDir}"/annotation"
export phaster=""${baseDir}"/phaster"
export amr=""${baseDir}"/amr"

#create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastqc" ] || mkdir -p "$fastqc"
[ -d "$kat" ] || mkdir -p "$kat"
[ -d "$genomescope" ] || mkdir -p "$genomescope"
[ -d "$blob" ] || mkdir -p "$blob"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$merged" ] || mkdir -p "$merged"
[ -d "$corrected" ] || mkdir -p "$corrected"
[ -d "$assembly" ] || mkdir -p "$assembly"
[ -d "$ordered" ] || mkdir -p "$ordered"
[ -d "$annotation" ] || mkdir -p "$annotation"
[ -d "$phaster" ] || mkdir -p "$phaster"
[ -d "$amr" ] || mkdir -p "$amr"


######################
#                    #
#   Initiating log   #
#                    #
######################


#Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#scipt version
echo -e "\n"$0" version "$version"\n" | tee -a "${logs}"/log.txt  # $0


#check if depenencies are installed
#if so, log version, if not exit

#bowtie2

# Centrifuge
if hash centrifuge 2>/dev/null; then  # if installed
    centrifuge --version | grep "centrifuge-class" | sed -e 's%^.*/%%' -e 's/-class//' | tee -a "${logs}"/log.txt
else
    echo >&2 "centrifuge was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Centrifuge database
if [ -s "${centrifuge_db}".1.cf ]; then
    echo "Centrifuge database: "$(basename "$centrifuge_db")"" | tee -a "${logs}"/log.txt
else
    echo "Could no find the provided Centrifude database. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#BBDuk
if hash bbduk.sh 2>/dev/null; then 
    bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbduk.sh was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#SPAdes
if hash spades.py 2>/dev/null; then
    spades.py -v | tee -a "${logs}"/log.txt
else
    echo >&2 "spades.py was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#BWA mem
if hash bwa 2>/dev/null; then
    version=$(bwa 2>&1 1>/dev/null | grep "Version" | sed 's/Version: //')
    echo "bwa v"$version"" | tee -a "${logs}"/log.txt
else
    echo >&2 "bwa was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi   

# Pilon
java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar
if [ $# -eq 0 ]; then
    java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar --version | tee -a "${logs}"/log.txt
else
    echo >&2 "pilon was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#QUAST
if hash spades.py 2>/dev/null; then
    quast.py -v | tee -a "${logs}"/log.txt
else
    echo >&2 "quast.py was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#CD-HIT-EST
if hash cd-hit-est 2>/dev/null; then
    cd-hit-est -h | head -n 1 | tr -d "=" | sed 's/^[ \t]*//;s/[ \t]*$//' | tee -a "${logs}"/log.txt
else
    echo >&2 "cd-hit-est was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Pilon
java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar
if [ $# -eq 0 ]; then
    java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar --version | tee -a "${logs}"/log.txt
else
    echo >&2 "pilon was not found. Aborting." | tee -a "${logs}"/log.txt
    #exit 1
fi

# Prokka
if hash prokka 2>/dev/null; then
    prokka --version 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "prokka was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Mash
if hash mash 2>/dev/null; then
    version=$(mash --version)
    echo "mash v"$version"" | tee -a "${logs}"/log.txt
else
    echo >&2 "mash was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#Mauve
echo "Mauve version \"Snapshot 2015-02-13\"" | tee -a "${logs}"/log.txt

# ProgessiveMauve
if hash "${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve 2>/dev/null; then
    "${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve --version 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "progressiveMauve was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# KAT
if hash kat 2>/dev/null; then
    version=$(kat --version | cut -d " " -f 2)
    echo "kat v"$version"" | tee -a "${logs}"/log.txt
else
    echo >&2 "kat was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#blobtools
if hash blobtools 2>/dev/null; then
    blobtools --version | tee -a "${logs}"/log.txt
else
    echo >&2 "blobtools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Jellyfish
if hash jellyfish 2>/dev/null; then
    jellyfish --version | tee -a "${logs}"/log.txt
else
    echo >&2 "jellyfish was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# GenomeScope
# it is a R script
# ????

# R
if hash R 2>/dev/null; then
    R --version | grep -F "R version" | tee -a "${logs}"/log.txt
else
    echo >&2 "R was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# ResFinder
perl "${prog}"/resfinder/resfinder.pl -h
if [ $# -eq 0 ]; then
    version=$(perl "${prog}"/resfinder/resfinder.pl -h | grep "Current" | cut -d ":" -f 2 | tr -d " ")
    echo "resfinder v"${version}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "resfinder was not found. Aborting." | tee -a "${logs}"/log.txt
    #exit 1
fi

# PHASTER
# Online phaster.ca

# Kraken
if hash kraken 2>/dev/null; then
    kraken --version | grep -F "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "kraken was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


#add space after prog version
echo -e "\n" | tee -a "${logs}"/log.txt


###################
#                 #
#   Fastq files   #
#                 #
###################


function compress_fastq ()
{
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')
    sample=$(basename "$r1" | cut -d '_' -f 1)

    clumpify.sh "$memJava" \
        in="$r1" \
        in2="$r2" \
        out="${fastq}"/"${sample}"/"${sample}"_R1.fastq.gz \
        out2="${fastq}"/"${sample}"/"${sample}"_R2.fastq.gz \
        reorder \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        2> >(tee "${logs}"/clump/"${sample}".txt)
}

export -f compress_fastq

[ -d "${logs}"/clump ] || mkdir -p "${logs}"/clump

find "$reads" -type f -name "*_R1*fastq.gz" \
    | parallel  --bar \
                --env compress_fastq \
                --env fastq \
                --env memJava \
                --env logs \
                --env cpu \
                --env maxProc \
                --jobs "$maxProc" \
                'compress_fastq {}'


###############
#             #
#   Read QC   #
#             #
###############


### FastQC ###

function run_fastqc()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))

    [ -d "${2}"/"$sample" ] || mkdir -p "${2}"/"$sample"

    fastqc \
        --o "${2}"/"$sample" \
        --noextract \
        --threads $((cpu/maxProc)) \
        "$r1" "$r2"
}

# Create folder to store report
[ -d "${qc}"/fastqc/raw ] || mkdir -p "${qc}"/fastqc/raw

export -f run_fastqc

find "$fastq" -type f -name "*_R1*fastq.gz" \
    | parallel  --bar \
                --env run_fastqc \
                --env qc \
                --env cpu \
                --env maxProc \
                --jobs "$maxProc" \
                "run_fastqc {} "${qc}"/fastqc/raw"

#Merge all FastQC reports together
multiqc \
    -o "${qc}"/fastqc/raw \
    -n merged_reports.html \
    "${qc}"/fastqc/raw


### Read taxonomic assignement ###

# Run Centrifuge is only one sample?
# Run Centrifuge is want to use nt database (about 40min per sample)
# Run Kraken if multiple samples and small database (e.g. Bacteria+viruses)


### Centrifuge ###

function run_centrifuge()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))

    [ -d "${2}"/"${sample}" ] || mkdir -p "${2}"/"${sample}"

    #build the command
    centrifuge \
        -p "$cpu" \
        -t \
        --seed "$RANDOM" \
        -x "$centrifuge_db" \
        -1 "$r1" \
        -2 "$r2" \
        --report-file "${2}"/"${sample}"/"${sample}"_report.tsv \
        > "${2}"/"${sample}"/"${sample}".tsv

    cat "${2}"/"${sample}"/"${sample}".tsv | \
        cut -f 1,3 | \
        ktImportTaxonomy /dev/stdin -o "${2}"/"${sample}"/"${sample}".html
}

for i in $(find "$fastq" -type f -name "*_R1*fastq.gz"); do
    run_centrifuge "$i" "${qc}"/centrifuge/raw
done


### Kraken ###

function run_kraken ()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))

    # Create folder to store report
    [ -d "${2}"/"${sample}" ] || mkdir -p "${2}"/"${sample}"

    #run Kraken
    kraken \
        --db "$kraken_db" \
        --output "${2}"/"${sample}"/"${sample}".kraken \
        --threads "$cpu" \
        --gzip-compressed \
        --check-names \
        --fastq-input \
        --paired "$r1" "$r2" \
        &> >(tee "${2}"/"${sample}"/"${sample}".kraken.log)

    #Create Kraken report
    kraken-report \
        --db "$kraken_db" \
        "${2}"/"${sample}"/"${sample}".kraken \
        > "${2}"/"${sample}"/"${sample}".report.tsv

    #Prepare result for display with Krona
    cat "${2}"/"${sample}"/"${sample}".kraken | \
        cut -f 2-3 | \
        ktImportTaxonomy /dev/stdin -o "${2}"/"${sample}"/"${sample}".html

}

# Put database into memory (cheat)
cat "${kraken_db}"/database.kdb > /dev/null  # About 6min to run

for i in $(find "$fastq" -type f -name "*_R1*fastq.gz"); do
    run_kraken "$i" "${qc}"/kraken/raw
done

### KAT ###

function run_kat()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))
    
    # GC plot on reads
    kat gcp \
        -t $((cpu/maxProc)) \
        -o "${2}"/"${sample}"/"$sample"_gcp \
        "$r1" "$r2"

    # Redo gpc plot to add labels to axes
    kat plot density \
        -o "${2}"/"${sample}"/"$sample"_gcp.mx.png \
        --y_label "GC count" \
        --x_label "K-mer multiplicity (R1 and R2 combined)" \
        --z_label "Distinct K-mers per bin" \
        "${2}"/"${sample}"/"$sample"_gcp.mx

    # k-mer frequency
    kat hist \
        -t $((cpu/maxProc)) \
        -o "${2}"/"${sample}"/"$sample"_histo \
        "$r1" "$r2"

    # Redo histo to add labels to axes
    kat plot spectra-hist \
        -o "${2}"/"${sample}"/"$sample"_histo.png \
        --y_label "Count" \
        --x_label "K-mer multiplicity" \
        "${2}"/"${sample}"/"$sample"_histo

    # Compare R1 vs R2
    kat comp \
        -t $((cpu/maxProc)) \
        -o "${2}"/"${sample}"/"$sample" \
        "$r1" "$r2"

    # density plot R1 vs R2
    kat plot density \
        -o "${2}"/"${sample}"/"${sample}"_r1_vs_r2.png \
        --y_label "K-mer multiplicity for $(basename "${r2%.*}")" \
        --x_label "K-mer multiplicity for $(basename "${r1%.*}")" \
        --z_label "Distinct K-mers per bin" \
        "${2}"/"${sample}"/"${sample}"-main.mx

    kat plot spectra-mx \
        -o "${2}"/"${sample}"/"${sample}"_spectra-mx.png \
        --intersection \
        --y_label "27-mer multiplicity for $(basename "${r2%.*}")" \
        --x_label "27-mer multiplicity for $(basename "${r1%.*}")" \
        "${2}"/"${sample}"/"${sample}"-main.mx
}

#make function available to parallel
export -f run_kat

[ -d "${kat}"/raw ] || mkdir -p "${kat}"/raw

source activate kat
# run samples in parallel
find "$fastq" -type f -name "*R1*.fastq.gz" \
    | parallel  --bar \
                --env run_kat \
                --env cpu \
                --env kat \
                --env maxProc \
                --jobs "$maxProc" \
                "run_kat {} "${kat}"/raw"

source deactivate

# Cleanup
find "$kat" -type f \
    | grep -vE ".png|.stats" \
    | xargs rm -r


### GenomeScope ###

function genomeStats()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))

    [ -d "${2}"/"$sample" ] || mkdir -p "${2}"/"$sample"

    zcat "$r1" "$r2"\
        | jellyfish count \
            -s "$mem" \
            -m 21 \
            -C \
            -t $((cpu/maxProc)) \
            -o "${2}"/"${sample}"/"${sample}".jf \
            -L 3 \
            /dev/stdin

    jellyfish histo \
        -t $((cpu/maxProc)) \
        -o "${2}"/"${sample}"/"${sample}".histo \
        "${2}"/"${sample}"/"${sample}".jf

    read_length=$(zcat "$1" \
        | head -n 4 \
        | sed -n '2p' \
        | tr -d "\n" \
        | wc -m)

    Rscript "${prog}"/genomescope/genomescope.R \
        "${2}"/"${sample}"/"${sample}".histo \
        21 \
        "$read_length" \
        "${2}"/"${sample}"
}

export -f genomeStats

[ -d "${genomescope}"/raw ] || mkdir -p "${genomescope}"/raw

find "$fastq" -type f -name "*_R1*" |
    parallel    --bar \
                --env genomeStats \
                --env output \
                --env mem \
                --env cpu \
                --env maxProc \
                --env genomescope \
                --jobs "$maxProc" \
                "genomeStats {} "${genomescope}"/raw"

# Make coverage report
echo -e "Sample\tEst_Coverage(X)\tEst_Genome_Size(bp)" > "${genomescope}"/estimated_coverages.tsv

for i in $(find "$genomescope" -type f -name "model.txt"); do
    name="$(dirname "$i")"
    sample=$(basename "$name")

    coverage=$(cat "$i" \
                | grep -E "^kmercov" \
                | tr -s " " \
                | cut -d " " -f 2)
    length=$(cat "$i" \
                | grep -E "^length" \
                | tr -s " " \
                | cut -d " " -f 2)
    cov=$(printf '%.0f' "$coverage")
    len=$(printf '%.0f' "$length")

    echo -e ""$sample"\t"$cov"\t"$len"" >> "${genomescope}"/estimated_coverages.tsv
done


###########################
#                         #
#   Read Pre-Processing   #
#                         #
###########################


#Create folders for logs for each steps
[ -d "${logs}"/tile_filtering ] || mkdir -p "${logs}"/tile_filtering
[ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming
[ -d "${logs}"/contaminant_removal ] || mkdir -p "${logs}"/contaminant_removal
[ -d "${logs}"/correction_step1 ] || mkdir -p "${logs}"/correction_step1
[ -d "${logs}"/correction_step2 ] || mkdir -p "${logs}"/correction_step2
[ -d "${logs}"/correction_step3 ] || mkdir -p "${logs}"/correction_step3
[ -d "${logs}"/merging ] || mkdir -p "${logs}"/merging

function preprocess_reads ()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))

    #Removing low-quality regions
    filterbytile.sh "$memJava" \
        in="$r1" \
        in2="$r2" \
        out="${trimmed}"/"${sample}"/"${sample}"_Filtered_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"/"${sample}"_Filtered_2P.fastq.gz \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        2> >(tee "${logs}"/tile_filtering/"${sample}".txt)

    rm "$r1" "$r2"

    # #QC
    # [ -d "${qc}"/fastqc/tile ] || mkdir -p "${qc}"/fastqc/tile
    # run_fastqc "$1" "${qc}"/fastqc/til

    #Quality trimming and adapter trimming
    bbduk.sh "$memJava" \
        in="${trimmed}"/"${sample}"/"${sample}"_Filtered_1P.fastq.gz \
        in2="${trimmed}"/"${sample}"/"${sample}"_Filtered_2P.fastq.gz \
        ref="${prog}"/bbmap/resources/adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \
        qtrim=lr trimq=10 \
        minlen=64 \
        out="${trimmed}"/"${sample}"/"${sample}"_Trimmed_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"/"${sample}"_Trimmed_2P.fastq.gz \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        ordered=t \
        2> >(tee "${logs}"/trimming/"${sample}".txt)

    rm "${trimmed}"/"${sample}"/"${sample}"_Filtered_?P.fastq.gz

    # #QC
    # [ -d "${qc}"/fastqc/trimmed ] || mkdir -p "${qc}"/fastqc/trimmed
    # run_fastqc \
    #     "${trimmed}"/"${sample}"/"${sample}"_Trimmed_1P.fastq.gz \
    #     "${qc}"/fastqc/trimmed
        

    #Removing synthetic artifacts and spike-ins
    #Add to "ref" any contaminants found with centrifuge
    bbduk.sh "$memJava" \
        in="${trimmed}"/"${sample}"/"${sample}"_Trimmed_1P.fastq.gz \
        in2="${trimmed}"/"${sample}"/"${sample}"_Trimmed_2P.fastq.gz \
        ref="${prog}"/bbmap/resources/phix174_ill.ref.fa.gz \
        k=31 \
        out="${trimmed}"/"${sample}"/"${sample}"_Cleaned_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"/"${sample}"_Cleaned_2P.fastq.gz \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        ordered=t \
        2> >(tee "${logs}"/contaminant_removal/"${sample}".txt)

    rm "${trimmed}"/"${sample}"/"${sample}"_Trimmed_?P.fastq.gz

    # #QC
    # [ -d "${qc}"/fastqc/cleaned ] || mkdir -p "${qc}"/fastqc/cleaned
    # run_fastqc \
    #     "${trimmed}"/"${sample}"/"${sample}"_Cleaned_1P.fastq.gz \
    #     "${qc}"/fastqc/cleaned
        

    #Correcting Illumina paired-end reads
    #Phase 1
    bbmerge.sh "$memJava" \
        in="${trimmed}"/"${sample}"/"${sample}"_Cleaned_1P.fastq.gz \
        in2="${trimmed}"/"${sample}"/"${sample}"_Cleaned_2P.fastq.gz \
        out="${corrected}"/"${sample}"/"${sample}"_Cor1_1P.fastq.gz \
        out2="${corrected}"/"${sample}"/"${sample}"_Cor1_2P.fastq.gz \
        ecco=t \
        mix=t \
        verystrict=t \
        ordered=t \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        ordered ihist="${logs}"/correction_step1/"${sample}"_ihist_corr_merge.txt \
        2> >(tee "${logs}"/correction_step1/"${sample}".txt)

    rm "${trimmed}"/"${sample}"/"${sample}"_Cleaned_?P.fastq.gz

    # Create folder to store insert size results
    [ -d "${qc}"/insert_size/"${sample}" ] || mkdir -p "${qc}"/insert_size/"${sample}"

    # Insert size
    cat "${logs}"/correction_step1/ihist_corr_merge.txt \
        | grep -vE "^#" \
        > "${qc}"/insert_size/"${sample}"/"${sample}"_insert_bbtools.tsv

    # Plot
    gnuplot -e "reset; \
        set term png; \
        set output '"${qc}"/insert_size/"${sample}"/"${sample}"_insert_bbtools.png'; \
        set xlabel 'Insert Size (bp)'; \
        set ylabel 'Read Count'; \
        set title 'Paired-end reads insert size distribution (BBmerge)'; \
        stats '"${qc}"/insert_size/"${sample}"/"${sample}"_insert_bbtools.tsv' u 1:2 nooutput; \
        set xrange [STATS_min_x:STATS_max_x]; \
        set yrange [STATS_min_y:STATS_max_y]; \
        set key inside top left; \
        plot '"${qc}"/insert_size/"${sample}"/"${sample}"_insert_bbtools.tsv' using 1:2 title \""$sample"\" with linespoints; \
        set output"

    rm "${qc}"/insert_size/"${sample}"/"${sample}"_insert_bbtools.tsv

    #Phase2
    clumpify.sh "$memJava" \
        in="${corrected}"/"${sample}"/"${sample}"_Cor1_1P.fastq.gz \
        in2="${corrected}"/"${sample}"/"${sample}"_Cor1_2P.fastq.gz \
        out="${corrected}"/"${sample}"/"${sample}"_Cor2_1P.fastq.gz \
        out2="${corrected}"/"${sample}"/"${sample}"_Cor2_2P.fastq.gz \
        ecc=t \
        passes=4 \
        reorder=t \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        2> >(tee "${logs}"/correction_step2/"${sample}".txt)

    rm "${corrected}"/"${sample}"/"${sample}"_Cor1_?P.fastq.gz

    #Phase3
    tadpole.sh "$memJava" \
        in="${corrected}"/"${sample}"/"${sample}"_Cor2_1P.fastq.gz \
        in2="${corrected}"/"${sample}"/"${sample}"_Cor2_2P.fastq.gz \
        out="${corrected}"/"${sample}"/"${sample}"_Cor3_1P.fastq.gz \
        out2="${corrected}"/"${sample}"/"${sample}"_Cor3_2P.fastq.gz \
        ecc=t \
        k=62 \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        mode=correct \
        2> >(tee "${logs}"/correction_step2/"${sample}".txt)

    rm "${corrected}"/"${sample}"/"${sample}"_Cor2_?P.fastq.gz

    #Merging Illumina overlapping paried-end reads
    bbmerge.sh "$memJava" \
        in="${corrected}"/"${sample}"/"${sample}"_Cor3_1P.fastq.gz \
        in2="${corrected}"/"${sample}"/"${sample}"_Cor3_2P.fastq.gz \
        out="${merged}"/"${sample}"/"${sample}"_merged.fastq.gz \
        outu="${merged}"/"${sample}"/"${sample}"_unmerged_1P.fastq.gz \
        outu2="${merged}"/"${sample}"/"${sample}"_unmerged_2P.fastq.gz \
        strict=t \
        k=93 \
        extend2=80 \
        rem=t \
        ordered=t \
        threads=$((cpu/maxProc)) \
        ziplevel=9 \
        ihist="${logs}"/merging/"${sample}"_ihist_merge.txt \
        2> >(tee "${logs}"/merging/"${sample}".txt)

    # rm "${corrected}"/"${sample}"/"${sample}"_Cor3_?P.fastq.gz

    # #QC
    [ -d "${qc}"/fastqc/merged/"$sample" ] || mkdir -p "${qc}"/fastqc/merged/"$sample"
    fastqc \
        --o "${qc}"/fastqc/merged/"$sample" \
        --noextract \
        --threads $((cpu/maxProc)) \
        "${merged}"/"${sample}"/"${sample}"_merged.fastq.gz \
        "${merged}"/"${sample}"/"${sample}"_unmerged_1P.fastq.gz \
        "${merged}"/"${sample}"/"${sample}"_unmerged_2P.fastq.gz
}

export -f preprocess_reads

find "$fastq" -type f -name "*_R1*" |
    parallel    --bar \
                --env preprocess_reads \
                --env memJava \
                --env cpu \
                --env maxProc \
                --env trimmed \
                --env logs \
                --env prog \
                --env corrected \
                --env merged \
                --jobs "$maxProc" \
                'preprocess_reads {}'

#Merge all FastQC reports together
multiqc \
    -o "${qc}"/fastqc/merged \
    -n merged_reports.html \
    "${qc}"/fastqc/merged

rm -rf "$fastq" "$trimmed"


################
#              #
#   Assembly   #
#              #
################


function assemble ()
{
    m="$1"
    mu1="${m%_merged.fastq.gz}"_unmerged_1P.fastq.gz
    mu2="${m%_merged.fastq.gz}"_unmerged_2P.fastq.gz
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${assembly}"/"$sample" ] || mkdir -p "${assembly}"/"$sample"
    
    #use merged paired-end as single-end reads
    python3 "${prog}"/Unicycler/unicycler-runner.py \
        -1 "${m%_merged.fastq.gz}"_unmerged_1P.fastq.gz \
        -2 "${m%_merged.fastq.gz}"_unmerged_2P.fastq.gz \
        -s "$1" \
        -o "${assembly}"/"$sample" \
        -t $((cpu/maxProc)) \
        --keep 3 \
        --no_correct \
        --verbosity 2 \
        --mode normal \
        --pilon_path "${prog}"/pilon/pilon-dev.jar

    cp "${assembly}"/"${sample}"/assembly.fasta \
        "${assembly}"/"${sample}"/"${sample}".fasta

    find "${assembly}"/"${sample}"/pilon_polish ! -name "*.out" ! -name "*.changes" -exec rm -rf {} \;
}

export -f assemble

find "$merged" -type f -name "*_merged.fastq.gz" | \
    parallel    --bar \
                --env assemble \
                --env merged \
                --env cpu \
                --env maxProc \
                --env assembly \
                --env prog \
                --jobs "$maxProc" \
                "assemble {}"

function trimAssembly()
{
    path=$(dirname "$1")
    sample=$(basename "$path")

    # Remove contigs smaller than 1000bp
    perl "${prog}"/phage_typing/removesmallscontigs.pl \
        "${smallest_contig}" \
        "$1" \
        > "${assembly}"/"${sample}"/"${sample}"_trimmed"${smallest_contig}".fasta
}

export -f trimAssembly

find "${assembly}" -type f -maxdepth 2 -name "*.fasta" \
    | parallel  --bar \
                --env prog \
                --env smallest_contig \
                --env assembly \
                'trimAssembly {}'


#######################
#                     #
#   Contig ordering   #
#                     #
#######################


# Folder to store refseq genomes
[ -d "${ordered}"/refseq ] || mkdir -p "${ordered}"/refseq
[ -d "${qc}"/distance ] || mkdir -p "${qc}"/distance

# Download all assemblies from species
bash "${scripts}"/get_assemblies.sh \
    -q ""$genus" "$species"" \
    -t 'refseq' \
    -a 'Complete Genome' \
    -o ""${ordered}"/refseq"

function order_contigs ()
{
    sample=$(cut -d "." -f 1 <<< $(basename "$1"))

    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ];then
        # Use mash to find closest genome from refseq
        bash "${scripts}"/findClosest.sh \
            "${ordered}"/refseq \
            "$1" \
            $((cpu/maxProc))

        distance=""${ordered}"/refseq/"${sample}".distance.tsv"
        closest=$(cat "$distance" | head -n 1 | cut -f 1)  # The closest is the top one
        closest_acc=$(zcat "$closest" | head -n 1 | cut -d " " -f 1 | tr -d ">")
        closest_ID=$(zcat "$closest" | head -n 1 | cut -d " " -f 2- | cut -d "," -f 1)
        score=$(cat "$distance" | head -n 1 | cut -f 6)  #its score out of 1000
        acc=$(cut -d "_" -f 1,2 <<< $(basename "$closest"))  # accession number of the closest match

        echo ""$sample" closest genome is \""${closest_ID}" ("${closest_acc}")\" with a score of "${score}"/1000" | tee -a "${logs}"/log.txt  #Add information to log

        #move distance file
        mv "${ordered}"/refseq/"${sample}".distance.tsv "${qc}"/distance

        # only order contigs if closest genome is at least 80% similar
        if [ "$score" -ge 800 ]; then
            # uncompress downloaded fasta file for Mauve
            [ -s "${closest%.gz}" ] || pigz -p $((cpu/maxProc)) -d -k "$closest" # decompress if not present

            # Use Mauve in batch mode to order contigs with closest genome
            java "$memJava" -cp "${prog}"/mauve_snapshot_2015-02-13/Mauve.jar \
                org.gel.mauve.contigs.ContigOrderer \
                -output "${ordered}"/mauve/"$sample" \
                -ref "${closest%.gz}" \
                -draft "$1"

            #fix formating of Mauve output
            #find out how many "alignment" folder there is
            n=$(find "${ordered}"/mauve/"$sample" -maxdepth 1 -type d | sed '1d' | wc -l)

            # reformat with no sorting
            perl "${scripts}"/formatFasta.pl \
                -i "${ordered}"/mauve/"${sample}"/alignment"${n}"/"$(basename "$1")" \
                -o "${ordered}"/"${sample}"_ordered.fasta \
                -w 80

            #align with progessiveMauve
            [ -d "${qc}"/mauve/"${sample}" ] || mkdir -p "${qc}"/mauve/"${sample}"

            "${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve \
                --output="${qc}"/mauve/"${sample}"/"${sample}"_ordered.xmfa \
                "${closest%.gz}" \
                "${ordered}"/"${sample}"_ordered.fasta

            # View alignemnt with mauve
            # "${prog}"/mauve_snapshot_2015-02-13/./Mauve \
            #     "${qc}"/mauve/"${sample}"_ordered.xmfa &
         else
            perl "${scripts}"/formatFasta.pl \
                -i "$1" \
                -o "${ordered}"/"${sample}"_ordered.fasta \
                -w 80
        fi
    else
        circlator fixstart \
            --verbose \
            "$1" \
            "${ordered}"/"${sample}"_ordered
    fi
}

export -f order_contigs

find "$assembly" -type f -name "*trimmed"${smallest_contig}"*" | \
    parallel    --bar \
                --env order_contigs \
                --env cpu \
                --env maxProc \
                --env prog \
                --env scripts \
                --env ordered \
                --env memJava \
                --jobs "$maxProc" \
                "order_contigs {}"

#cleanup
rm -rf "${ordered}"/mauve
rm -rf "${ordered}"/refseq
find "$ordered" -type f -name "*.sslist" -exec rm {} \;


###################
#                 #
#   Assembly QC   #
#                 #
###################


### blast ###

#TODO -> Make output reusable by blobtools
# keep more output
# create a filtered version to replicate current with only best hit per contig

function blast()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "${1%.fasta}"))

    # blobtools has the following requirements for the blast output
        # 1st column: sequenceID (must be part of the assembly)
        # 2nd column: TaxID (a NCBI TaxID)
        # 3rd column: score (a numerical score)
    blastn \
        -task megablast \
        -db nt \
        -query "$1" \
        -out "${qc}"/blast/"${sample}".all.blastn.tsv \
        -evalue "1e-25" \
        -outfmt '6 qseqid staxids bitscore sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue sscinames sskingdoms' \
        -num_threads $((cpu/maxProc)) \
        -culling_limit 5

    # Best hit only
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsscinames\tsskingdoms\tstaxids" \
        > "${qc}"/blast/"${sample}".blastn.tsv.tmp

    cat "${qc}"/blast/"${sample}".all.blastn.tsv \
        | sort -k1,1n -k3,3gr \
        | sort -uk1,1n \
        | awk -F $'\t' 'BEGIN {OFS = FS} {print $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $3, $15, $16, $2}' \
        >> "${qc}"/blast/"${sample}".blastn.tsv.tmp

    mv "${qc}"/blast/"${sample}".blastn.tsv.tmp \
        "${qc}"/blast/"${sample}".bestHit.blastn.tsv

    # Formated for Blobtools
    [ -d "${blob}"/"$sample" ] || mkdir -p "${blob}"/"$sample"
    cat "${qc}"/blast/"${sample}".all.blastn.tsv \
        | awk -F $'\t' 'BEGIN {OFS = FS} {print $1, $2, $3}' \
        > "${blob}"/"${sample}"/"${sample}".blast_nt.tsv
}

export -f blast

[ -d "${qc}"/blast ] || mkdir -p "${qc}"/blast

find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel  --bar \
                --env blast \
                --env cpu \
                --env maxProc \
                --env qc \
                --env blob \
                --jobs "$maxProc" \
                "blast {}"

#clean blast index files
find "$ordered" -type f ! -name "*.fasta" -exec rm {} \;


########
########  Have a look at the IDs and remove contaminant contigs from assembly
########


# ### Plasmid detection ###

# # Detect contigs part of a plasmid
# [ -d "${qc}"/plasmid ] || mkdir -p "${qc}"/plasmid

# function detect_plasmid ()
# {
#     # should be more elaborate than this...
#     # should look at the blast result to identify the plasmid-related contigs
#     # should look at the assembly graph to see which contigs belong to the same component (or plasmid)
#     # Should make a nice report.

#     sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

#     if [[ -n $(cat "${qc}"/blast/$(basename "${1%.fasta}").blastn.tsv | grep -i "plasmid") ]]; then
#         echo "Plasmid(s) detected in "$sample" assembly" \
#             | tee -a "${logs}"/log.txt

#         # Create report for plasmid detection
#         echo -e "Contig(s) from "$sample" assembly identified as plasmid:\n" \
#             | tee -a "${qc}"/plasmid/"${sample}".txt
        
#         cat "${qc}"/blast/$(basename "${1%.fasta}").blastn.tsv \
#             | grep -i "plasmid" \
#             | cut -f 1,3 \
#             | sort -u -k1,1 \
#             | tee -a "${qc}"/plasmid/"${sample}".txt
#     else
#         echo "No plasmid detected in "$sample" assembly" \
#             | tee -a "${logs}"/log.txt
#     fi
# }

# export -f detect_plasmid

# find "$ordered" -type f -name "*_ordered.fasta" \
#     | parallel  --bar \
#                 --env detect_plasmid \
#                 --env qc \
#                 --jobs "$maxProc" \
#                 'detect_plasmid {}'


### coverage ###

function get_coverage()  # unsing unmerged reads only
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    [ -d "${qc}"/coverage/"$sample" ] || mkdir -p "${qc}"/coverage/"$sample"

    bwa index "$1"

    #Align corrected paired-end
    rg_pe="@RG\tID:"${sample}"\tCN:"${centre}"\tLB:NexteraXT\tPL:ILLUMINA\tSM:"${sample}""
    r1="${corrected}"/"${sample}"/"${sample}"_Cor3_1P.fastq.gz
    r2="${corrected}"/"${sample}"/"${sample}"_Cor3_2P.fastq.gz

    bwa mem -t $((cpu/maxProc)) -M -R "$rg_pe" "$1" "$r1" "$r2" | \
    samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
    samtools sort -@ $((cpu/maxProc)) -m 10G - | \
    samtools rmdup - "${qc}"/coverage/"${sample}"/"${sample}".bam

    samtools index "${qc}"/coverage/"${sample}"/"${sample}".bam

    #Average genome depth of coverage
    average_cov=$(samtools depth \
        "${qc}"/coverage/"${sample}"/"${sample}".bam  \
        | awk '{sum+=$3} END { print sum/NR}')

    printf "%s\t%.*f\n" "$sample" 0 "$average_cov" | tee -a "${qc}"/coverage/average_cov.tsv
}

export -f get_coverage

[ -d "${qc}"/coverage ] || mkdir -p "${qc}"/coverage
echo -e "Sample\tAverage_Cov" > "${qc}"/coverage/average_cov.tsv

find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" | \
    parallel    --bar \
                --env get_coverage \
                --env cpu \
                --env maxProc \
                --env qc \
                --env corrected \
                --env centre \
                --jobs "$maxProc"  \
                "get_coverage {}"

#clean bwa index files
find "$ordered" -type f ! -name "*.fasta" -exec rm {} \;


### quast ###

# #Merge all bam files
# declare -a bams=()
# for i in $(find "${qc}"/coverage -type f -name "*.bam"); do 
#     bams+=("$i")
# done

# samtools merge -@ "$cpu" - ${bams[@]} | \
# samtools rmdup - - | \
# samtools sort -@ "$cpu" -m 10G -o "${qc}"/coverage/all.bam -
# samtools index "${qc}"/coverage/all.bam
#All the genomes compared
declare -a genomes=()
for i in $(find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta"); do 
    genomes+=("$i")
done

source activate quast

quast.py \
    --output-dir "${qc}"/quast/all \
    --threads "$cpu" \
    --min-contig "$smallest_contig" \
    --est-ref-size "$size" \
    ${genomes[@]}

# Make quast report on individual assembly
function run_quast()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    quast.py \
        -m "$smallest_contig" \
        -t $((cpu/maxProc)) \
        -o "${qc}"/quast/"$sample" \
        --min-contig "$smallest_contig" \
        --est-ref-size "$size" \
        $1  # don't put in quotes
}

#make function available to parallel
export -f run_quast  # -f is to export functions

#run paired-end merging on multiple samples in parallel
find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel  --bar \
                --env run_quast \
                --env maxProc \
                --env cpu \
                --env qc \
                --env smallest_contig \
                --jobs "$maxProc" \
                "run_quast {}"

source deactivate

### kat ###

# Create output folder
[ -d "${kat}"/assembly ] || mkdir -p "${kat}"/assembly

# unzipping the file is no longer required in Kat v2.4.1 

# genome assembly analysis
# both paired-end compared as one group
function katAssembly ()
{
    sample=$(basename "$1" | cut -d "_" -f 1)

    kat comp \
        -t $((cpu/maxProc)) \
        -o "${kat}"/assembly/"${sample}"/"${sample}"_pe_vs_asm \
        ""${corrected}"/"${sample}"/"${sample}"_Cor3_?P.fastq.gz" \
        "$1"

    # Redo the plot to rename the axes
    kat plot spectra-cn \
        -o "${kat}"/assembly/"${sample}"/"${sample}"_pe_vs_asm-main.mx.spectra-cn.png \
        --y_label "Number of distinct K-mers" \
        --x_label "K-mer multiplicity" \
        "${kat}"/assembly/"${sample}"/"${sample}"_pe_vs_asm-main.mx

    # kat plot spectra-mx \
    #     -o "${kat}"/assembly/"${sample}"/"${sample}"_spectra-mx.png \
    #     --y_label "# Distinct Kmers" \
    #     --x_label "Kmer Frequency" \
    #     --list c1,r1 \
    #     "${kat}"/assembly/"${sample}"/"${sample}"_pe_vs_asm-main.mx

    kat sect  \
        -t $((cpu/maxProc)) \
        -o "${kat}"/assembly/"${sample}"/"${sample}"_gpc \
        "$1" \
        "${corrected}"/"${sample}"/"${sample}"_Cor3_1P.fastq.gz \
        "${corrected}"/"${sample}"/"${sample}"_Cor3_2P.fastq.gz

    # Plot coverage for individual contigs
    contigs=($(cat "$1" | grep -F ">" | cut -d " " -f 1 | tr -d ">" | tr "\n" " "))
    # echo "${contigs[@]}"

    [ -d "${kat}"/assembly/"${sample}"/coverage_plots ] || mkdir -p "${kat}"/assembly/"${sample}"/coverage_plots

    for i in "${contigs[@]}"; do
        kat plot profile \
            -o "${kat}"/assembly/"${sample}"/coverage_plots/"${sample}"_"${i}".png \
            --y_label "Coverage" \
            --x_label "Size (bp)" \
            --title "Coverage of contig "${i}"" \
            --index "$i" \
            "${kat}"/assembly/"${sample}"/"${sample}"_gpc-counts.cvg
    done

    kat cold \
        -t $((cpu/maxProc)) \
        -o "${kat}"/assembly/"${sample}"/"${sample}"_cold \
        "$1" \
        "${corrected}"/"${sample}"/"${sample}"_Cor3_1P.fastq.gz \
        "${corrected}"/"${sample}"/"${sample}"_Cor3_2P.fastq.gz

    python3 "${CONDA_PREFIX}"/lib/python3.6/local/kat/distanalysis.py \
        -o "${kat}"/assembly/"${sample}"/"${sample}" \
        "${kat}"/assembly/"${sample}"/"${sample}"_pe_vs_asm-main.mx \
        --plot \
        --format 'png' \
        --verbose \
        | tee "${kat}"/assembly/"${sample}"/"${sample}"_decompostion_analysis.stats
}

export -f katAssembly

source activate kat

find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel  --bar \
                --env katAssembly \
                --env cpu \
                --env maxProc \
                --env kat \
                --jobs "$maxProc" \
                'katAssembly {}'

source deactivate

# Cleanup
find "${kat}"/assembly -type f \
    | grep -vE ".png|stats" \
    | xargs rm -f


### blobtools ###

# Activate blobtools virtual environment
source "${HOME}"/my_virtualenv/blob_env/bin/activate

function run_blobtools ()
{
    sample=$(basename "$1" | cut -d "_" -f 1)
    genome="$1"

    [ -d "${blob}"/"$sample" ] || mkdir -p "${blob}"/"$sample"

    # Make coverage file from bam
    blobtools map2cov \
        -i "$genome" \
        -b "${qc}"/coverage/"${sample}"/"${sample}".bam \
        -o "${blob}"/"${sample}"/

    #blast assembly on nt
    # 1st column: sequenceID (must be part of the assembly)
    # 2nd column: TaxID (a NCBI TaxID)
    # 3rd column: score (a numerical score)
    # blastn \
    #     -task megablast \
    #     -query "$genome" \
    #     -db nt \
    #     -out "${blob}"/"${sample}"/"${sample}".blast_nt.tsv \
    #     -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    #     -culling_limit 5 \
    #     -evalue 1e-25 \
    #     -num_threads $((cpu/maxProc))

    # Create BlobDB
    blobtools create \
        -i "$genome" \
        -c "${blob}"/"${sample}"/"${sample}".bam.cov \
        -t "${blob}"/"${sample}"/"${sample}".blast_nt.tsv \
        -o "${blob}"/"${sample}"/"${sample}"

    # Print BlobDB as a table
    blobtools view \
        -i "${blob}"/"${sample}"/"${sample}".blobDB.json \
        -o "${blob}"/"${sample}"/

    # Colour a covplot by GC-categories
    # Extract sequenceID and GC
    grep -v '^#' "${blob}"/"${sample}"/"${sample}".blobDB.table.txt \
        | cut -f 1,3 \
        > "${blob}"/"${sample}"/"${sample}".blobDB.id.gc.txt

    # Divide into GC categories 
    i=0
    division=5
    while [ "$i" -lt $((100+division)) ]; do
        # echo "$i"
        cat "${blob}"/"${sample}"/"${sample}".blobDB.id.gc.txt \
            | awk -v count="$i" -v div="$division" '
            BEGIN {OFS = ","}
            {
                if (($2 * 100) >= count && ($2 * 100) < (count + div))
                {
                    $2 = count"-"(count + (div - 1))"%"
                    print
                }
            }' >> "${blob}"/"${sample}"/"${sample}".blobDB.id.gc.catcolour.txt
        let i=$i+division
    done

    # Generate covplot with colour by %GC categories
    blobtools covplot \
        -i "${blob}"/"${sample}"/"${sample}".blobDB.json \
        -c "${blob}"/"${sample}"/"${sample}".bam.cov \
        --catcolour "${blob}"/"${sample}"/"${sample}".blobDB.id.gc.catcolour.txt \
        -o "${blob}"/"${sample}"/

    # Plot BlobDB as a blobplot
    blobtools blobplot \
      -i "${blob}"/"${sample}"/"${sample}".blobDB.json \
      -o "${blob}"/"${sample}"/

    #cleanup (only keep the image files)
    find "${blob}"/"$sample" -type f | grep -vF ".png" | xargs rm  # blobtools temp files
    # find "${ordered}"/"$sample" -type f -name "*.fasta.*" -exec rm {} +  # blast index files
}

export -f run_blobtools

find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel  --bar \
                --env run_blobtools \
                --env cpu \
                --env maxProc \
                --env blob \
                --env qc \
                --env corrected \
                --env centre \
                --jobs "$maxProc" \
                'run_blobtools {}'

# Remove bam files
find "${qc}"/coverage -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} \;

# Deactivate the virtual environment
deactivate


###################
#                 #
#   Annotation    #
#                 #
###################


function annotate()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    #Prokka
    prokka  --outdir "${annotation}"/"$sample" \
            --force \
            --prefix "$sample" \
            --kingdom "$kingdom" \
            --genus "$genus" \
            --species "$species" \
            --strain "$sample" \
            --gram "$gram" \
            --locustag "$locustag" \
            --compliant \
            --centre "$centre" \
            --cpus $((cpu/maxProc)) \
            --rfam \
            "$1"

    #extract hypothetical proteins
    cat "${annotation}"/"${sample}"/"${sample}".faa | \
        awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
        grep --no-group-separator -A 1 -F "hypothetical protein" \
        > "${annotation}"/"${sample}"/"${sample}"_hypoth.faa

    echo -e ""$sample" hypothetical proteins (round1): $(cat "${annotation}"/"${sample}"/"${sample}".faa | grep -ic "hypothetical")" \
        | tee -a "${logs}"/log.txt

    #make sure $BLASTDB is set in environment variables
    # export BLASTDB=/media/3tb_hdd/db/nr:/media/3tb_hdd/db/nt
    blastp -query "${annotation}"/"${sample}"/"${sample}"_hypoth.faa \
        -db nr \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -evalue 1e-30 \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -num_threads $((cpu/maxProc)) \
        > "${annotation}"/"${sample}"/"${sample}"_hypoth.blastp

    #Fetch the fasta entry of the hits that do not contain "hypothetical"
    #Re-filter for evalues
    cat "${annotation}"/"${sample}"/"${sample}"_hypoth.blastp | \
        grep -viF "hypothetical" | \
        awk '{if($12 < 1e-30) {print}}' | \
        cut -f 2 | \
        cut -d "|" -f4 \
        > "${annotation}"/"${sample}"/accession.list

    #Download the sequences
    perl "${scripts}"/http_post.pl \
        "${annotation}"/"${sample}"/accession.list \
        "${annotation}"/"${sample}"/extra_hits.fasta
    # esearch -db protein -query "$(cat "${spadesOut}"/annotation/accession.list)" | \
    #     efetch -db protein -format fasta \
    #     > "${spadesOut}"/annotation/extra_hits.fasta

    #TODO -> Cleanup sequence titles. E.g. "MULTISPECIES: "
    # Make first letter uppercase
    # remove duplicate entries
    cat "${annotation}"/"${sample}"/extra_hits.fasta \
        | sed 's/MULTISPECIES: //' \
        | sed 's/ ./\U&/' \
        | awk 'BEGIN {RS=">"} NR>1 {sub("\n","\t"); gsub("\n",""); print RS$0}' \
        | sort -uk1,1 \
        | tr "\t" "\n" \
        > "${annotation}"/"${sample}"/extra_hits_nodup.fasta

    #relaunch Prokka annotation with the new positive blast hit fasta file as reference
    prokka  --outdir "${annotation}"/"$sample" \
            --force \
            --prefix "$sample" \
            --kingdom "$kingdom" \
            --genus "$genus" \
            --species "$species" \
            --strain "$sample" \
            --gram "$gram" \
            --locustag "$locustag" \
            --compliant \
            --centre "$centre" \
            --cpus $((cpu/maxProc)) \
            --rfam \
            --proteins "${annotation}"/"${sample}"/extra_hits_nodup.fasta \
            "$1"

    echo -e ""$sample" hypothetical proteins(round2): $(cat "${annotation}"/"${sample}"/"${sample}".faa | grep -ic "hypothetical")" \
        | tee -a "${logs}"/log.txt
}

export -f annotate

find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" | \
    parallel    --bar \
                --env annotate \
                --env annotation \
                --env kingdom \
                --env genus \
                --env species \
                --env gram \
                --env locustag \
                --env centre \
                --env cpu \
                --env maxProc \
                --env scripts \
                --jobs "$maxProc" \
                "annotate {}"


### Resfinder

function run_resfinder ()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${amr}"/"${sample}"/"$resfinder_db" ] || mkdir -p "${amr}"/"${sample}"/"$resfinder_db"

    perl "${prog}"/resfinder/resfinder.pl \
        -d "${prog}"/resfinder/resfinder_db/ \
        -a "$resfinder_db" \
        -i "$1" \
        -o "${amr}"/"$sample"/"$resfinder_db" \
        -k 90 \
        -l 60
}

export -f run_resfinder

for h in $(find "${prog}"/resfinder/resfinder_db -type f -name "*.fsa"); do
    export resfinder_db=$(sed 's/\.fsa//' <<< $(basename "$h"))

    find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" | \
        parallel --env run_resfinder \
            --env resfinder_db \
            "run_resfinder {}"
done

#Check if any hit
find "${amr}" -type f -name "results_tab.txt" \
    -exec cat {} \; | sed -n '1d' | tee "${amr}"/resfinder_hits.txt


### Phaster

#trim assemblies
function phaster_trim()
{
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

    # http://phaster.ca/instructions
    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # if more than one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${phaster}"/assemblies/"${sample}"_trimmed2000.fasta
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # if only one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            1500 \
            "$1" \
            > "${phaster}"/assemblies/"${sample}"_trimmed1500.fasta
    else
        echo "No assembly for "$sample""  # Should not get here!
        exit 1
    fi
}

#make function available to parallel
export -f phaster_trim  # -f is to export functions

 # To store trimmed assemblies for phaster submission
[ -d "${phaster}"/assemblies ] || mkdir -p "${phaster}"/assemblies 

#run trimming on multiple assemblies in parallel
find "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel  --bar \
                --env phaster_trim \
                --env prog \
                'phaster_trim {}'

function phasterSubmit ()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
    wget --post-file="$i" \
        http://phaster.ca/phaster_api?contigs=1 \
        -O "${phaster}"/"${sample}".json \
        -o "${phaster}"/"${sample}"_wget.log
}

# Submit to phaster sequencially
for i in $(find "${phaster}"/assemblies -type f -name "*.fasta"); do
    phasterSubmit "$i"
done

python3 ~/scripts/checkPhasterServer.py -f "$phaster"


#TODO -> wrap all relevant information in a nice PDF report
