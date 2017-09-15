#!/bin/bash

: <<'END'

# Worth having a look:
https://github.com/tseemann/shovill

END


######################
#                    #
#    User Defined    #
#                    #
######################


#Analysis folder
export baseDir=""${HOME}"/analyses/SE_jiewen_virulence"

# export sample="H-H-1"

#reads
export reads="/media/6tb_raid10/data/SE_jiewen"

#program location
export prog=""${HOME}"/prog"
export picard=""${prog}"/picard-tools/picard.jar"

#script location
export scripts=""${HOME}"/scripts"

# Centrifuge DB to use
export db="/media/6tb_raid10/db/centrifuge/p_compressed+h+v"
# db="/media/6tb_raid10/db/centrifuge/nt"

#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=6

#k-mer size for SPAdes assembler (must be odd number(s))
#Should be smaller that minimum trimmed read length
export kmer="21,33,55,77,99,127"

# Assembly trimming
smallest_contig=200

# Annotation
kingdom="Bacteria"
genus="Salmonella "
species="enterica subsp. enterica serovar Heidelberg"
strain="HH1"
tag="XX000"
centre="NCBI"
gram="neg"
export bioproject="PRJNX000000"
export biosample="SAMX00000000"

# Genomic island
export token="21a9f2b6-b6a4-b07a-6b55-9e065e9aa362"  # Make sure it's still active. Must login first. http://www.pathogenomics.sfu.ca/islandviewer/user/token/
export email=marc-olivier.duceppe@inspection.gc.ca



#######################
#                     #
#   Data Stucture     #
#                     #
#######################


#Folder structure
export fastq=""${baseDir}"/fastq"
export logs=""${baseDir}"/logs"
export qc=""${baseDir}"/QC"
export centrifugeOut=""${qc}"/centrifuge/raw"
export kat=""${qc}"/kat"
export genomescope=""${qc}"/genomescope/raw"
export trimmed=""${baseDir}"/trimmed"
export corrected=""${baseDir}"/corrected"
export merged=""${baseDir}"/merged"
export assembly=""${baseDir}"/assembly"
# export spadesOut=""${assembly}"/"${sample}""
export polished="${baseDir}"/polished
blob=""${qc}"/blobtools"
# scaffolded="${baseDir}"/scaffolded
export annotation="${baseDir}"/annotation
export phaster=""${baseDir}"/phaster"
export islandviewer=""${baseDir}"/islandviewer"


#create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$centrifugeOut" ] || mkdir -p "$centrifugeOut"
[ -d "$kat" ] || mkdir -p "$kat"
[ -d "$genomescope" ] || mkdir -p "$genomescope"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$corrected" ] || mkdir -p "$corrected"
[ -d "$merged" ] || mkdir -p "$merged"
[ -d "$assembly" ] || mkdir -p "$assembly"
# [ -d "$spadesOut" ] || mkdir -p "$spadesOut"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$blob" ] || mkdir -p "$blob"
# [ -d "$scaffolded" ] || mkdir -p "$scaffolded"
[ -d "$annotation" ] || mkdir -p "$annotation"
[ -d "$phaster" ] || mkdir -p "$phaster"
[ -d "$islandviewer" ] || mkdir -p "$islandviewer"


######################
#                    #
#     Resources      #
#                    #
######################


#computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"


#######################
#                     #
#   Initiating log    #
#                     #
#######################


#Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#pipeline version
echo -e "\ngenome_assembly.sh version 0.1\n" | tee -a "${logs}"/log.txt  # $0

#check if depenencies are installed
#if so, log version

# Centrifuge
if hash centrifuge 2>/dev/null; then  # if installed
    centrifuge --version | grep "centrifuge-class" | sed -e 's%^.*/%%' -e 's/-class//' | tee -a "${logs}"/log.txt
else
    echo >&2 "centrifuge was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Centrifuge database
if [ -s "${db}".1.cf ]; then
    echo "Centrifuge database: "$(basename "$db")"" | tee -a "${logs}"/log.txt
else
    echo "Could no find the provided Centrifude database. Aborting." | tee -a "${logs}"/log.txt
fi

# FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# BBDuk
if hash bbduk.sh 2>/dev/null; then 
    bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbduk.sh was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# BBmerge
if hash bbmerge.sh 2>/dev/null; then 
    bbmerge.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbmerge.sh was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#Lighter (read error correction)
if hash lighter 2>/dev/null; then
    lighter -v | tee -a "${logs}"/log.txt
else
    echo >&2 "lighter was not found. Aborting." | tee -a "${logs}"/log.txt
    # exit 1
fi

# SPAdes
if hash spades.py 2>/dev/null; then
    spades.py -v 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "spades.py was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# QUAST
if hash quast.py 2>/dev/null; then
    quast.py -v  2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "quast.py was not found. Aborting." | tee -a "${logs}"/log.txt
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

# GenomeScope
# TODO


########################
#                      #
#     Fastq files      #
#                      #
########################


# Added manually to "fastq" folder

# loop
for i in $(find -L "$reads" -type f -name "*.fastq.gz"); do
    sample=$(basename "$i")
    ln -s "$i" "${fastq}"/"$sample"
done


####################
#                  #
#   FastQc - Raw   #
#                  #
####################


[ -d "${qc}"/fastqc ] || mkdir -p "${qc}"/fastqc
[ -e "${qc}"/fastqc/list.txt ] && rm  "${qc}"/fastqc/list.txt
for i in $(find -L "$fastq" -type f -name "*.fastq.gz"); do 
    echo "$i" >> "${qc}"/fastqc/list.txt
done

ARRAY=($( cat  "${qc}"/fastqc/list.txt ))
list=$(echo "${ARRAY[@]}")

[ -d  "${qc}"/fastqc/raw ] || mkdir -p  "${qc}"/fastqc/raw

fastqc \
    --o  "${qc}"/fastqc/raw \
    --noextract \
    --threads "$cpu" \
    $list  # don't put in quotes


########################
#                      #
#   Centrifuge - Raw   #
#                      #
########################


for i in $(find -L "$fastq" -maxdepth 1 -type f -name "*_R1*"); do
    r1="$i"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')
    sample=$(cut -d '_' -f 1 <<< $(basename "$r1"))

    [ -d "${qc}"/centrifuge/raw/"${sample}" ] || mkdir -p "${qc}"/centrifuge/raw/"${sample}"

    centrifuge \
        -p "$cpu" \
        -t \
        --seed "$RANDOM" \
        -x  "$db" \
        -1 "$r1" \
        -2 "$r2" \
        --report-file "${baseDir}"/centrifuge/raw/"${sample}"/"${sample}"_report.tsv \
        > "${qc}"/centrifuge/raw/"${sample}"/"${sample}".tsv

    cat "${qc}"/centrifuge/raw/"${sample}"/"${sample}".tsv | \
        cut -f 1,3 | \
        ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/raw/"${sample}"/"${sample}".html

    #visualize the resutls in Firefow browser
    # firefox file://"${qc}"/centrifuge/raw/"${sample}"/"${sample}".html &
done


###########
#         #
#   KAT   #
#         #
###########


#https://kat.readthedocs.io/en/latest/walkthrough.html

# Output folder
[ -d "${kat}"/reads ] || mkdir -p "${kat}"/reads

# KAT can't take compressed files as input
find -L "$fastq" -type f -name "*.fastq.gz" \
    | parallel --env qc "pigz -d -c {} > "${kat}"/reads/{/.}"
"${sample}"/


function runKat()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")

    sample=$(cut -d "_" -f 1 <<< $(basename $r1))
    
    # GC plot on reads
    kat gcp \
        -t $((cpu/maxProc)) \
        -o "${kat}"/"${sample}"/"$sample"_gcp \
        "$r1" "$r2"

    # Redo gpc plot to add labels to axes
    kat plot density \
        -o "${kat}"/"${sample}"/"$sample"_gcp.mx.png \
        --y_label "GC count" \
        --x_label "K-mer multiplicity (R1 and R2 combines)" \
        --z_label "Distinct K-mers per bin" \
        "${kat}"/"${sample}"/"$sample"_gcp.mx

    # k-mer frequency
    kat hist \
        -t $((cpu/maxProc)) \
        -o "${kat}"/"${sample}"/"$sample"_histo \
        "$r1" "$r2"

    # Redo histo to add labels to axes
    kat plot spectra-hist\
        -o "${kat}"/"${sample}"/"$sample"_histo.png \
        --y_label "Count" \
        --x_label "K-mer multiplicity" \
        "${kat}"/"${sample}"/"$sample"_histo

    # Compare R1 vs R2
    kat comp \
        -t $((cpu/maxProc)) \
        -o "${kat}"/"${sample}"/"$sample" \
        "$r1" "$r2"

    # density plot R1 vs R2
    kat plot density \
        -o "${kat}"/"${sample}"/"${sample}"_r1_vs_r2.png \
        --y_label "K-mer multiplicity for $(basename "${r2%.*}")" \
        --x_label "K-mer multiplicity for $(basename "${r1%.*}")" \
        --z_label "Distinct K-mers per bin" \
        "${kat}"/"${sample}"/"${sample}"-main.mx

    kat plot spectra-mx \
        -o "${kat}"/"${sample}"/"${sample}"_spectra-mx.png \
        --intersection \
        --y_label "27-mer multiplicity for $(basename "${r2%.*}")" \
        --x_label "27-mer multiplicity for $(basename "${r1%.*}")" \
        "${kat}"/"${sample}"/"${sample}"-main.mx
}

#make function available to parallel
export -f runKat

# run samples in parallel
find "${kat}"/reads -type f -name "*R1*.fastq" \
    | parallel  --bar \
                --env runKat \
                --env cpu \
                --env kat \
                --env maxProc \
                --jobs "$maxProc" \
                'runKat {}'

# Cleanup
rm -rf "$kat"/reads
find "$kat" -type f \
    | grep -vE ".png|.stats" \
    | xargs rm -r


####################
#                  #
#   GenomeScope    #
#                  #
####################


function genomeStats()
{
    r1="$1"
    r2=$(sed 's/_R1_/_R2_/' <<< "$r1")

    sample=$(cut -d "_" -f 1 <<< $(basename $r1))

    [ -d "${genomescope}"/"$sample" ] || mkdir -p "${genomescope}"/"$sample"

    zcat "$r1" "$r2"\
        | jellyfish count \
            -s "$mem" \
            -m 21 \
            -C \
            -t $((cpu/maxProc)) \
            -o "${genomescope}"/"${sample}"/"${sample}".jf \
            -L 3 \
            /dev/stdin

    jellyfish histo \
        -t $((cpu/maxProc)) \
        -o "${genomescope}"/"${sample}"/"${sample}".histo \
        "${genomescope}"/"${sample}"/"${sample}".jf

    read_length=$(zcat "$1" \
        | head -n 4 \
        | sed -n '2p' \
        | tr -d "\n" \
        | wc -m)

    Rscript "${prog}"/genomescope/genomescope.R \
        "${genomescope}"/"${sample}"/"${sample}".histo \
        21 \
        "$read_length" \
        "${genomescope}"/"${sample}"
}

export -f genomeStats

find -L "$fastq" -type f -name "*_R1*" |
    parallel    --bar \
                --env genomeStats \
                --env output \
                --env mem \
                --env cpu \
                --env maxProc \
                --env genomescope \
                --jobs "$maxProc" \
                'genomeStats {}'

# Make coverage report
for i in $(find "$genomescope" -type f -name "model.txt"); do
    name="$(dirname "$i")"
    sample=$(basename "$name")

    coverage=$(cat "$i" \
                | grep -E "^kmercov" \
                | tr -s " " \
                | cut -d " " -f 2)
    cov=$(printf '%.0f' "$coverage")
    echo -e ""$sample"\t"$cov"" >> "${genomescope}"/coverages.txt
done


################
#              #
#   Prinseq    #
#              #
################


# # output folder
# [ -d "${qc}"/prinseq/reads ] || mkdir -p "${qc}"/prinseq/reads

# # 
# zcat $(find "$fastq" -type f -name "*.fastq.gz") | \
# prinseq-lite \
#     -verbose \
#     -fastq stdin \
#     -graph_data "${qc}"/prinseq/reads/"${sample}".gd \
#     -out_good null \
#     -out_bad nul

# prinseq-graphs \
#     -i "${qc}"/prinseq/reads/"${sample}".gd \
#     -png_all \
#     -o "${qc}"/prinseq/reads


#################
#               #
#   Trimming    #
#               #
#################


function trim()
{
    #sequence nomenclature:
    # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')
    sample=$(basename "$r1" | cut -d '_' -f 1)

    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        ref="${prog}"/bbmap/resources/nextera.fa.gz \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \
        qtrim=lr trimq=10 \
        minlen=64 \
        out1="${trimmed}"/"${sample}"_Trimmed_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"_Trimmed_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/trimming/"${sample}".txt)
}

#make function available to parallel
export -f trim  # -f is to export functions

#Create report output directory
[ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming

#run trimming on multiple samples in parallel
find -L "$fastq" -type f -name "*.fastq.gz" -name "*_R1*" \
    | parallel  --env trim \
                --env cpu \
                --env maxProc \
                --env memJava \
                --env prog \
                --env trimmed \
                --env logs \
                --jobs "$maxProc" \
                'trim {}'


#############################
#                           #
#   Read Error Correction   #
#                           #
#############################


function correct()
{
    r1="$1"
    filename=$(basename "$r1")
    sample=$(cut -d '_' -f 1 <<< "$filename")

    # lighter
    lighter \
        -od "$corrected" \
        -r "$r1"\
        -K 23 3000000 \
        -t "$cpu" \
        2> >(tee "${logs}"/correct/"${filename%%.*}".txt)

    #rename as before
    mv "${corrected}"/"${filename%%.*}.cor.fq.gz" "${corrected}"/"$filename"
}

#make function available to parallel
export -f correct  # -f is to export functions

#Create report output directory
[ -d "${logs}"/correct ] || mkdir -p "${logs}"/correct

#run correction on multiple samples in parallel
find "$trimmed" -type f -name "*.fastq.gz" \
    | parallel  --env correct \
                --env maxProc \
                --env cpu \
                --env corrected \
                --env logs \
                --jobs "$maxProc" \
                'correct {}'


####################
#                  #
#     Merging      #
#                  #
####################


function merge()
{
    #sequence nomenclature:
    # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
    r1="$1"
    r2=$(echo "$r1" | sed 's/_1P/_2P/')
    sample=$(basename "$r1" | cut -d '_' -f 1)

    bbmerge.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        out="${merged}"/"${sample}"_merged.fastq.gz \
        outu1="${merged}"/"${sample}"_unmerged_1P.fastq.gz \
        outu2="${merged}"/"${sample}"_unmerged_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee -a "${logs}"/merging/"${sample}".txt)
}

#make function available to parallel
export -f merge  # -f is to export functions

#Create report output directory
[ -d "${logs}"/merging ] || mkdir -p "${logs}"/merging

#run paired-end merging on multiple samples in parallel
# find -L "$trimmed" -type f -name "*.fastq.gz" -name "*_1P*" \
find -L "$corrected" -type f -name "*.fastq.gz" -name "*_1P*" \
    | parallel  --env merge \
                --env maxProc \
                --env cpu \
                --env memJava \
                --env merged \
                --env logs \
                --jobs "$maxProc" \
                'merge {}'


#####################
#                   #
#     Assembly      #
#                   #
#####################


function assemble()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    m1=""${merged}"/"${sample}"_merged.fastq.gz"  # merged
    u11="$1" # ""${merged}"/"${sample}"_unmerged_1P.fastq.gz"  # unmerged R1
    u12=$(sed 's/_1P/_2P/' <<< "$u11")  #""${merged}"/"${sample}"_unmerged_2P.fastq.gz"  # unmerged R2

    spadesOut=""${assembly}"/"${sample}""

    spades.py \
        --only-assembler \
        -t $((cpu/maxProc)) \
        -m "$mem" \
        -k "$kmer" \
        --careful \
        --s1 "$m1" \
        --pe1-1 "$u11" \
        --pe1-2 "$u12" \
        -o "$spadesOut"
}

#make function available to parallel
export -f assemble  # -f is to export functions

find -L "$merged" -type f -name "*.fastq.gz" -name "*_1P*" \
    | parallel  --env assemble \
                --env maxProc \
                --env cpu \
                --env memJ \
                --env merged \
                --env kmer \
                --env logs \
                --jobs "$maxProc" \
                'assemble {}'


#################
#               #
#   Polishing   #
#               #
#################


# Using Pilon
function polish()
{
    genome="$1"

    path=$(dirname "$genome")
    sample=$(basename "$path")

    #remap reads onto assembly
    bwa index "$genome"

    [ -d "${polished}"/"$sample" ] || mkdir -p "${polished}"/"$sample"

    #map merged (single end) reads  
    bwa mem -x intractg -t $((cpu/maxProc)) -r 1 -a -M "$genome" "${merged}"/"${sample}"_merged.fastq.gz | \
        samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
        samtools sort -@ $((cpu/maxProc)) -m 10G -o "${polished}"/"${sample}"/"${sample}"_merged.bam -

    # remove duplicates for sinle-end reads
    java -Xmx64g -jar "$picard" MarkDuplicates \
        INPUT="${polished}"/"${sample}"/"${sample}"_merged.bam \
        OUTPUT="${polished}"/"${sample}"/"${sample}"_merged_nodup.bam \
        METRICS_FILE="${polished}"/"${sample}"/"${sample}"_merged_duplicates.txt \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true

    # Merge and sort all merged files
    samtools sort -@ $((cpu/maxProc)) -m 10G \
        -o "${polished}"/"${sample}"/"${sample}"_merged_sorted.bam \
        "${polished}"/"${sample}"/"${sample}"_merged_nodup.bam

    # index bam file
    samtools index "${polished}"/"${sample}"/"${sample}"_merged_sorted.bam

    #map unmerged (paired-end) reads
    bwa mem -x intractg -t $((cpu/maxProc)) -r 1 -a -M "$genome" "${merged}"/"${sample}"_unmerged_1P.fastq.gz "${merged}"/"${sample}"_unmerged_2P.fastq.gz | \
        samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
        samtools sort -@ $((cpu/maxProc)) -m 10G -o "${polished}"/"${sample}"/"${sample}"_unmerged.bam -

    # remove duplicates for paired-end reads
    samtools rmdup \
        "${polished}"/"${sample}"/"${sample}"_unmerged.bam \
        "${polished}"/"${sample}"/"${sample}"_unmerged_nodup.bam

    # Merge and sort all unmerged files from bbmerge
    samtools sort -@ $((cpu/maxProc)) -m 10G \
        -o "${polished}"/"${sample}"/"${sample}"_unmerged_sorted.bam \
        "${polished}"/"${sample}"/"${sample}"_unmerged_nodup.bam

    # index bam file
    samtools index "${polished}"/"${sample}"/"${sample}"_unmerged_sorted.bam

    #cleanup
    # find "${polished}"/"${sample}" -type f -name "*.nodup.bam" -delete
    # find "${polished}"/"${sample}" -type f -name "*.merged.bam" -delete

    #Correct contigs using pilon based on the Illumina reads
    # java "$memJava" -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit \
    java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar \
        --threads $((cpu/maxProc)) \
        --genome "$genome" \
        --unpaired "${polished}"/"${sample}"/"${sample}"_merged_sorted.bam \
        --frags "${polished}"/"${sample}"/"${sample}"_unmerged_sorted.bam \
        --outdir "${polished}"/"${sample}" \
        --output "${sample}"_pilon \
        --changes

    #clean up
    ls "${polished}"/"$sample" | grep -v "pilon" | xargs rm
}

export -f polish

find "$assembly" -maxdepth 2 -type f -name "scaffolds.fasta" \
    | parallel  --bar \
                --env polish \
                --env merged \
                --env cpu \
                --env maxProc \
                --env polished \
                --env picard \
                --env memJava \
                --env prog \
                --jobs "$maxProc" \
                'polish {}'


#########################
#                       #
#   Assembly trimming   #
#                       #
#########################


# Remove contigs smaller than 200bp
perl "${prog}"/phage_typing/removesmallscontigs.pl \
    "$smallest_contig" \
    "${polished}"/"${sample}"_pilon.fasta \
    > "${polished}"/"${sample}"_pilon"${smallest_contig}".fasta



:<<'END'

###################
#                 #
#   Scaffolding   #
#                 #
###################


# Map trimmed paired-end reads to assembly
for i in $(find "$trimmed" -type f -name "*Trimmed_1P.fastq.gz"); do
    r1="$i"
    r2=$(sed 's/_1P/_2P/' <<< "$r1")
    name=$(cut -d "_" -f 1 <<< $(basename "$i"))

    bwa mem -t "$cpu" -r 1 -a -M "$genome" "$r1" "$r2" | \
        samtools view -@ "$cpu" -b -h - | \
        samtools sort -@ "$cpu" -m 10G -o "${scaffolded}"/"${name}".bam -

    # remove duplicates for paired-end reads
    samtools rmdup \
        "${scaffolded}"/"${name}".bam \
        "${scaffolded}"/"${name}"_nodup.bam

    # Cleanup
    rm "${scaffolded}"/"${name}".bam
done

# Merge bam files
samtools merge -@ "$cpu" \
    - \
    $(find "$scaffolded" -type f -name "*_nodup.bam") | \
    samtools sort -@ "$cpu" -m 10G -o "${scaffolded}"/"${sample}"_all.bam -

# index bam files
samtools index "${scaffolded}"/"${sample}"_all.bam

#cleanup
rm "${scaffolded}"/*_nodup.bam

#average coverage 
coverage=$(samtools depth "${scaffolded}"/"${sample}"_all.bam \
    | awk '{sum+=$3} END { print sum/NR}')

echo ""${sample}" average coverage: "$coverage"x" | tee -a tee -a "${logs}"/log.txt

#SSPACE
echo -e "PE1-1\tbwa\t"$r1"\t"$r2"\t250\t0.5\tFR" > "${scaffolded}"/sample.list
perl "${prog}"/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl \
    -T "$cpu" \
    -l "${scaffolded}"/sample.list \
    -s "${assembly}"/metassembly/CISA_renamed.fasta \
    -x 1 -m 30 -o 20 \
    -p 1 \
    -b Lw_sspace_extension


# Using OPREA-LG
OPERA-LG \
    "${polished}"/"${sample}"_pilon.fasta \
    "${scaffolded}"/"${sample}"_all.bam \
    "$scaffolded"

# unsing BESST
python "${prog}"/BESST/runBESST \
    -c "${polished}"/"${sample}"_pilon.fasta \
    -f "${scaffolded}"/"${sample}"_all.bam \
    -orientation fr \
    -o "$scaffolded"


# Order conigs based on reference genome
python "${prog}"/pyScaf/pyScaf.py \
    -t "$cpu" \
    -r /home/bioinfo/Desktop/mauve/SLCC2372.fasta \
    -f "${polished}"/"${sample}"_pilon.fasta \
    -o "${scaffolded}"/"${sample}"_ordered.fasta \
    --dotplot png
END


#######################
#                     #
#   Contig ordering   #
#                     #
#######################


# Find ResSeq closest genome

# Folder to store refseq genomes
[ -d "${polished}"/refseq ] || mkdir -p "${polished}"/refseq

# Download all assemblies from species
bash "${scripts}"/get_assemblies.sh \
    -q ""$genus" "$species"" \
    -t refseq \
    -a "Complete Genome" \
    -o "${polished}"/refseq \
    -r

# Use mash to find closest genome from refseq
bash "${scripts}"/findClosest.sh \
    "${polished}"/refseq \
    "${polished}"/"${sample}"_pilon"${smallest_contig}".fasta

# The closest is the top one
closest=$(cat "${polished}"/refseq/"${sample}"_pilon"${smallest_contig}".distance.tsv | head -n 1 | cut -f 1)

#its score out of 1000
score=$(cat "${polished}"/refseq/"${sample}"_pilon"${smallest_contig}".distance.tsv | head -n 1 | cut -f 6)

# dont order contigs if not at least 80% similar
if [ "$score" -gt 800 ]; then

    # Use Mauve in batch mode to order contigs with closest genome
    pigz -d "$closest"
    java "$memJava" -cp "${prog}"/mauve_snapshot_2015-02-13/Mauve.jar \
        org.gel.mauve.contigs.ContigOrderer \
        -output "${polished}"/ordered \
        -ref "${closest%.gz}" \
        -draft "${polished}"/"${sample}"_pilon"${smallest_contig}".fasta

    #fix formating of Mauve output
    #find out how many "alignment" folder there is
    n=$(find "${polished}"/ordered -maxdepth 1 -type d | sed '1d' | wc -l)

    # Need to be fixed to avoid sorting
    perl "${scripts}"/formatFasta.pl \
        -i "${sample}"_pilon"${smallest_contig}".fasta \
        -o "${polished}"/"${sample}"_ordered.fasta \
        -w 80

    #align with progessiveMauve
    [ -d "${qc}"/mauve ] || mkdir -p "${qc}"/mauve

    "${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve \
        --output="${qc}"/mauve/"${sample}"_ordered.xmfa \
        "${closest%.gz}.fna" \
        "${polished}"/"${sample}"_ordered.fasta
else
    perl "${scripts}"/formatFasta.pl \
        -i "${polished}"/"${sample}"_pilon.fasta \
        -o "${polished}"/"${sample}"_ordered.fasta \
        -w 80
fi


###################
#                 #
#   Assembly QC   #
#                 #
###################


#QUAST

quast.py \
    -s \
    -o "${qc}"/quast \
    -t "$cpu" \
    "${polished}"/"${sample}"_ordered.fasta

# #All the genomes compared
# [ -e "${assembly}"/list.txt ] && rm "${assembly}"/list.txt
# for i in $(find "$assembly" -type f -name "*_assembly.fasta"); do 
#     echo "$i" >> "${assembly}"/list.txt
# done

# ARRAY=($(cat "${assembly}"/list.txt))
# list=$(echo "${ARRAY[@]}")

# quast.py \
#     -s \
#     -o "${qc}"/quast \
#     -t "$cpu" \
#     $list  # don't put in quotes

# # Make quast report on individual assembly
# function runQuast()
# {
#     sample=$(basename "$1" | cut -d '_' -f 1)

#     quast.py \
#         -t $((cpu/maxProc)) \
#         -o "${assembly}"/"${sample}"/quast \
#         -s \
#         -L \
#         "$1"  # don't put in quotes
# }

# #make function available to parallel
# export -f runQuast  # -f is to export functions

# #run paired-end merging on multiple samples in parallel
# find "$assembly" -type f -name "*_assembly_trimmed1000.fasta" \
#     | parallel  --env runQuast \
#                 --env maxProc \
#                 --env cpu \
#                 --env assembly \
#                 --env logs \
#                 --jobs "$maxProc" \
#                 "runQuast {}"


###########
#         #
#   KAT   #
#         #
###########


# Create output folder
[ -d "${kat}"/assembly ] || mkdir -p "${kat}"/assembly

find "$trimmed" -type f -name "*.fastq.gz" \
    | parallel --env qc "pigz -d -c {} > "${kat}"/assembly/{/.}"

# Stats per contig
kat sect \
    -t "$cpu" \
    -o "${kat}"/assembly/"$sample" \
    --output_gc_stats \
    "${polished}"/"${sample}"_ordered.fasta \
    "${kat}"/assembly/"${sample}"_Trimmed_?P.fastq

[ -d "${kat}"/assembly/coverage ] || mkdir -p "${kat}"/assembly/coverage
[ -d "${kat}"/assembly/gc ] || mkdir -p "${kat}"/assembly/gc

function kat_profile()
{
    node=$(cut -d "_" -f 2 <<< "$1")  # Works for SPAdes assemblies ("NODE_9_length_5363_cov_1554.3")
    
    kat plot profile \
        -o "${kat}"/assembly/coverage/"${sample}"_NODE_"${node}".png \
        --header=$(tr -d ">" <<< "$1") \
        --y_label "K-mer frequency" \
        --x_label "Position (bp)" \
        "${kat}"/assembly/"${sample}"-counts.cvg

    kat plot profile \
        -o "${kat}"/assembly/gc/"${sample}"_NODE_"${node}".png \
        --header=$(tr -d ">" <<< "$1") \
        --y_label "GC percentage" \
        --x_label "Position" \
        "${kat}"/assembly/"${sample}"-counts.gc
}

# make function available to parallel
export -f kat_profile

for i in $(cat "${polished}"/"${sample}"_ordered.fasta | grep -E "^>"); do
    echo "$i" >> "${kat}"/assembly/contig.list
done

cat "${kat}"/assembly/contig.list \
    | parallel  --env kat_profile \
                --env qc \
                --env sample \
                'kat_profile {}'

# kat plot spectra-mx \
#     --intersection \
#     -o "${kat}"/assembly/"${sample}"_spectra-mx.png \
#     "${kat}"/assembly/"${sample}"-main.mx

# genome assembly analysis
# both paired-end compared as one group
kat comp \
    -t "$cpu" \
    -o "${kat}"/assembly/"${sample}"_pe_vs_asm \
    ""${kat}"/assembly/"${sample}"_Trimmed_?P.fastq" \
    "${polished}"/"${sample}"_ordered.fasta

# Rename axes
kat plot spectra-cn \
    -o "${kat}"/assembly/"${sample}"_pe_vs_asm-main.mx.spectra-cn.png \
    --y_label "Number of distinct K-mers" \
    --x_label "K-mer multiplicity" \
    "${kat}"/assembly/"${sample}"_pe_vs_asm-main.mx

kat_distanalysis.py \
    --plot "${kat}"/assembly/"${sample}"_pe_vs_asm-main.mx \
    | tee "${kat}"/assembly/"${sample}"_decompostion_analysis.stats


# Cleanup
find "${kat}"/assembly -maxdepth 1 -type f \
    | grep -vE ".png|stats" \
    | xargs rm -r


#################
#               #
#   blobtools   #
#               #
#################


genome=""${polished}"/"${sample}"_ordered.fasta"
bwa index "$genome"

# Activate blobtools virtual environment
source "${HOME}"/my_virtualenv/blob_env/bin/activate

# map trimmed and corrected reads to assembly
for i in $(find "$corrected" -type f -name "*_Trimmed_1P.fastq.gz"); do
    r1="$i"
    r2=$(sed 's/_1P/_2P/' <<< "$r1")
    name=$(cut -d "_" -f 1 <<< $(basename "$i"))

    bwa mem -x intractg -t "$cpu" -r 1 -a -M "$genome" "$r1" "$r2" | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${blob}"/"${name}".bam -

    # remove duplicates for paired-end reads
    samtools rmdup \
        "${blob}"/"${name}".bam \
        "${blob}"/"${name}"_nodup.bam

    # index bam files
    samtools index "${blob}"/"${name}"_nodup.bam

    # Cleanup
    rm "${blob}"/"${name}".bam
done

# Merge and sort bam files
if [ $(find "$blob" -type f -name "*._nodup.bam" | wc -l) -gt 1 ]; then
    samtools merge -@ "$cpu" - $(find "$blob" -type f -name "*_nodup.bam" | tr "\n" " ") | \
    samtools sort -@ "$cpu" -m 10G -o "${blob}"/"${sample}".bam -

    #Cleanup
    rm "${blob}"/*nodup.bam*
else
    samtools sort -@ "$cpu" -m 10G -o "${blob}"/"${sample}".bam "${blob}"/"${sample}"_nodup.bam
fi

#Cleanup
rm "${blob}"/*nodup.bam*

#index bam file
samtools index "${blob}"/"${sample}".bam

# Make coverage file from bam
blobtools map2cov \
    -i "$genome" \
    -b "${blob}"/"${sample}".bam \
    -o "${blob}"/

#blast assembly on nt
blastn \
    -task megablast \
    -query "$genome" \
    -db nt \
    -out "${blob}"/"${sample}".blast_nt.tsv \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -culling_limit 5 \
    -evalue 1e-25 \
    -num_threads "$cpu"

# Create BlobDB
blobtools create \
    -i "$genome" \
    -c "${blob}"/"${sample}".bam.cov \
    -t "${blob}"/"${sample}".blast_nt.tsv \
    -o "${blob}"/"${sample}"

# Print BlobDB as a table
blobtools view \
    -i "${blob}"/"${sample}".blobDB.json \
    -o "${blob}"/

# Colour a covplot by GC-categories
# Extract sequenceID and GC
grep -v '^#' "${blob}"/"${sample}".blobDB.table.txt \
    | cut -f 1,3 \
    > "${blob}"/"${sample}".blobDB.id.gc.txt

# Divide into GC categories 
i=0
division=5
while [ "$i" -lt 105 ]; do
    # echo "$i"
    cat "${blob}"/"${sample}".blobDB.id.gc.txt \
        | awk -v count="$i" -v div="$division" '
        BEGIN {OFS = ","}
        {
            if (($2 * 100) >= count && ($2 * 100) < (count + div))
            {
                $2 = count"-"(count + (div - 1))"%"
                print
            }
        }' >> "${blob}"/"${sample}".blobDB.id.gc.catcolour.txt
    let i=$i+division
done

# Generate covplot with colour by %GC categories
blobtools covplot \
    -i "${blob}"/"${sample}".blobDB.json \
    -c "${blob}"/"${sample}".bam.cov \
    --catcolour "${blob}"/"${sample}".blobDB.id.gc.catcolour.txt \
    -o "${blob}"/

# Plot BlobDB as a blobplot
blobtools blobplot \
  -i "${blob}"/"${sample}".blobDB.json \
  -o "${blob}"/

# Deactivate the virtual environment
deactivate


###############
#             #
#   PHASTER   #
#             #
###############


function assemblyTrimm()
{
    name=$(basename "$1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    # http://phaster.ca/instructions
    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # if more than one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${1%.fasta}"_trimmed2000.fasta
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # if only one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            1500 \
            "$1" \
            > "${1%.fasta}"_trimmed1500.fasta
    else
        echo "No assembly for "$sample""  # Should not get here!
        exit 1
    fi
}

#make function available to parallel
export -f assemblyTrimm  # -f is to export functions

#run trimming on multiple assemblies in parallel
find "$polished" -type f -name "*_ordered.fasta" \
    | parallel --env assemblyTrimm --env prog 'assemblyTrimm {}'

#submit to phaster
for i in $(find "$polished" -type f | grep -E "trimmed1500|trimmed2000"); do
    sample=$(basename "$i" | cut -d '_' -f 1)

    # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
    wget --post-file="$i" \
        http://phaster.ca/phaster_api?contigs=1 \
        -O "${phaster}"/"${sample}"_query.json \
        -o "${phaster}"/"${sample}"_wget.log
done


function phasterResults()
{
    name=$(basename "$1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    #Retrieve job ID from json file
    jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')
    # echo ""${sample}": "$jobID""  # debug

    #Check if all submission were successful
    while [ -z "$jobID" ]; do  # if no jobID (unsuccessful submission)
        phasterSubmit "$1"  # resubmit the sample
        jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')
    done

    
    #get status
    wget http://phaster.ca/phaster_api?acc="$jobID" -O "${phaster}"/"${sample}"_status.json
    status=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 2 | cut -d ":" -f 2 | tr -d '"')

    # echo "PHASTER analysis of "$sample" is "$status""

    #check if PHASTER job is finished running
    while [ "$status" != "Complete" ]; do
        waitTime="1m" # sleep

        #check status every X time
        echo "Job status of "$sample" is "$status". Checking status back in "$waitTime"." 
        sleep "$waitTime"

        #get status
        wget http://phaster.ca/phaster_api?acc="$jobID" -O "${phaster}"/"${sample}"_status.json

        #check job status
        status=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 2 | cut -d ":" -f 2 | tr -d '"')
    done

    echo "PHASTER analysis of "$sample" is "$status""

    #get the PHASTER output file
    phasterZip=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 4 | cut -d ":" -f 2 | tr -d '"')

    # Check if was already downloaded
    if [ ! -s "${phaster}"/"${sample}"_phaster.zip ]; then
        wget "$phasterZip" -O "${phaster}"/"${sample}"_phaster.zip

        #Only get the fasta file out of the zip
        unzip -p \
            -j "${phaster}"/"${sample}"_phaster.zip \
            "phage_regions.fna" \
            > "${phaster}"/"${sample}"_phages.fasta

        #Add sample name and entry number to fasta header
        sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"_phages.fasta
    fi
}

#make function available to parallel
export -f phasterResults  # -f is to export functions

find "$polished" -type f | grep -E "trimmed1500|trimmed2000" \
    | parallel --bar --delay 0.3 --env phasterResults --env phaster 'phasterResults {}'

#unzip phaster result file
unzip "${phaster}"/"${sample}"_phaster.zip -d "$phaster"

#Convert summary.txt to tsv format
cat "${phaster}"/summary.txt \
    | sed '1,/gi|/d' \
    | awk '{$1=$1;print}' \
    | sed '/^-/d' \
    | tr " " "\t" \
    > "${phaster}"/summary.tsv

#Convert detail.txt to tsv format
cat "${phaster}"/detail.txt \
    | sed -e '/^-/d' -e 's/  /@/g' \
    | tr -s '@' \
    | sed 's/@ /@/' \
    | tr -s '[:space:]+' \
    | tr "@" "\t" \
    > "${phaster}"/summary.tsv


#########################
#                       #
#   Plasmid detection   #
#                       #
#########################


# Detect contigs part of a plasmid

# # PlasmidFinder -> https://cge.cbs.dtu.dk/services/PlasmidFinder/  ???
[ -d "${qc}"/plasmid ] || mkdir -p "${qc}"/plasmid

# perl  "${prog}"/plasmidfinder/plasmidfinder.pl \
#     -d "${prog}"/plasmidfinder/database \
#     -p enterobacteriaceae,gram_positive \
#     -o "${qc}"/plasmid \
#     -k 95.00 \
#     -i "${polished}"/"${sample}"_assembly_trimmed"${smallest_contig}".fasta

# -task megablast \
# Using blast
blastn \
    -query "${polished}"/"${sample}"_ordered.fasta \
    -db nt \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -evalue 1e-25 \
    -num_threads "$cpu" \
    -culling_limit 5 \
    -max_target_seqs 5 \
    -out "${qc}"/plasmid/"${sample}".blastn

if [[ -n $(cat "${qc}"/plasmid/"${sample}".blastn | grep -i "plasmid") ]]; then
    echo -e "\nPlasmid(s) detected in "$sample" assembly\n" | tee -a "${logs}"/log.txt

    # Create report for plasmid detection
    echo -e "Contig(s) from "$sample" assembly identified as plasmid:\n" | tee -a "${qc}"/plasmid/plasmid.txt
    
    cat "${qc}"/plasmid/"${sample}".blastn \
        | grep -i "plasmid" \
        | cut -f 1,3 \
        | sort -u -k1,1 \
        | tee -a "${qc}"/plasmid/plasmid.txt
else
    echo -e "\nNo plasmid detected in "$sample" assembly\n" | tee -a "${logs}"/log.txt
fi


##################
#                #
#   Annotation   #
#                #
##################


# https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/
prokka  --outdir "${annotation}"/"$bioproject" \
        --force \
        --compliant \
        --prefix "$biosample" \
        --locustag "$tag" \
        --centre "$centre" \
        --kingdom "$kingdom" \
        --genus "$genus" \
        --species "$species" \
        --strain "$strain" \
        --gram "$gram" \
        --cpus "$cpu" \
        "${polished}"/"${sample}"_ordered.fasta



#extract hypothetical proteins
cat "${annotation}"/"${bioproject}"/"${biosample}".faa | \
    awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
    grep --no-group-separator -A 1 -F "hypothetical protein" \
    > "${annotation}"/"${bioproject}"/"${biosample}"_hypoth.faa

echo -e "Number of hypothetical proteins found by Prokka: $(cat "${annotation}"/"${bioproject}"/"${biosample}"_hypoth.faa | grep -ic "hypothetical")" \
    | tee -a "${logs}"/log.txt

echo -n "Blasting hypothetical proteins on NR..."

#make sure $BLASTDB is set in environment variables
# export BLASTDB=/media/3tb_hdd/db/nr:/media/3tb_hdd/db/nt
blastp -query "${annotation}"/"${bioproject}"/"${biosample}"_hypoth.faa \
    -db nr \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -evalue 1e-30 \
    -max_target_seqs 1 \
    -num_threads "$cpu" \
    > "${annotation}"/"${bioproject}"/"${biosample}"_hypoth.blastp

#Fetch the fasta entry of the hits that do not contain "hypothetical"
#Re-filter for evalues
cat "${annotation}"/"${bioproject}"/"${biosample}"_hypoth.blastp | \
    grep -viF "hypothetical" | \
    awk '{if($12 < 1e-30) {print}}' | \
    cut -f 2 | \
    cut -d "|" -f4 \
    > "${annotation}"/"${bioproject}"/accession.list

#Download the sequences
perl "${scripts}"/http_post.pl \
    "${annotation}"/"${bioproject}"/accession.list \
    "${annotation}"/"${bioproject}"/extra_hits.fasta
# esearch -db protein -query "$(cat "${spadesOut}"/annotation/accession.list)" | \
#     efetch -db protein -format fasta \
#     > "${spadesOut}"/annotation/extra_hits.fasta


#relaunch Prokka annotation with the new positive blast hit fasta file as reference
prokka  --outdir "${annotation}"/"$bioproject" \
        --force \
        --compliant \
        --prefix "$biosample" \
        --locustag "$tag" \
        --centre "$centre" \
        --kingdom "$kingdom" \
        --genus "$genus" \
        --species "$species" \
        --strain "$strain" \
        --gram "$gram" \
        --cpus "$cpu" \
        --rfam \
        --proteins "${annotation}"/"${bioproject}"/extra_hits.fasta \
        "${polished}"/"${sample}"_ordered.fasta

echo -e "Number of hypothetical proteins remaining after the BLAST (1e-30): $(cat "${annotation}"/"${bioproject}"/"${biosample}".faa | grep -ic "hypothetical")" \
    | tee -a "${logs}"/log.txt


#######################
#                     #
#   Genomic islands   #
#                     #
#######################


# http://www.pathogenomics.sfu.ca/islandviewer/http_api/

# Query the availability of reference and precomputed genomes
curl http://www.pathogenomics.sfu.ca/islandviewer/rest/genomes/ \
    -H "x-authtoken:"${token}"" \
    > "${islandviewer}"/genomes.list
# # Closest genome available in islandviewer database according to mash
# # "Listeria monocytogenes FSL R2-561, complete genome." -> NC_017546.1

#Find closest genome from their database

#Search Islandviewer database for species of interest and get accession number(s)
cat "${islandviewer}"/genomes.list \
    | grep --no-group-separator -B 1 "$genus $species" \
    | grep "ref_accnum" \
    | cut -d '"' -f 4 \
    > "${islandviewer}"/accession.list

#check if species found in Islandviewer database
if [ -s "${islandviewer}"/accession.list ]; then
    #creat folder to hold related genomes
    [ -d "${islandviewer}"/genomes ] || mkdir -p "${islandviewer}"/genomes

    #Download the genomes
    perl "${scripts}"/get_assemblies.sh \
        -l "${islandviewer}"/accession.list \
        -t "refseq" \
        -o "${islandviewer}"/genomes

    #find closest genome using Mash
    bash "${scripts}"/findClosest.sh \
        "${islandviewer}"/genomes \
        "${polished}"/"${sample}"_ordered.fasta

    # The closest is the top one
    closest=$(cat "${islandviewer}"/genomes/"${sample}"_ordered.distance.tsv | head -n 1 | cut -f 1)

    #its score out of 1000
    score=$(cat "${islandviewer}"/genomes/"${sample}"_ordered.distance.tsv | head -n 1 | cut -f 6)

    export acc_closest=$(zcat "$closest" | head -n 1 | tr -d ">" | cut -d " " -f 1)

    echo -e "\nClosest genome from Islandviewer database is: $(basename "$closest" | cut -d "_" -f 1,2) with a scrore of "$score"/1000" \
        | tee -a "${logs}"/log.txt
else
    echo -e "\n"$genus" "$species" not found in Islandviewer database" \
        | tee -a "${logs}"/log.txt
fi



#BASH

#Submit genome
curl -X POST \
    -H "x-authtoken:${token}" \
    -Fref_accnum="$acc_closest" \
    -Fgenome_file=""${annotation}"/"${bioproject}"/"${biosample}".gbk" \
    -Fgenome_name=""$genus"_"$species"_"$strain"" \
    -Femail_addr="$email" \
    -Fformat_type="GENBANK" \
    http://www.pathogenomics.sfu.ca/islandviewer/rest/submit/

#Retrieve information about already submitted genomes
curl http://www.pathogenomics.sfu.ca/islandviewer/rest/jobs/ \
    -H "x-authtoken:${token}" \
    -o "${islandviewer}"/jobs.list

job_token="ysUgsMwSLnreuhxkrei7gp"

#Query the status specific job
curl http://www.pathogenomics.sfu.ca/islandviewer/rest/job/"${job_token}"/ \
    -H "x-authtoken:${token}" \
    -o "${islandviewer}"/"${job_token}".status

#Download the Genbank including genomic island predictions, and the reordered concatenated contigs for draft genomes
curl http://www.pathogenomics.sfu.ca/islandviewer/rest/job/"${job_token}"/download/genbank/ \
    -H "x-authtoken:${token}" \
    -o "${islandviewer}"/"${biosample}"_islandviewer.gbk

#Download the table summarizing genomic island predictions
curl http://www.pathogenomics.sfu.ca/islandviewer/rest/job/"${job_token}"/download/tab/ \
    -H "x-authtoken:${token}" \
    -o "${islandviewer}"/"${biosample}"_islandviewer.tsv


#Python

#### NOT WORKING ####

#Submit genome
python << END

import requests
import sys
from requests_toolbelt.multipart.encoder import MultipartEncoder

server = "http://www.pathogenomics.sfu.ca/islandviewer"
ext = "/rest/submit/"

mygenome = "$annotation" + '/' + "$bioproject" + '/' + "$biosample" + '.gbk'

multipart_data = MultipartEncoder(
    fields={ 'format_type': 'GENBANK',
             'email_addr': "$email",
             'genome_file': ('filename', open(mygenome, 'rb'), 'text/plain'),
             'ref_accnum': "$acc_closest"}
)

headers={ 'Content-Type': multipart_data.content_type,
          'x-authtoken': "$token"}

r = requests.post(server+ext, headers=headers, data=multipart_data)

if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
status=repr(decoded)

handle=(open $islandviewer/job_id.json)
handle.write(status)

END


#get results

python << END

import requests
import sys
from requests_toolbelt.multipart.encoder import MultipartEncoder

id_file = "$islandviewer/job_id.json"
handle = open(id_file, 'r')
token_id = handle.read()

token_id = token_id.split(',')[2]
token_id = token_id.split("'")[len(token_id.split("'")) - 2]

server = 'http://www.pathogenomics.sfu.ca/islandviewer'
ext1 = '/rest/job/'
ext2 = '/download/genbank/'
ext3 = '/download/tab/'

headers={ 'x-authtoken': "$token" }

r1 = requests.post(server+ext1+token_id+ext2, headers=headers)
r2 = requests.post(server+ext1+token_id+ext3, headers=headers)

if not r1.ok or not r2.ok:
    r1.raise_for_status()
    r2.raise_for_status()
    sys.exit()
else:
    with open("$islandviewer/$sample.gbk", 'w') as f:
        for chunk in r1:
            f.write(chunk)
    with open("$islandviewer/$sample.tsv", 'w') as f:
        for chunk in r2:
            f.write(chunk)
END

