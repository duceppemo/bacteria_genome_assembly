#!/bin/bash

version="0.3.1"


######################
#                    #
#    User Defined    #
#                    #
######################


#Analysis folder
export baseDir=""${HOME}"/analyses/salmonella_pigeon"

#reads
export reads="/media/30tb_raid10/data/salmonella_pigeon/merged"

#program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

# Centrifuge DB to use
# export centrifuge_db="/media/30tb_raid10/db/centrifuge/2017-10-12_bact_vir_h"

# #Kraken DB to use
export kraken2_db="/media/30tb_raid10/db/kraken2/standard"  # refseq
# export kraken2_db="/media/30tb_raid10/db/kraken2/nt"

#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=8

#Annotation
export kingdom="Bacteria"
export genus="Listeria"
export species="monocytogenes"
export gram="pos"
export locus_tag="$RANDOM"
export centre="OLF"

# For assembly trimming
export smallest_contig=1000
export size=3000000


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
export qc=""${baseDir}"/qc"
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
conda activate quast
if hash quast 2>/dev/null; then
    quast --version | tee -a "${logs}"/log.txt
    conda deactivate
else
    echo >&2 "QUAST was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#CD-HIT-EST
if hash cd-hit-est 2>/dev/null; then
    cd-hit-est -h | head -n 1 | tr -d "=" | sed 's/^[ \t]*//;s/[ \t]*$//' | tee -a "${logs}"/log.txt
else
    echo >&2 "cd-hit-est was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
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

# RGI/CARD

# PHASTER
# local version

# Kraken2
if hash kraken2 2>/dev/null; then
    kraken2 --version | grep -F "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "kraken2 was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# kraken2 database
# if [ -s "${centrifuge_db}".1.cf ]; then
#     echo "Centrifuge database: "$(basename "$centrifuge_db")"" | tee -a "${logs}"/log.txt
# else
#     echo "Could no find the provided Centrifude database. Aborting." | tee -a "${logs}"/log.txt
#     exit 1
# fi

#add space after prog version
echo -e "\n" | tee -a "${logs}"/log.txt


###################
#                 #
#   Fastq files   #
#                 #
###################


# create symbolic links
function create_symlink()
{
    name=$(basename "$1")
    sample=$(basename "$1" | cut -d '_' -f 1)
    [ -d "${fastq}"/"$sample" ] || mkdir -p "${fastq}"/"$sample"
    ln -s "$1" "${fastq}"/"${sample}"/"$name"
}
export -f create_symlink
find "$reads" -type f -name "*fastq.gz" |
    parallel    --bar \
                --env create_symlink \
                --env fastq \
                'create_symlink {}'


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

find -L "$fastq" -type f -name "*_R1*fastq.gz" \
    | parallel  --bar \
                --env run_fastqc \
                --env qc \
                --env cpu \
                --env maxProc \
                --jobs "$maxProc" \
                "run_fastqc {} "${qc}"/fastqc/raw"

conda activate multiqc

#Merge all FastQC reports together
multiqc \
    -o "${qc}"/fastqc/raw \
    -n merged_reports.html \
    "${qc}"/fastqc/raw

conda deactivate


### Contamination finder ###

[ -d  "${qc}"/confindr ] || mkdir -p "${qc}"/confindr

conda activate confindr
confindr.py \
    -i "$reads" \
    -o "${qc}"/confindr \
    -d /media/bioinfo/30tb_raid10/db/rMLST \
    -t "$cpu" \
    -Xmx "${mem}"g

conda deactivate


################
#              #
#   Trimming   #
#              #
################


# Using fastp
function run_fastp()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))    

    # --merge \
    # --merged_out "${merged}"/"${sample}"/"${sample}"_merged.fastq.gz \
    # --correction \
    fastp \
        --in1 "$r1" \
        --in2 "$r2" \
        --out1 "${trimmed}"/"${sample}"_R1.fastq.gz \
        --out2 "${trimmed}"/"${sample}"_R2.fastq.gz \
        --length_required 64 \
        --html "${logs}"/fastp/"${sample}".html \
        --thread $((cpu/maxProc))
}

export -f run_fastp

conda activate fastp  # conda create -n fastp -c bioconda fastp parallel

[ -d "${logs}"/fastp ] || mkdir -p "${logs}"/fastp

find -L "$fastq" -type f -name "*_R1*" |
    parallel    --bar \
                --env cpu \
                --env maxProc \
                --env trimmed \
                --env merged \
                --env logs \
                --jobs "$maxProc" \
                'run_fastp {}'

conda deactivate


# QC trimmed reads

[ -d "${qc}"/fastqc/trimmed ] || mkdir -p "${qc}"/fastqc/trimmed

find "$trimmed" -type f -name "*_R1*fastq.gz" \
    | parallel  --bar \
                --env run_fastqc \
                --env qc \
                --env cpu \
                --env maxProc \
                --jobs "$maxProc" \
                "run_fastqc {} "${qc}"/fastqc/trimmed"

conda activate multiqc

#Merge all FastQC reports together
multiqc \
    -o "${qc}"/fastqc/trimmed \
    -n merged_reports.html \
    "${qc}"/fastqc/trimmed

conda deactivate



################
#              #
#   Assembly   #
#              #
################


function assemble_skesa()
{

    r1="$1"
    r2=$(echo "$1" | sed 's/_R1/_R2/')

    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    skesa \
        --cores $((cpu/maxProc)) \
        --memory "$mem" \
        --fastq "${r1}","${r2}" \
        --gz \
        --contigs_out "${assembly}"/"${sample}".fasta \
        --min_contig 1000

    # Reformat fasta
}

export -f assemble_skesa

find "$trimmed" -type f -name "*_R1*" | \
    parallel    --bar \
                --env assemble_skesa \
                --env cpu \
                --env maxProc \
                --env assembly \
                --env scripts \
                --jobs "$maxProc" \
                "assemble_skesa {}"


###################
#                 #
#   Assembly QC   #
#                 #
###################


# coverage
function get_coverage_skesa()  # unsing unmerged reads only
{
    r1="$1"
    r2=$(echo "$r1" | sed s'/_R1/_R2/')
    sample=$(basename "$1" | cut -d '_' -f 1)

    [ -d "${qc}"/coverage/"$sample" ] || mkdir -p "${qc}"/coverage/"$sample"

    genome="${assembly}"/"${sample}".fasta

    bwa index "$genome"

    #Align corrected paired-end
    rg_pe="@RG\tID:"${sample}"\tCN:"${centre}"\tLB:NexteraXT\tPL:ILLUMINA\tSM:"${sample}""

    bwa mem -t $((cpu/maxProc)) -M -R "$rg_pe" "$genome" "$r1" "$r2" | \
    samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
    samtools sort -@ $((cpu/maxProc)) - | \
    samtools rmdup - "${qc}"/coverage/"${sample}"/"${sample}".bam

    samtools index "${qc}"/coverage/"${sample}"/"${sample}".bam

    #Average genome depth of coverage
    average_cov=$(samtools depth \
        "${qc}"/coverage/"${sample}"/"${sample}".bam  \
        | awk '{sum+=$3} END { print sum/NR}')

    printf "%s\t%.*f\n" "$sample" 0 "$average_cov" | tee -a "${qc}"/coverage/average_cov.tsv
}

export -f get_coverage_skesa

[ -d "${qc}"/coverage ] || mkdir -p "${qc}"/coverage
echo -e "Sample\tAverage_Cov" > "${qc}"/coverage/average_cov.tsv

find "$trimmed" -type f -name "*_R1*fastq.gz" | \
    parallel    --bar \
                --env get_coverage_skesa \
                --env cpu \
                --env maxProc \
                --env qc \
                --env assembly \
                --env centre \
                --jobs "$maxProc"  \
                "get_coverage_skesa {}"

#clean bwa index files
find "$assembly" -type f ! -name "*.fasta" -exec rm {} \;


# Qualimap
function run_qualimap()
{
    sample=$(cut -d "." -f 1 <<< $(basename "$1"))
    
    [ -d "${qc}"/qualimap/"$sample" ] || mkdir -p "${qc}"/qualimap/"$sample"

    qualimap bamqc \
        --paint-chromosome-limits \
        -bam "$1" \
        --java-mem-size="${mem}"G \
        -nt $((cpu/maxProc)) \
        -outdir "${qc}"/qualimap/"$sample" \
        -outfile "${sample}" \
        -outformat HTML

    # Remove bam files
    # rm -rf "${qc}"/coverage/"$sample"
}

export -f run_qualimap

conda activate qualimap  # conda create -n qualimap -c bioconda qualimap parallel

find "${qc}"/coverage -type f -name "*.bam" |
parallel    --bar \
            --env run_qualimap \
            --env qc \
            --env mem \
            --env cpu \
            --env maxProc \
            --jobs "$maxProc" \
            'run_qualimap {}'

conda deactivate

# delete mapping files
rm -rf "${qc}"/coverage


### quast ###

#All the genomes compared
declare -a genomes=()
for i in $(find "$assembly" -type f -name "*.fasta"); do 
    genomes+=("$i")
done

conda activate quast

quast.py \
    --output-dir "${qc}"/quast/all \
    --threads "$cpu" \
    --min-contig "$smallest_contig" \
    --est-ref-size "$size" \
    --no-icarus \
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
        --no-icarus \
        $1  # don't put in quotes
}

#make function available to parallel
export -f run_quast  # -f is to export functions

#run paired-end merging on multiple samples in parallel
find "$assembly" -type f -name "*.fasta" \
    | parallel  --bar \
                --env run_quast \
                --env maxProc \
                --env cpu \
                --env qc \
                --env smallest_contig \
                --jobs "$maxProc" \
                "run_quast {}"

conda deactivate


###################
#                 #
#   Annotation    #
#                 #
###################


function annotate()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    # #Prokka
    # prokka  --outdir "${annotation}"/"$sample" \
    #         --force \
    #         --prefix "$sample" \
    #         --kingdom "$kingdom" \
    #         --genus "$genus" \
    #         --species "$species" \
    #         --strain "$sample" \
    #         --gram "$gram" \
    #         --locustag "$locustag" \
    #         --compliant \
    #         --centre "$centre" \
    #         --cpus $((cpu/maxProc)) \
    #         --rfam \
    #         "$1"

    # #extract hypothetical proteins
    # cat "${annotation}"/"${sample}"/"${sample}".faa | \
    #     awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
    #     grep --no-group-separator -A 1 -F "hypothetical protein" \
    #     > "${annotation}"/"${sample}"/"${sample}"_hypoth.faa

    # echo -e ""$sample" hypothetical proteins (round1): $(cat "${annotation}"/"${sample}"/"${sample}".faa | grep -ic "hypothetical")" \
    #     | tee -a "${logs}"/log.txt


    # export BLASTDB=$BLASTDB:/media/bioinfo/30tb_raid10/db/nr:/media/bioinfo/30tb_raid10/db/nt

    # diamond makedb \
    #     --in /media/bioinfo/30tb_raid10/db/diamond/nr/nr.fasta \
    #     --db /media/bioinfo/30tb_raid10/db/diamond/nr/nr

    # diamond blastp \
    #     --db /media/bioinfo/30tb_raid10/db/diamond/nr/nr \
    #     --query "${annotation}"/"${sample}"/"${sample}"_hypoth.faa \
    #     --outfmt 5 \
    #     --out "${annotation}"/"${sample}"/"${sample}"_hypoth.faa.xml \
    #     --threads "$cpu" \
    #     --evalue 1e-10 \
    #     --tmpdir /media/bioinfo/2TB_NVMe/tmp \
    #     --block-size 12 \
    #     --index-chunks 1

    # python3 "${scripts}"/bda.py \
    #     -i "${annotation}"/"${sample}"/"${sample}"_hypoth.faa \
    #     -x "${annotation}"/"${sample}"/"${sample}"_hypoth.faa.xml
    #     -o "${annotation}"/"${sample}"/extra_hits.fasta

    # # Make first letter uppercase
    # # remove duplicate entries
    # cat "${annotation}"/"${sample}"/extra_hits.fasta \
    #     | sed 's/MULTISPECIES: //' \
    #     | sed 's/ ./\U&/' \
    #     | awk 'BEGIN {RS=">"} NR>1 {sub("\n","\t"); gsub("\n",""); print RS$0}' \
    #     | sort -uk1,1 \
    #     | tr "\t" "\n" \
    #     > "${annotation}"/"${sample}"/extra_hits_nodup.fasta

    #relaunch Prokka annotation with the new positive blast hit fasta file as reference
    prokka  --outdir "${annotation}"/"$sample" \
            --force \
            --prefix "$sample" \
            --kingdom "$kingdom" \
            --genus "$genus" \
            --species "$species" \
            --strain "$sample" \
            --gram "$gram" \
            --locustag "$RANDOM" \
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

conda activate prokka

find "$assembly" -type f -name "*.fasta" | \
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

conda deactivate


###########
#         #
#   AMR   #
#         #
###########


### Resfinder

function run_resfinder ()
{
    # https://bitbucket.org/genomicepidemiology/resfinder/overview
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${amr}"/resfinder/"$sample" ] || mkdir -p "${amr}"/resfinder/"$sample"

    python3 "${prog}"/resfinder/resfinder.py \
        -i "$1" \
        -o "${amr}"/resfinder/"$sample" \
        -p "${prog}"/resfinder/resfinder_db/ \
        -t 0.9 \
        -l 0.6 #\
        # 1> >(tee "${amr}"/resfinder/"${sample}"/"${sample}"_resfinder.txt)

    rm -rf "${amr}"/resfinder/"${sample}"/tmp
}

export -f run_resfinder

find -L "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" | \
    parallel    --env run_resfinder \
                --env resfinder_db \
                "run_resfinder {}"

# Create merged report
echo -e 'Sample\tResistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\tPosition in reference\tContig\tPosition in contig\tPhenotype\tAccession no.' \
    > "${amr}"/resfinder/resfinder_merged.tsv.tmp
    
for i in $(find "${amr}"/resfinder -name "*results_tab.txt"); do
    # sample name is folder name
    sample=$(basename $(dirname "$i"))

    # Add a leading column with sample name
    cat "$i" \
        | sed -e '1d' \
        | awk -F $'\t' -v s="$sample" 'BEGIN {OFS = FS} {print s,$0}' \
        >> "${amr}"/resfinder/resfinder_merged.tsv.tmp
done

# sort by sample name (column 1), then by Identity (column 3)
(cat "${amr}"/resfinder/resfinder_merged.tsv.tmp | head -n 1;
    cat "${amr}"/resfinder/resfinder_merged.tsv.tmp | sed -e '1d' | sort -t $'\t' -k1,1 -k3,3) \
    > "${amr}"/resfinder/resfinder_merged.tsv

rm "${amr}"/resfinder/resfinder_merged.tsv.tmp

### RGI (CARD)

function run_rgi()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${amr}"/rgi/"$sample" ] || mkdir -p "${amr}"/rgi/"$sample"

    "${prog}"/rgi-4.2.2/./rgi main \
        -i "$1" \
        -o  "${amr}"/rgi/"${sample}"/"$sample" \
        -t 'contig' \
        -a 'BLAST' \
        -n $((cpu/maxProc)) \
        --clean \
        -d "chromosome"
}

export -f run_rgi

find -L "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel    --bar \
                --env run_rgi \
                --env resfinder_db \
                "run_rgi {}"

# Create merged report
echo -e 'Sample\tORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_typeSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug\tClass\tResistance\tMechanism\tAMR\tGene\tFamily\tPredicted_DNA\tPredicted_Protein\tCARD_Protein_Sequence\tPercentage_Length_of_Reference_Sequence\tID\tModel_ID' \
    > "${amr}"/rgi/rgi_merged.tsv.tmp
    
for i in $(find "${amr}"/rgi -name "*.txt"); do
    # sample name is folder name
    sample=$(basename $(dirname "$i"))

    # Add a leading column with sample name
    cat "$i" \
        | sed -e '1d' \
        | awk -F $'\t' -v s="$sample" 'BEGIN {OFS = FS} {print s,$0}' \
        >> "${amr}"/rgi/rgi_merged.tsv.tmp
done

# sort by Sample name (column 1), then by Cutt_off (column 7)
(cat "${amr}"/rgi/rgi_merged.tsv.tmp | head -n 1;
    cat "${amr}"/rgi/rgi_merged.tsv.tmp | sed -e '1d' | sort -t $'\t' -k1,1 -k7,7) \
    > "${amr}"/rgi/rgi_merged.tsv

rm "${amr}"/rgi/rgi_merged.tsv.tmp


################
#              #
#   Prophage   #
#              #
################


### Phaster

#trim assemblies
function phaster_trim()
{
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

    # http://phaster.ca/instructions
    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # If more than one contig
        # Remove contigs smaller than 2000 bp from assembly
        perl "${scripts}"/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${phaster}"/assemblies/"${sample}"_trimmed2000.fasta
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # If only one contig
        # Check if contig is at least 1500 bp
        seqlen=$(cat "$1" | awk '!/^>/ {l+=length($0)} END {print l}')
        if [ "$seqlen" -lt 1500 ]; then
            echo "Assembly is one contig, but smaller than 1500bp! Skipping."
        else
            ln -s "$1" "${phaster}"/assemblies/"${sample}".fasta
        fi
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
find -L "$ordered" -type f -name "*trimmed"${smallest_contig}"_ordered.fasta" \
    | parallel  --bar \
                --env phaster_trim \
                --env prog \
                'phaster_trim {}'

# function phasterSubmit ()
# {
#     sample=$(basename "$1" | cut -d '_' -f 1)

#     # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
#     wget --post-file="$i" \
#         http://phaster.ca/phaster_api?contigs=1 \
#         -O "${phaster}"/"${sample}".json \
#         -o "${phaster}"/"${sample}"_wget.log
# }

# # Submit to phaster sequencially
# counter=0
# total=$(find "${phaster}"/assemblies -type f -name "*.fasta" | wc -l)
# for i in $(find "${phaster}"/assemblies -type f -name "*.fasta"); do
#     let counter+=1
#     sample=$(cut -d '_' -f 1 <<< $(basename "$1"))
#     echo -ne "Submitting "${sample}" ("${counter}"/"${total}")"\\r
#     phasterSubmit "$i"
# done

# # TODO -> check that all submission went OK
# #         if not, resubmit

# python3 "${scripts}"/checkPhasterServer.py \
#     --check \
#     -o "$phaster"


# function extract_fasta()
# {
#     sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

#     #Only get the fasta file out of the zip
#     unzip -p \
#         -j "${phaster}"/"${sample}"_phaster.zip \
#         "phage_regions.fna" \
#         > "${phaster}"/"${sample}"_phages.fasta

#     #Add sample name and entry number to fasta header
#     sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"_phages.fasta
# }

# export -f extract_fasta

# find "$phaster" -type f -name "*_phaster.zip" |
# parallel    --bar \
#             --env extract_fasta \
#             --env phaster \
#             'extract_fasta {}'


# Run local version of phaster
# /media/30tb_raid10/db/blast/prophages/prophage_virus_filtered_fixedheaders.db
for i in $(find -L "${phaster}"/assemblies -type f -name "*.fasta"); do
    sample=$(basename "$1" "_trimmed2000.fasta")
    [ -d "${phaster}"/"$sample" ] || mkdir -p "${phaster}"/"$sample"

    perl /home/bioinfo/prog/phaster-app/scripts/phaster.pl \
        -c \
        -i "$i" \
        -o "${phaster}"/"${sample}"/"$sample"
done

# Cleanup
rm -rf "${phaster}"/assemblies

#reformat phaster output
function reformat_output ()
{
    # #unzip phaster result file
    # [ -d "${phaster}"/"$sample" ] || mkdir -p "${phaster}"/"$sample"
    # unzip "$1" -d "${phaster}"/"${sample}"

    #Convert summary.txt to tsv format
    cat "${1}"/summary.txt \
        | sed '1,/gi|/d' \
        | awk '{$1=$1;print}' \
        | sed '/^-/d' \
        | tr " " "\t" \
        > "${1}"/summary.tsv

    #Convert detail.txt to tsv format
    cat "${1}"/detail.txt \
        | sed -e '/^-/d' -e 's/  /@/g' \
        | tr -s '@' \
        | sed 's/@ /@/' \
        | tr -s '[:space:]+' \
        | tr "@" "\t" \
        > "${1}"/detail.tsv
}

for i in $(find "$phaster" -mindepth 1 -type d ! -name assemblies); do
    reformat_output "$i"
done

# Rename output file
for i in $(find "${phaster}" -type f -name "region_DNA.txt"); do
    d=$(dirname "$i")
    sample=$(basename "$d")

    mv "$i" "${d}"/"${sample}"_phages.fasta
done

#Add sample name and entry number to fasta header
for i in $(find "${phaster}" -name "*_phages.fasta"); do
    sample=$(basename "$i" "_phages.fasta")
    sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"/"${sample}"_phages.fasta
done


#TODO -> wrap all relevant information in a nice PDF report
