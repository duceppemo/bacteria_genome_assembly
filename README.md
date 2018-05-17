# Description

This pipeline uses gzipped paired-end fastq files as input. Sample names are taken from everything that is before the first underscore in the file name. So filenames must have an underscore for the pipeline to work as is. Also, there are a few hard coded paths in the first section of the script (User Defined). These paths must be updated to fit your installation.

# Dependencies

Software:
- Centrifuge
- Java
- FastQC
- BBTools
- SPAdes assembler
- BWA-mem
- Samtools
- Pilon
- QUAST
- CD-HIT-EST
- Prokka
- Mash
- Mauve
- BLAST
- KAT
- Jellyfish
- GenomeScope
- Conda
- Blobtools
- ResFinder
- Unicycler

Note that some of these tools also have dependencies that are not covered here. An active internet connection is also required.

Scripts:
- removesmallscontigs.pl
- findClosest.sh
- http_post.pl
- checkPhasterServer.py

Perl modules:
TODO

Python modules:
TODO

# Workflow
1. Fastq files are reordered and compressed to increase speed of downstream processes. The script can be modified to skip this step if disk space is an issue (set the "fastq" folder to the "reads" folder).
2. Quality checks (QC) on raw data. Estimates genome coverage based on kmer frequency.
3. Read pre-processing (quality and adapter trimming, contaminant removal, error-correction and overlapping paired-end merging).
4. *de novo* assembly (with Unicycler)
5. Assembly trimming
6. Assembly QC and statistics
7. Assembly annotation
8. Prophage detection

