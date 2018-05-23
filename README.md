# Description

*illuminaPE_parallel_assemblies.sh* uses gzipped paired-end fastq files as input. Sample names are taken from everything that is before the first underscore in the file name. So filenames must have an underscore for the pipeline to work as is. Also, there are a few hard coded paths in the first section of the script (User Defined). These paths must be updated to fit your installation.

# Dependencies

Software:
- Centrifuge (https://ccb.jhu.edu/software/centrifuge/manual.shtml)
- Java
- FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- BBTools (https://jgi.doe.gov/data-and-tools/bbtools/)
- SPAdes assembler (cab.spbu.ru/software/spades/)
- BWA-mem (http://bio-bwa.sourceforge.net/)
- SAMtools/BCFtools/HTSlib (http://www.htslib.org/download/)
- Pilon (https://github.com/broadinstitute/pilon/wiki)
- QUAST (http://quast.sourceforge.net/quast)
- Prokka (https://github.com/tseemann/prokka)
- Mash (http://mash.readthedocs.io/en/latest/)
- Mauve (http://darlinglab.org/mauve/mauve.html)
- BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi_
- KAT 2.4+ (http://kat.readthedocs.io/en/latest/)
- GenomeScope (https://github.com/schatzlab/genomescope)
- Conda (https://conda.io/miniconda.html)
- Blobtools (https://blobtools.readme.io/docs)
- ResFinder (https://bitbucket.org/genomicepidemiology/resfinder)
- Unicycler (https://github.com/rrwick/Unicycler)

Note that some of these tools also have dependencies that are not covered here. An active internet connection is also required. Many of these tools can now be installed through Conda too. For such tool, I prefer to create a virtual environment. In this pipeline, blobtools and KAT are installed this way.

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

