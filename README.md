Allopolyploid AD1 RNA-seq Mapping and Analysis
---
## Introduction
The genome sequences of allopolyploid cotton species have been completed. I would like to use the [Zhang et al.] (http://www.nature.com/nbt/journal/v33/n5/full/nbt.3207.html/) G. hirsutum var. TM1 genome sequence as reference to conduct RNA-seq mapping and downstream analysis.

RNA-seq of 56 TM1 tissue samples were sequenced using Illumina Hiseq2000: developing seeds and cotyledon from several stages, roots and stems of 2-week-old plants; petals, torus, pistils, stamens and lower sepals dissected from whole mature flowers; ovules from −3, −1, 0, 1, 3, 5, 10, 20, 25 and 35 days DPA; fibers from 5, 10, 20 and 25 DPA; true leaves of the seedlings treated with salt, PEG, heat and cold. A total of 297.3 Gb of raw RNA-seq data were generated from the 56 libraries, and deposited together with genomic sequences in [PRJNA248163] (http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA248163/). 

## Prerequisites
First, repare reference genome.

    cd jfw-lab/Projects
    mkdir AD1_mapping
    cd AD1_mapping
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/TM1.fasta
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/Gossypium_hirsutum_v1.1.gene.gff3

Then go to NCBI Biobroject (http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA248163/), select all 56 RNA-seq sequences to get SRR IDs. Run batch file `downloadSRR.sh` and `prepareFastq.sh`.

    # download SRA file
    module load sratoolkit/2.5.2
    perfetch SRR1695160
    # convert sra to fastq files
    fastq-dump fastq-dump -I --split-files SRR1695160.sra
    # Checking read quality with FastQC
    module load fastqc/0.11.3
    fastqc SRR1695160_1.fastq
    fastqc SRR1695160_2.fastq
    # trim fastq files
    module load sickle/1.33
    sickle pe -t sanger -f SRR1695160_1.fastq -r SRR1695160_2.fastq -o SRR1695160_1.trimmed.fastq -p SRR1695160_2.trimmed.fastq -s SRR1695160_S.trimmed.fastq
    
To save some storage space, consider compress or delete used fastq files

    rm *sra
    rm *fastqc.zip
    gzip -9 *fastq
    
Double check and record the number of reads in raw fastq and trimmed fastq files
    grep "@SRR" fastq_raw/*fastq -c > count.raw.txt
    grep "@SRR" fastq_trimmed/*fastq -c > count.trimmed.txt
                                                                                                                                               
## RSEM quantification
RSEM ([Lit et al. 2010](http://bioinformatics.oxfordjournals.org/content/26/4/493.long); [Li&Dewey 2011](http://bioinformatics.oxfordjournals.org/content/26/4/493.long/)) is an RNA-Seq transcript quantification program developed in 2009. 

Besides [manual](http://deweylab.github.io/RSEM/README.html/), this tutorial is helpful [here](https://github.com/bli25ucb/RSEM_tutorial/).

RSEM can directly call Bowtie, Bowtie2 or STAR for read mapping, however, RNA-seq reads need to by mapped to transcripts not genomic sequences. `TM1.transcripts.fa` contains all extracted transcript sequences from "exon" features, and `*.bt2` are Bowtie2 indices.
    
    module load rsem/1.2.22
    module load bowtie2/2.2.6
    module load cufflinks/2.2.1 
    
    # convert gff3 to gtf
    gffread Gossypium_hirsutum_v1.1.gene.gff3 -T -o TM1.gtf
    
    # built transcript reference from genome
    rsem-prepare-reference --gtf TM1.gtf --bowtie2 TM1.fasta TM1
    
    # estimate gene and isoform expression from RNA-seq data
    ## rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name
    ## rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 
    ## rsem-calculate-expression [options] --alignments [--paired-end] input reference_name sample_name
    rsem-calculate-expression -p 2 --paired-end --output-genome-bam --bowtie2 --time SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq TM1 SRA/fastq_trimmed/rsem/SRR1695160
    
 
In the quantification command, `-p 8` tells to use 8 threads. By default, RSEM generates an annotated BAM file in transcript coordinates. If we also want the alignments in genomic coordinates, turn on the `--output-genome-bam` option. Note that this option is only available when you build references from a genome. The option `--time` report time consumed by each step of RSEM to "sample_name.time". At the end we provide RSEM with two FASTQ files, which contain the first and second mates of the paired-end reads, and also tell RSEM where the references locate and where to output the results.

    module load gsnap
    gsnap -n 1 -N 1 -Q -t 2 --merge-distant-samechr -d AD1TM1 -D /home/jfw-lab-local/gmapdb/AD1TM1/AD1_TM1/ -A sam SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq > SRA/fastq_trimmed/gsnap/SRR1695160.sam 2>> log

`-n` one path to show, `-N` looks for novel splicing, `-t 8` tells to use 8 threads, `-Q` output protein seq, 

# 3	seeds from 0, 5, and 10 h

