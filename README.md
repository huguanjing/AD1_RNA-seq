Allopolyploid AD1 RNA-seq Mapping and Analysis
---
## Introduction
The genome sequences of allopolyploid cotton species have been completed. __I would like to use the [Zhang et al.](http://www.nature.com/nbt/journal/v33/n5/full/nbt.3207.html/) G. hirsutum var. TM1 genome sequence as reference to conduct RNA-seq mapping and downstream analysis, and compare the mapping results to the strategy using diploid reference genome with SNP index.__

RNA-seq of 56 TM1 tissue samples were sequenced using Illumina Hiseq2000: developing seeds and cotyledon from several stages, roots and stems of 2-week-old plants; petals, torus, pistils, stamens and lower sepals dissected from whole mature flowers; ovules from −3, −1, 0, 1, 3, 5, 10, 20, 25 and 35 days DPA; fibers from 5, 10, 20 and 25 DPA; true leaves of the seedlings treated with salt, PEG, heat and cold. A total of 297.3 Gb of raw RNA-seq data were generated from the 56 libraries, and deposited together with genomic sequences in [PRJNA248163](http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA248163/).

A note from Udall lab (Page et al., 2016 in press) is highly relevant: 

"Mapping to the diploid sequences for this report is tenable because **1) the AT- and DT-genomes in tetraploids do not have common loci positions and 2) >25% of the draft tetraploid sequences remain formally unanchored to either AT- or DT-genomes.** Much of our study included comparisons between A and D (or AT vs. DT), and the comparisons are only possible in regions present in both A and D genomes, making the draft tetraploid sequences less informative. To improve results based on diploid sequences, we account for the differences between the respective diploid and tetraploid genomes by adjusting the diploid reference sequences to the genotypes observed in the tetraploid species."

"We also mapped reads to the tetraploid TM-1 reference sequence [20]. The numbers of mapped and categorized reads were less than those obtained with PolyDog using the diploid reference sequences. In addition, a significant percentage of the tetraploid sequence was unanchored to either an AT- or DT-genome. Unanchored scaffolds could be due to either partial assembly or mis-assembly. Thus, further analyses did not use the tetraploid sequence as a genome reference (S1 Table). Eventually, additional improvement of the reference tetraploid sequences may provide better rates of read mapping than PolyDog, but **PolyDog is currently the most thorough method of mapping polyploid reads in cotton**."

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
RSEM ([Lit et al. 2010](http://bioinformatics.oxfordjournals.org/content/26/4/493.long); [Li&Dewey 2011](http://bioinformatics.oxfordjournals.org/content/26/4/493.long/)) is an RNA-Seq transcript quantification program developed in 2009. When RNA-seq reads do not map uniquely to a single location, it applies the Expectation-Maximization (EM) algorithm to consider mapping uncertainty for accurate abundance estimates.

Besides [manual](http://deweylab.github.io/RSEM/README.html/), this tutorial is helpful [here](https://github.com/bli25ucb/RSEM_tutorial/).

### Run RSEM

RSEM can directly call Bowtie, Bowtie2 or STAR for read mapping, or use provided .sam alignment files from any aligner. It is critical to have aligners to report all valid alignment, not just the single "best". Also, RSEM does not handle indel, local and discordant alignments, so RNA-seq reads need to by mapped to transcripts not genomic sequences, and discordant pair end reads are not allowed. `TM1.transcripts.fa` contains all extracted transcript sequences from "exon" features, and `*.bt2` are Bowtie2 indices. In addition, RSEM does not support estimating using a mix of single-end and 
paired-end reads, YET.
    
    module load rsem/1.2.22
    module load bowtie2/2.2.6
    module load cufflinks/2.2.1 
    module load samtools/1.2
    
    # convert gff3 to gtf
    gffread Gossypium_hirsutum_v1.1.gene.gff3 -T -o TM1.gtf
    
    # built transcript reference from genome
    rsem-prepare-reference --gtf TM1.gtf --bowtie2 TM1.fasta TM1
    
    # estimate gene and isoform expression from RNA-seq data
    ## rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name
    ## rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 
    ## rsem-calculate-expression [options] --alignments [--paired-end] input reference_name sample_name
    rsem-calculate-expression -p 2 --paired-end --output-genome-bam --bowtie2 --time SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq TM1 SRA/fastq_trimmed/rsem/SRR1695160  > SRA/fastq_trimmed/rsem/SRR1695160.log 2>&1
    
In the quantification command, `-p 8` tells to use 8 threads. By default, RSEM generates an annotated BAM file in transcript coordinates. If we also want the alignments in genomic coordinates, turn on the `--output-genome-bam` option. Note that this option is only available when you build references from a genome. The option `--time` report time consumed by each step of RSEM to "sample_name.time". At the end we provide RSEM with two FASTQ files, which contain the first and second mates of the paired-end reads, and also tell RSEM where the references locate and where to output the results.

### Explore the data
The raw fastq file of SRR1695160 contains 32,728,693 spots, and after trimming, 24,584,115 pairs and 8,023,391 singltons were left. 

First the log file recorded bowtie2 command and mapping results. Input option `-q` takes fastq files, `--dpad 0 --gbar 99999999` allows no gap, `--mp 1,1` sets maximum (MX) and minimum (MN) mismatch penalties to 1, `--np 1` sets position penalty to 1, `--score-min L,0,-0.1` sets the minimum-score function f to f(x) = 0 + -0.1 * x, where x is the read length, `-I 1` sets he minimum fragment length for valid paired-end alignments to 1, `-X 1000` sets maximum fragment length for valid paired-end alignments to 1000 bp, `--no-mixed --no-discordant` disable looking for unpaired alignment, `-p 2` uses two threads, `-k 200` searched for at most 200 distinct valid alignment for each read.

    bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -p 2 -k 200 -x TM1 -1 SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq -2 SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq | samtools view -S -b -o SRA/fastq_trimmed/rsem/SRR1695160.temp/SRR1695160.bam -
    [samopen] SAM header is present: 70478 sequences.
    24584115 reads; of these:
    24584115 (100.00%) were paired; of these:
    12447305 (50.63%) aligned concordantly 0 times
    1202393 (4.89%) aligned concordantly exactly 1 time
    10934417 (44.48%) aligned concordantly >1 times
    49.37% overall alignment rate

Only 4.89% of reads are uniquely mapped, and 44.48% are **multireads**. But 50.63% of reads not mapped to reference is QUITE high!!!!

We can also use Samtools `view` to examine mapping results, but option `flagstat` counts reads mapping to multiple locations multiple times, which makes it very hard to make sense of.

    # count unmapped reads
    samtools view -f4 -c SRR1695160.transcript.sorted.bam
 
Then let us focus on the data statistics learned by RSEM. Inside the folder "SRR1695160.stat", '.cnt' contains alignment statistics, with the format and meanings of each field described in '[cnt_file_description.txt] (https://github.com/deweylab/RSEM/blob/master/cnt_file_description.txt/)' under RSEM directory; '.model' stores RNA-Seq model parameters learned from the data, with described format in '[model_file_description.txt] (https://github.com/deweylab/RSEM/blob/master/model_file_description.txt/)' under RSEM directory.

Command rsem-plot-model can also plots the model statistics RSEM learned from the data. The resulting file, SRR1695160.pdf contains plots of learned fragment length distribution, read length distribution, read start position distribution, quality score information and alignment statistics.  

    rsem-plot-model SRR1695160 SRR1695160.pdf

### Visualize aligments
I got tens sets of A2-D5-TM1.At-TM1.Dt ortholog groups from Justin, listed in "gene_ids.txt"), which can be used to examine the alignment results and quantification.
    
    rsem-plot-transcript-wiggles --show-unique SRR1695160 gene_ids.txt genes.pdf

In the generated figure (genes.pdf), black refers to uniquely aligned reads and red refers to the expected depth from multi-mapping reads. Everything is red as multireads, and the pair of homoeolog genes appear to share almost same sets of mapped reads, so I assume RSEM assign reads counts based on mapping quality.

## GSNAP with PolyDog

## GSNAP with PolyCat
Use GSNAP to conduct maping, option `-n` tells to report one best alignment only, `-N` looks for novel splicing, `-t 2` tells to use 2 threads, `-Q` output protein seq.

`-m, --max-mismatches=FLOAT` defines maximum  number  of  mismatches  allowed (if not specified, then defaults to the ultrafast level of ((readlength+index_interval-1)/kmer-2)). If specified between 0.0 and 1.0, then treated as a fraction of each read length.  Otherwise, treated as an  integral  number of  mismatches  (including  indel  and  splicing  penalties) For RNA-Seq, you may need to increase this value slightly  to  align reads extending past the ends of an exon.


    module load gsnap
    module load samtools/1.2
    
    gsnap -n 1 -N 1 -Q -t 2 --merge-distant-samechr -d AD1TM1 -D /home/jfw-lab-local/gmapdb/AD1TM1/AD1_TM1/ -A sam SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq > SRA/fastq_trimmed/gsnap/SRR1695160.sam 2>> log

	samtools view -Sb SRR1695160.sam > SRR1695160.gsnap.bam
	
	samtools flagstat *bam
49174198 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
37921323 + 0 mapped (77.12%:-nan%)
49174198 + 0 paired in sequencing
24589161 + 0 read1
24585037 + 0 read2
36745117 + 0 properly paired (74.72%:-nan%)
37164969 + 0 with itself and mate mapped
756354 + 0 singletons (1.54%:-nan%)
422119 + 0 with mate mapped to a different chr
422119 + 0 with mate mapped to a different chr (mapQ>=5)

	
	
	samtools sort $j.bam $j.sort
	samtools index $j.sort.bam
	
	counter -g Dgenome2_13.gene.gff *fiber.sort.bam > fiber.total.count.GH040716.txt
		
                # polycat is not necessary for diploids, below command for AD1
                echo ''
		polyCat -x 1 -p 1 -s /home/jfw-lab/gmapdb/D5/snpindex4.1/D13.snp4.1  $j.sort.bam
Use samtoo
    samtools flagstat SRR1695160.genome.bam

# 3	seeds from 0, 5, and 10 h
