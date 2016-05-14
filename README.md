Allopolyploid AD1 RNA-seq Mapping and Analysis
---
## Introduction
The genome sequences of allopolyploid cotton species have been completed. I would like to use the [Zhang et al.](http://www.nature.com/nbt/journal/v33/n5/full/nbt.3207.html/) G. hirsutum var. TM1 genome sequence as reference to conduct RNA-seq mapping and downstream analysis, and compare the mapping results to the strategy using diploid reference genome with SNP index.

RNA-seq of 56 TM1 tissue samples were sequenced using Illumina Hiseq2000: developing seeds and cotyledon from several stages, roots and stems of 2-week-old plants; petals, torus, pistils, stamens and lower sepals dissected from whole mature flowers; ovules from −3, −1, 0, 1, 3, 5, 10, 20, 25 and 35 days DPA; fibers from 5, 10, 20 and 25 DPA; true leaves of the seedlings treated with salt, PEG, heat and cold. A total of 297.3 Gb of raw RNA-seq data were generated from the 56 libraries, and deposited together with genomic sequences in [PRJNA248163](http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA248163/).

A note from Udall lab (Page et al., 2016 in press) is highly relevant: 

"Mapping to the diploid sequences for this report is tenable because **1) the AT- and DT-genomes in tetraploids do not have common loci positions and 2) >25% of the draft tetraploid sequences remain formally unanchored to either AT- or DT-genomes.** Much of our study included comparisons between A and D (or AT vs. DT), and the comparisons are only possible in regions present in both A and D genomes, making the draft tetraploid sequences less informative. To improve results based on diploid sequences, we account for the differences between the respective diploid and tetraploid genomes by adjusting the diploid reference sequences to the genotypes observed in the tetraploid species."

"We also mapped reads to the tetraploid TM-1 reference sequence [20]. The numbers of mapped and categorized reads were less than those obtained with PolyDog using the diploid reference sequences. In addition, a significant percentage of the tetraploid sequence was unanchored to either an AT- or DT-genome. Unanchored scaffolds could be due to either partial assembly or mis-assembly. Thus, further analyses did not use the tetraploid sequence as a genome reference (S1 Table). Eventually, additional improvement of the reference tetraploid sequences may provide better rates of read mapping than PolyDog, but **PolyDog is currently the most thorough method of mapping polyploid reads in cotton**."


## Conclusions
1. Both A2 and D5 genomes provide better mapping results than TM1 genome (Zhang et al.) for polyploid RNA-seq reads
2. For extracting read counts from mapping results, especially for paired-end reads, counter (bambam) always counts as signle-ends (and without quality filter?), while HTSeq-count counts alignment pair (concerdant or only one mate mapped) as 1.
3. GSNAP performs quality filtering, and all paired and solo alignment after filtering are accecpted for mapping; using gsnap -n 1 followed by htseq-count should reach maximum RNA-seq usage (countable_readsout of sequencer_raw_output).
4. In order to test the use of RSEM, I will use seed transcriptiom of A2, D5 and the in silico synthetic ADs to check the parition of At and Dt reads in `[seed_eflen_RSEM.md](https://github.com/huguanjing/AD1_RNA-seq/blob/master/seed_eflen_RSEM.md/)`.


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

Only 4.89% of reads are uniquely mapped, and 44.48% are **multireads**. But 50.63% of reads not mapped to reference seems QUITE high!!!! Corresponding to the 12,136,810 (49.36%) pairs with concordant alignment, **the total read count ("SRR1695160.genes.results" column "expected_count") is 11,840,146, accounting for 48.16% of trimmed pairs, or 36.18% of all sequenced pairs.** This result is comparable to GSNAP genomic mapping results (total sequenced => passing quality check => mappable reads => countable reads), where countable reads are about 18% of all sequenced, and 27% of all mapped.  Note that multireads were discarded by counter or HTSeq.

We can also use Samtools `view` to examine mapping results, but option `flagstat` counts alignments not reads, so reads mapping to multiple locations are counted multiple times, which makes it very hard to make sense of.

    samtools view -f 4 -c filename.bam	# count unmapped reads
    samtools view -f 2 -c filename.bam	# count reads mapped in proper pair
    samtools view -F 4 -c filename.bam | cut -f1 | sort | uniq | wc -l	# count mapped reads, after removing unmapped
    samtools view -f 0x40 -F 0x4 filename.bam | cut -f1 | sort | uniq | wc -l #left mate
    samtools view -f 0x80 -F 0x4 filename.bam | cut -f1 | sort | uniq  | wc -l #right mate
 
Then let us focus on the data statistics learned by RSEM. Inside the folder "SRR1695160.stat", '.cnt' contains alignment statistics, with the format and meanings of each field described in '[cnt_file_description.txt] (https://github.com/deweylab/RSEM/blob/master/cnt_file_description.txt/)' under RSEM directory; '.model' stores RNA-Seq model parameters learned from the data, with described format in '[model_file_description.txt] (https://github.com/deweylab/RSEM/blob/master/model_file_description.txt/)' under RSEM directory.

Command rsem-plot-model can also plots the model statistics RSEM learned from the data. The resulting file, SRR1695160.pdf contains plots of learned fragment length distribution, read length distribution, read start position distribution, quality score information and alignment statistics.  

    rsem-plot-model SRR1695160 SRR1695160.pdf

### Visualize aligments
I got tens sets of A2-D5-TM1.At-TM1.Dt ortholog groups from Justin, listed in "gene_ids.txt"), which can be used to examine the alignment results and quantification.

    rsem-plot-transcript-wiggles --show-unique SRR1695160 gene_ids.txt genes.pdf

In the generated figure (genes.pdf), black refers to uniquely aligned reads and red refers to the expected depth from multi-mapping reads. Everything is red as multireads, and the pair of homoeolog genes appear to share almost same sets of mapped reads, so I assume RSEM assign reads counts based on mapping quality.

## GSNAP with PolyDog (failed)

Use GSNAP to conduct maping, option `-n` tells to report one best alignment only, `-N` looks for novel splicing, `-t 2` tells to use 2 threads, `-Q` output protein seq.

`-m, --max-mismatches=FLOAT` defines maximum  number  of  mismatches  allowed (if not specified, then defaults to the ultrafast level of ((readlength+index_interval-1)/kmer-2)). If specified between 0.0 and 1.0, then treated as a fraction of each read length.  Otherwise, treated as an  integral  number of  mismatches  (including  indel  and  splicing  penalties) For RNA-Seq, you may need to increase this value slightly  to  align reads extending past the ends of an exon.

    module load gsnap/20151120
    module load samtools/1.2
    
    gsnap -n 1 -N 1 -Q -t 2 --merge-distant-samechr -d AD1TM1 -D /home/jfw-lab-local/gmapdb/AD1TM1/AD1_TM1/ -A sam SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq > SRA/fastq_trimmed/gsnap/SRR1695160.TM1.sam 2>> log
    
    gsnap -n 1 -N 1 -Q -t 2 --merge-distant-samechr -d D5 -D /home/jfw-lab-local/gmapdb/D5 -A sam SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq > SRA/fastq_trimmed/gsnap/SRR1695160.D5.sam 2>> log
    
    gsnap -n 1 -N 1 -Q -t 2 --merge-distant-samechr -d A2Li -D /home/jfw-lab-local/gmapdb/A2Li/A2_Li/A2Li -A sam SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq > SRA/fastq_trimmed/gsnap/SRR1695160.A2.sam 2>> log
    
    # count pairs with mapped reads
    samtools view -F 4 filename.bam | cut -f1 | sort | uniq | wc -l
    # count concordantly mapped pairs
    samtools view -f 2 filename.bam | cut -f1 | sort | uniq | wc -l
    
    # convert sam to bam, sort and index for PolyCat/PolyDog
    samtools view -Sb SRR1695160.A2.sam > SRR1695160.A2.bam
    samtools sort -n SRR1695160.A2.bam SRR1695160.A2.nsort
    
    # check counter read table
    samtools sort SRR1695160.A2.bam SRR1695160.A2.sort
    samtools index SRR1695160.A2.sort.bam
    module load bambam/1.3
    counter -g /home/jfw-lab-local/gmapdb/A2Li/A2Li.exons.gff SRR1695160.A2.sort.bam > counter.A2.txt
    counter -g ~/jfw-lab/Projects/Eflen/eflen_recheck/D5.13chrs.exon_unnamed.gff SRR1695160.D5.sort.bam > counter.D5.txt
    counter -g ~/jfw-lab/Projects/AD1_mapping/Gossypium_hirsutum_v1.1.exon.gff3 SRR1695160.TM1.sort.bam > counter.TM1.txt

    # check HTSeq read table
    module load python
    htseq-count -f bam --stranded=no -r pos SRR1695160.D5.sort.bam ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > htseq.D5.txt
    tail htseq.D5.txt
    # Among 24,589,150 SAM alignment pairs processed.
    # __no_feature    1,948,227(7.92%) ===== reads not assigned to any feature, aka those mapped to non-genic regions here
    # __ambiguous       416,480(1.69%) ===== assigned to more than one feature, probably reads assigned to region with multiple (anti-sense?) annotation 
    # __too_low_aQual 1,930,408(7.85%) ===== skip all reads with alignment quality lower than the given minimum value, default: 10
    # __not_aligned   2,276,510(9.26%) ===== in the SAM file without alignment
    # __alignment_not_unique  35,698 (0.15%) ==== reads (or read pairs) with more than one reported alignment, recognized from the NH optional SAM field tag. maybe the pair mapped to different chromosome??
    # The rest of counts assigned to genes sum up to 17,981,827 accounting for 73.13% of trimmed pairs.


Clearly, approximately 5% more TM1 reads can be mapped to A2 and D5 genomes than the TM1 genome, which makes TM1 the less favorable reference genome. **Note that Counter doesn't work for PE reads.** For some reasons, counter command failed for A2, and D5 counter result is weird, resulted total read count bigger than total mapped alignment according to samtools flagstat.

Genome                     |      TM1	     |       A2        |      D5
---------------------------|-----------------|-----------------|---------------:
pair with alignment        |19334995 (78.65%)|22254614 (90.52%)|22307621(90.74%)
pair mapped concordantly   |18370740 (74.73%)|20138464 (81.92%)|20175471(82.07%)
Counter total reads (SE)   |33823921 (68.79%)|    error        |46808845(95.20%)
HTSeq total reads (PE)     |        NA       |       NA        |17981827(73.13%)

    module load bambam/1.3
    polyDog -o test -A /home/jfw-lab-local/gmapdb/A2Li/A2genome_13.fasta -B /home/jfw-lab-local/gmapdb/D5/Dgenome2_13.fasta SRR1695160.A2.sort.bam SRR1695160.D5.sort.bam
	
I cannot get polyDog work, got error message:

    polyDog: src/polyDog.cpp:199: int main(int, char**): Assertion `!pair || aln2f.Name == aln2r.Name' failed.
    Aborted
    
Log into [babycrunch.las.iastate.edu] and try the bambam/1.4
module
    ln -s /net/my.files.iastate.edu/ifs/isu/las/research/jfw-lab
    cd jfw-lab/Projects/AD1_mapping/SRA/fastq_trimmed/gsnap/
    
    module load bambam/1.4
    # On biocrunch, 
    polyDog -o test -A A2genome_13.fasta -B Dgenome2_13.fasta SRR1695160.A2.sort.bam SRR1695160.D5.sort.bam


## Homoeolog read partition: GSNAP with PolyCat v.s. RSEM
A fair comparison between these two requires the same reference genome or transcriptome to be used. For TM1 RNA-seq reads, GSNAP mapping can be conducted againest D5 reference with snp index 4.1; while RSEM can quantify homoeolog reads based on read alignment against a psuedo-AD1 transcriptome constructed from D5 genome and snp index 4.1. 

### Prepare psuedo-transcriptomes
First, pseudogenomes of A2, AD1 At and Dt can be constructed using `pseudogenome_by_snp.py`.
         
    cd /jfw-lab/GenomicResources/pseudogenomes
    python pseudogenome_by_snp.py /home/jfw-lab-local/gmapdb/D5/Dgenome2_13.fasta /home/jfw-lab-local/gmapdb/D5/snpindex4.0/D13.snp4.0 pseudo4.0
    python pseudogenome_by_snp.py /home/jfw-lab-local/gmapdb/D5/Dgenome2_13.fasta /home/jfw-lab-local/gmapdb/D5/snpindex4.1/D13.snp4.1 pseudo4.1
    
And reversely, `get_snps_from_genomes.py` generate SNP index from the two resulted fasta files.

    python get_snps_from_genomes.py pseudo4.0.A.fasta pseudo4.0.D.fasta snp40
    python get_snps_from_genomes.py pseudo4.1.A.fasta pseudo4.1.D.fasta snp41
    
Then it doesn't hurt to compare re-constructed snp index to the original one.

    diff snp40 /home/jfw-lab-local/gmapdb/D5/snpindex4.0/D13.snp4.0
    diff snp41 /home/jfw-lab-local/gmapdb/D5/snpindex4.1/D13.snp4.1

Next, the subgenomes need to be combined for build transcriptome references.

    # combine subgenomes
    cat pseudo4.0.A.fasta pseudo4.0.D.fasta > A2D5.fasta
    cat pseudo4.1.A.fasta pseudo4.1.D.fasta > AtDt.fasta
    
    # convert gff3 to gtf
    module load cufflinks/2.2.1
    ln ~/jfw-lab/Projects/Eflen/eflen_recheck/D5.gff
    gffread D5.gff -T -o D5.gtf
    
    # append subgenome tag to Chrs and Gorai IDs
    sed "s/-JGI_221_v2.1/.A/g" D5.gtf | sed 's/\tphy/_A\tphy/' | grep "[.]1[.]A" >At.gtf
    sed "s/-JGI_221_v2.1/.D/g" D5.gtf | sed 's/\tphy/_D\tphy/' | grep "[.]1[.]D" >Dt.gtf
    cat At.gtf Dt.gtf > pseudoAD.gtf
    
    # check to make sure there are the right number (37223 x 2) genes on chromosomes for unnamed
    sed '/Gorai[.]N/d' pseudoAD.gtf | cut -f9 | sort | uniq |wc -l

    # built transcript reference from genome
    module load rsem/1.2.22
    module load bowtie2/2.2.6
    rsem-prepare-reference --gtf pseudoAD.gtf --bowtie2 A2D5.fasta A2D5
    rsem-prepare-reference --gtf pseudoAD.gtf --bowtie2 AtDt.fasta AtDt

    

### RSEM against A2D5 and AtDt psuedo-transcriptome using bowtie2. 

    module load bowtie2/2.2.6
    module load rsem/1.2.22
    module load samtools/1.2
    
    
    mkdir ~/jfw-lab/Projects/AD1_mapping/SRA/fastq_trimmed/rsem/pseudo
    cd ~/jfw-lab/Projects/AD1_mapping/
    rsem-calculate-expression -p 2 --paired-end --output-genome-bam --bowtie2 --time SRA/fastq_trimmed/SRR1695160_1.trimmed.fastq SRA/fastq_trimmed/SRR1695160_2.trimmed.fastq ~/jfw-lab/GenomicResources/pseudogenomes/A2D5 SRA/fastq_trimmed/rsem/pseudo/SRR1695160.A2D5  > SRA/fastq_trimmed/rsem/pseudo/SRR1695160.A2D5.log 2>&1


Together with previous RSEM mapping results for SRR1695160 (raw 32,728,693 spots, and after trimming, 24,584,115 pairs and 8,023,391 singltons were left), approximately 10% more TM1 reads can be mapped to psuedotranscriptome based on D5 reference and SNP index than that of TM1, consistent with GSNAP mappping results

Transcriptome           |      TM1	  |     A2D5        |      AtDt
------------------------|-----------------|-----------------|---------------:
pair mapped concordantly|12136810 (49.36%)|15151380 (61.63%)|15161659(61.67%)
RSEM total read count   |11840146 (48.16%)|14771572 (60.09%)|14781099(60.12%)   
    
    
    
### RSEM against A2D5 and AtDt psuedo-transcriptome using GSNAP.


