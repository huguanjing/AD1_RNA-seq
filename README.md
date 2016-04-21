# Allopolyploid AD1 RNA-seq Mapping and Analysis
---
## Introduction
The genome sequences of allopolyploid cotton species have been completed. I would like to use the [Zhang et al.] (http://www.nature.com/nbt/journal/v33/n5/full/nbt.3207.html/) G. hirsutum var. TM1 genome sequence as reference to conduct RNA-seq mapping and downstream analysis.

RNA-seq of 56 TM1 tissue samples were sequenced using Illumina Hiseq2000: leaves, roots and stems of 2-week-old plants; petals, torus, pistils, stamens and lower sepals dissected from whole mature flowers; ovules from −3, −1, 0, 1, 3, 5, 10, 20 and 25 days DPA; fibers from 5, 10, 20 and 25 DPA; true leaves of the seedlings treated with salt, PEG, heat and cold. A total of 297.3 Gb of raw RNA-seq data were generated from the 57 libraries, and deposited together with genomic sequences in [PRJNA248163] (http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA248163/). 

## Prerequisites
First, repare reference genome

    cd jfw-lab/Projects
    mkdir AD1_mapping
    cd AD1_mapping
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/TM1.fasta

Then get to NCBI Biobroject (http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA248163/), select all 56 RNA-seq sequences to get SRR IDs.

    bash downloadSRR.sh

## RSEM quantification
