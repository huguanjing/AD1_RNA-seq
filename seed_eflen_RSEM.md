Eflen Analysis: bowtie2-RSEM datasets
---
    
## Run biowtie-RSEM
    
    cd ~/jfw-lab/Projects/Eflen/seed_for_eflen_paper/
    mkdir bowtie2RSEM
    
    # fastq files 100bp SE 
    ls eflen_seed/ | grep 'cut.fq.gz'
    # A2-10-R1.cut.fq.gz A2-10-R2.cut.fq.gz A2-10-R3.cut.fq.gz
    # A2-20-R1.cut.fq.gz A2-20-R2.cut.fq.gz A2-20-R3.cut.fq.gz
    # A2-30-R1.cut.fq.gz A2-30-R2.cut.fq.gz A2-30-R3.cut.fq.gz
    # A2-40-R1.cut.fq.gz A2-40-R2.cut.fq.gz A2-40-R3.cut.fq.gz
    # AD-10-R1.cut.fq.gz AD-10-R2.cut.fq.gz AD-10-R3.cut.fq.gz
    # AD-20-R1.cut.fq.gz AD-20-R3.cut.fq.gz
    # AD-30-R1.cut.fq.gz AD-30-R2.cut.fq.gz AD-30-R3.cut.fq.gz
    # AD-40-R1.cut.fq.gz AD-40-R2.cut.fq.gz AD-40-R3.cut.fq.gz
    # D5-10-R1.cut.fq.gz D5-10-R2.cut.fq.gz D5-10-R3.cut.fq.gz
    # D5-20-R1.cut.fq.gz D5-20-R2.cut.fq.gz D5-20-R3.cut.fq.gz
    # D5-30-R1.cut.fq.gz D5-30-R2.cut.fq.gz D5-30-R3.cut.fq.gz
    # D5-40-R1.cut.fq.gz D5-40-R2.cut.fq.gz D5-40-R3.cut.fq.gz
    
    # note that I don't need D5-20-R2 in future
    module load bowtie2/2.2.6
    module load rsem/1.2.22
    # -p 4 threads
    
    rsem-calculate-expression -p 4 --bowtie2 --time eflen_seed/A2-10-R1.cut.fq.gz ~/jfw-lab/GenomicResources/pseudogenomes/A2D5 bowtie2RSEM/A2-10-R1 >bowtie2RSEM/A2-10-R1.log 2>&1
    
