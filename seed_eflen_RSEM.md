Eflen Analysis: bowtie2-RSEM datasets
---

In this analysis, I will use seed transcriptiom of A2, D5 and the in silico synthetic ADs to check the parition of At and Dt reads. The partition results from Bowtie2-RSEM will be compared tp GSNAP-Polycat results.

This analysis is part of the Effective Length project.
    
## bowtie2-RSEM

### Exploratory analysis
Run biowtie-RSEM against A2D5 psuedo transcriptome reference(batch file `bowtie2rsem.sh`).
    
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
    
Check mapping results from log files.
    
    # get total sequenced read number
    grep 'reads; of these:' *log
    # get RSEM alignment rate
    grep 'alignment rate' *log

The following R Analysis explore the results At and Dt reads compared to true A2 and D5 reads.
   
    fileL<- grep("genes.results", list.files(),value=TRUE)
    # make count table for expected read counts
    rm(count)
    for ( file in fileL ) {
        temp <- read.table(file, header=TRUE, sep="\t")
        temp <- temp[,c(1,5)]
        names(temp)[2]<- gsub(".genes.results","",file)
        if (!exists("count")) { count<-temp }
        else count <- merge(count, temp, by="gene_id")
        }
    write.table(count,file="count.txt", row.names=FASLE, sep="\t")
    
    count<-read.table("count.txt", header=TRUE, sep="\t")
    totalS <- colSums(count[,-1])
    AtS <- colSums(count[grep(".A",count$gene_id),-1])
    DtS <- colSums(count[grep(".D",count$gene_id),-1])
    pdf("rsem.pdf")
    barplot(rbind(AtS,DtS), main="Total read counts: At vs Dt",col=c("brown","blue"),legend=c("At","Dt"), las=2)
    
    A2t<-count[grep(".A",count$gene_id), c(2:5,7:13)]
    D5t<-count[grep(".D",count$gene_id), c(25:28, 30:36)]
    A2 <-count[grep(".D",count$gene_id), c(2:5,7:13)] + A2t
    D5 <-count[grep(".A",count$gene_id), c(25:28, 30:36)]    + D5t
    At <-count[grep(".A",count$gene_id), c(14:24)]
    Dt <-count[grep(".D",count$gene_id), c(14:24)]
    names(At) <- gsub("AD","At",names(At))
    names(Dt) <- gsub("AD","Dt",names(Dt))
    
    source('~/Dropbox/Scripts/addTrans.r', chdir = TRUE)
    par(mfrow=c(2,2))
    plot(log2(as.matrix(A2)),log2(as.matrix(D5)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(At)),log2(as.matrix(Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2)),log2(as.matrix(At)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(D5)),log2(as.matrix(Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2/D5)),log2(as.matrix(At/Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2t/D5t)),log2(as.matrix(At/Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2/D5)),log2(as.matrix(A2t/D5t)), pch=16, col=addTrans("black",10))

    dev.off()

### Differential expression analysis
Next, I will conduct differential expression analysis between expected (At vs Dt) and true homoeolog partitions (A2 vs D5).    
    
    # differential analysis between homoeologs
    goraiIDs <- gsub(".A", "",grep(".A",count$gene_id, value=TRUE) )
    ## true exprected
    expr.A2D5 <- cbind(A2, D5)     
    expr.ADs <- cbind(At, Dt)
 
[EBSeq](https://www.biostat.wisc.edu/~kendzior/EBSEQ/) is an empirical Bayesian DE analysis tool developed in UW-Madison, can take variance due to read mapping ambiguity into consideration by grouping isoforms with parent gene’s number of isoforms. In addition, it is more robust to outliers. RSEM also includes EBSeq in its folder named ‘EBSeq’.

    library(EBSeq)
    pairwiseDE.ebseq<-function(expr, contrast) {
        # for 2 conditions, e.g.
        # geneMat <- as.matrix( expr[,coldata$sample %in% c("A2.10", "D5.10" )] )
        geneMat <- as.matrix( expr[,gsub(".R.*","",colnames(expr)) %in% contrast])
        rownames(geneMat)<-goraiIDs    
        # library size factor
        Sizes = MedianNorm(geneMat)
        
        EBOut = EBTest(Data = geneMat, Conditions =as.factor(gsub(".R.","",colnames(geneMat))), sizeFactors=Sizes, maxround=5 )
        # checking convergency, the differences between the 4th and 5th iterations need to be less than 0.01 or 0.001.
        print(EBOut$Alpha)
        print(EBOut$Beta)
        print(EBOut$P)
        #Checking the model fit and other diagnostics: QQP data points should lie on the y = x line for both conditions; DenNHist check the density plot of empirical q's vs the simulated q's, good model fit is expected.
        par(mfrow=c(1,2))
        QQP(EBOut)
        DenNHist(EBOut)
                
        EBDERes=GetDEResults(EBOut, FDR=0.05)  
        # DEfound - DE gene IDs
        # PPEE and PPDE - the posterior probabilities of being EE or DE for each gene. 
        # Status - each gene's status called by EBSeq. 
        print( table(EBDERes$Status)) 
        # calculate fold change
        GeneFC=PostFC(EBOut)
        # str(GeneFC) # PostFC , RealFC , Direction
        # PlotPostVsRawFC(EBOut,GeneFC)
        res= data.frame( Status = EBDERes$Status,PostFC=NA, RealFC=NA, Direction=NA )
        rownames(res) <- goraiIDs
        res[names(GeneFC$PostFC),"PostFC"] <- GeneFC$PostFC
        res[names(GeneFC$PostFC),"RealFC"] <- GeneFC$RealFC
        res[names(GeneFC$PostFC),"Direction"] <- GeneFC$Direction
        
        write.table(res, file=paste("DE/ebseq.",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
    }
    batch<- rbind( c("A2.10", "D5.10" ), c("A2.20", "D5.20" ), c("A2.30", "D5.30" ), c("A2.40", "D5.40" ) )
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.A2D5,x))
    batch<- rbind( c("At.10", "Dt.10" ), c("At.20", "Dt.20" ), c("At.30", "Dt.30" ), c("At.40", "Dt.40" ) )
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.ADs,x))

Because RSEM output expected_count is not integer, the counts need to be rounded for DESeq2 analysis. 

    library(DESeq2)
    pairwiseDE.deseq2<-function(dds, contrast,savePath)
    {
        # DE analysis
        print(contrast)
        ddsPW <-dds[,dds$sample %in% contrast]
        ddsPW$sample<-droplevels(ddsPW$sample)
        res <- results(DESeq(ddsPW))
        print( summary(res,alpha=.05) ) # print results
        write.table(res, file=paste(savePath,"DE/deseq2.",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
    }
    
    expr <- round(expr.A2D5)
    rownames(expr) <- goraiIDs  
    coldata<-data.frame( sample= gsub("[.]R.","", names(expr)), homoeolog = gsub("[.].*", "", names(expr)), dpa = gsub("D.[.]|A.[.]|[.]R.","", names(expr)), rep = gsub(".*[.]", "", names(expr) ) )
    dds <- DESeqDataSetFromMatrix( countData = expr, colData = coldata, design = ~ sample)
    batch<- rbind( c("A2.10", "D5.10" ), c("A2.20", "D5.20" ), c("A2.30", "D5.30" ), c("A2.40", "D5.40" ) )
    apply(batch,1,function(x) pairwiseDE.deseq2(dds,x,savePath = ""))
    
    expr <- round(expr.ADs)
    rownames(expr) <- goraiIDs  
    coldata<-data.frame( sample= gsub("[.]R.","", names(expr)), homoeolog = gsub("[.].*", "", names(expr)), dpa = gsub("D.[.]|A.[.]|[.]R.","", names(expr)), rep = gsub(".*[.]", "", names(expr) ) )
    dds <- DESeqDataSetFromMatrix( countData = expr, colData = coldata, design = ~ sample)
    batch<- rbind( c("At.10", "Dt.10" ), c("At.20", "Dt.20" ), c("At.30", "Dt.30" ), c("At.40", "Dt.40" ) )
    apply(batch,1,function(x) pairwiseDE.deseq2(dds,x,savePath = ""))

Compare EBSeq and DEseq2 results.

    fileL<-list.files("DE/")
    fileL.ebseq <- paste0("DE/", grep("ebseq",fileL, value=TRUE) )
    fileL.deseq2 <-paste0("DE/", grep("deseq2",fileL, value=TRUE) )
    comp <- gsub("DE/ebseq.|.txt","",fileL.ebseq)
    results <-data.frame(deseq2.sig=NA, ebseq.sig=NA, both.sig =NA)
    pdf("EBSeqvsDese2.pdf")
    par(mfrow=c(2,2))
    for( i in comp)
    {
        ebseq  <- read.table(file=grep(i,fileL.ebseq, value=TRUE), header=TRUE)
        deseq2 <- read.table(file=grep(i,fileL.deseq2,value=TRUE), header=TRUE)
        # table(row.names(deseq2) ==row.names(ebseq))
        results[i,] <- c(table(deseq2$padj < 0.05)["TRUE"], table(ebseq$Status == "DE")["TRUE"], table(deseq2$padj < 0.05 & ebseq$Status == "DE")["TRUE"])
        source("~/Dropbox/Scripts/addTrans.r")
        colors <- rep(addTrans("grey",10), length(row.names(ebseq)))
        # color DE genes: both - red, deseq2 - green, ebs - blue, esle=grey
        colors[ebseq$Status == "DE" ] <- addTrans("blue",10)
        colors[deseq2$padj < 0.05]    <- addTrans("green",10)
        colors[ebseq$Status == "DE" & deseq2$padj < 0.05 ] <- addTrans("red",10)
        plot(-deseq2$log2FoldChange, log2(ebseq$PostFC), type="p", pch=16, col=colors, main=i, )
        legend("bottomright", inset=.05, title="signicant DE", c("ebseq","deseq2","both"), fill=c("blue","green","red") )
    }
    dev.off()
    print(results)
    #                  deseq2.sig ebseq.sig both.sig
    # A2.10vsD5.10      12911      9504     8950
    # A2.20vsD5.20       8382      8458     6463
    # A2.30vsD5.30      13429      9832     9288
    # A2.40vsD5.40       6216      4667     3334
    # At.10vsDt.10      11974      9101     8450
    # At.20vsDt.20       7530      7890     5792
    # At.30vsDt.30      12269      9285     8562
    # At.40vsDt.40       5446      4310     2885

**EBSeq detected fewer DEs than DESeq2, and the overlaps make a much bigger proportion, which suggests that EBSeq is more stringent.**

Compare A2vsD5 v.s. AtvsDt and in DEseq2 and ebseq each 

    comp <- c(10,20,30,40)
    results <-data.frame(A2D5.sig=NA, ADs.sig=NA, both.sig =NA, both.not.sig=NA)
    pdf("TruevsExpected.pdf")
    par(mfrow=c(2,2))
    for( i in comp)
    {
        # ebseq results
        ebseq.A2D5  <- read.table(file=grep("A2",grep(i,fileL.ebseq, value=TRUE),value=TRUE), header=TRUE)
        ebseq.ADs   <- read.table(file=grep("At",grep(i,fileL.ebseq, value=TRUE),value=TRUE), header=TRUE)
        results[paste0("ebseq.",i,"dpa"),] <- c(table(ebseq.A2D5$Status == "DE")["TRUE"], table(ebseq.ADs$Status == "DE")["TRUE"], table(ebseq.A2D5$Status == "DE" & ebseq.ADs$Status == "DE")["TRUE"], table(ebseq.A2D5$Status == "EE" & ebseq.ADs$Status == "EE")["TRUE"] )
        colors <- rep(addTrans("grey",10), length(row.names(ebseq.A2D5)))
        # color DE genes: both - red, A2D5 - green, ADs - blue, esle=grey
        colors[ebseq.A2D5$Status == "DE" ] <- addTrans("green",10)
        colors[ebseq.ADs$Status == "DE" ] <- addTrans("blue",10)
        colors[ebseq.ADs$Status == "DE"  & ebseq.A2D5$Status == "DE" ] <- addTrans("red",10)
        plot(log2(ebseq.ADs$PostFC), log2(ebseq.A2D5$PostFC), type="p", pch=16, col=colors, main=paste0("ebseq.",i,"dpa"))
        legend("bottomright", inset=.05, title="signicant DE", c("ADs","A2D5","both"), fill=c("blue","green","red") )
        
        #deseq2
        deseq2.A2D5  <- read.table(file=grep("A2", grep(i,fileL.deseq2,value=TRUE),value=TRUE), header=TRUE)
        deseq2.ADs   <- read.table(file=grep("At", grep(i,fileL.deseq2,value=TRUE),value=TRUE), header=TRUE)
        results[paste0("deseq2.",i,"dpa"),]  <- c(table(deseq2.A2D5$padj < 0.05)["TRUE"], table(deseq2.ADs$padj < 0.05)["TRUE"], table(deseq2.ADs$padj < 0.05 & deseq2.A2D5$padj < 0.05)["TRUE"], table(deseq2.ADs$padj >= 0.05 & deseq2.A2D5$padj >= 0.05)["TRUE"])
        
        colors <- rep(addTrans("grey",10), length(row.names(deseq2.A2D5)))
        # color DE genes: both - red, A2D5 - green, ADs - blue, esle=grey
        colors[deseq2.A2D5$padj < 0.05 ] <- addTrans("green",10)
        colors[deseq2.ADs$padj  < 0.05 ] <- addTrans("blue",10)
        colors[deseq2.A2D5$padj < 0.05  & deseq2.ADs$padj  < 0.05 ] <- addTrans("red",10)
        plot(deseq2.ADs$log2FoldChange , deseq2.A2D5$log2FoldChange, type="p", pch=16, col=colors, main=paste0("deseq2.",i,"dpa"))
        legend("bottomright", inset=.05, title="signicant DE", c("ADs","A2D5","both"), fill=c("blue","green","red") )
    }
    dev.off()
    print(results)
    
    # sensitivity, true-positive rate (TPR), both.sig / A2D5.sig
    # specificity),true-negative rate (TNR)  both.not.sig / (37223-A2D5.not)
    results$sensitivity <- results$both.sig/results$A2D5.sig
    results$specificity <- results$both.not.sig/(37223-results$A2D5.sig)

To detect differential homoeolog expression, do we want to recover as many as possible TRUE DEs (high sensitivity), or we want to maximize the correct number of not DEs (high specificity)?? **DEseq2 shows slightly higher sensitivity, while EBSeq shows much higher specificity.**

results       |A2D5.sig|ADs.sig|both.sig|both.not.sig|sensitivity|specificity 
--------------|--------|-------|--------|------------|-----------|-----------|:
ebseq.10dpa   |  9504  |  9101 |  7942  |    19905   | 0.8356481 |  0.7180995
deseq2.10dpa  | 12911  | 11974 | 11071  |    15557   | 0.8574859 |  0.6398898
ebseq.20dpa   |  8458  |  7890 |  6720  |    21232   | 0.7945141 |  0.7381192
deseq2.20dpa  |  8382  |  7530 |  6819  |    19205   | 0.8135290 |  0.6658923
ebseq.30dpa   |  9832  |  9285 |  8180  |    19502   | 0.8319772 |  0.7119857
deseq2.30dpa  | 13429  | 12269 | 11411  |    16125   | 0.8497282 |  0.6776919
ebseq.40dpa   |  4667  |  4310 |  3591  |    25803   | 0.7694450 |  0.7925728
deseq2.40dpa  |  6216  |  5446 |  4877  |    21152   | 0.7845882 |  0.6821685

## GSNAP-polyCat
Gather read count tables from HTseq-count results

    load('~/jfw-lab/Projects/Eflen_networks/R-01-countDatasets.RData')
    # "A2.Total"  "D5.Total"  "ADs.Total" "A2.At"     "D5.At"     "ADs.At"  "A2.Dt"     "D5.Dt"     "ADs.Dt"    "ADs.AtN"   "ADs.DtN"   "A2.N"  "D5.N"      "ADs.N"     "A2.AD"     "D5.AD"     "ADs.AD"   
    AtS <- colSums(cbind(A2.At,ADs.At,D5.At) )
    DtS <- colSums(cbind(A2.Dt,ADs.Dt,D5.Dt) )
    NS <- colSums(cbind(A2.N,ADs.N,D5.N) )
    XS <- colSums(cbind(A2.AD,ADs.AD,D5.AD) )
    
    pdf("polycat.pdf")
    barplot(rbind(AtS,DtS,NS,XS), main="Total read counts: At vs Dt",col=c("brown","blue", "grey", "pink"),legend=c("At","Dt", "N" ,"X"), las=2)
    
    source('~/Dropbox/Scripts/addTrans.r', chdir = TRUE)
    par(mfrow=c(2,2))
    plot(log2(as.matrix(A2.Total)),log2(as.matrix(D5.Total)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(ADs.At)),log2(as.matrix(ADs.Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.Total)),log2(as.matrix(ADs.At)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(D5.Total)),log2(as.matrix(ADs.Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.Total/D5.Total)),log2(as.matrix(ADs.At/ADs.Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.At/D5.Dt)),log2(as.matrix(ADs.At/ADs.Dt)), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.Total/D5.Total)),log2(as.matrix(A2.At/D5.Dt)), pch=16, col=addTrans("black",10))

    dev.off()

### Differential expression analysis
Run EBseq and DESeq2 and make comparisons as above 

    ## true exprected
    expr.A2D5 <- cbind(A2.Total, D5.Total)     
    expr.ADs <- cbind(ADs.At, ADs.Dt)
    names(expr.ADs)<-gsub("D5","Dt",gsub("A2","At",names(expr.A2D5)))
    goraiIDs <- row.names(A2.Total)
 
    library(EBSeq)
    batch<- rbind( c("A2-10", "D5-10" ), c("A2-20", "D5-20" ), c("A2-30", "D5-30" ), c("A2-40", "D5-40" ) )
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.A2D5,x))
    batch<- rbind( c("At-10", "Dt-10" ), c("At-20", "Dt-20" ), c("At-30", "Dt-30" ), c("At-40", "Dt-40" ) )
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.ADs,x))

    library(DESeq2)
    expr <- expr.A2D5
    coldata<-data.frame( sample= gsub(".R.*","", names(expr)), homoeolog = gsub("-.*", "", names(expr)), dpa = gsub("D.-|A.-|.R.*","", names(expr))) 
    dds <- DESeqDataSetFromMatrix( countData = expr, colData = coldata, design = ~ sample)
    batch<- rbind( c("A2-10", "D5-10" ), c("A2-20", "D5-20" ), c("A2-30", "D5-30" ), c("A2-40", "D5-40" ) )
    apply(batch,1,function(x) pairwiseDE.deseq2(dds,x,savePath = ""))    
    expr <- expr.ADs
    coldata<-data.frame( sample= gsub(".R.*","", names(expr)), homoeolog = gsub("-.*", "", names(expr)), dpa = gsub("D.-|A.-|.R.*","", names(expr))) 
    dds <- DESeqDataSetFromMatrix( countData = expr, colData = coldata, design = ~ sample)
    batch<- rbind( c("At-10", "Dt-10" ), c("At-20", "Dt-20" ), c("At-30", "Dt-30" ), c("At-40", "Dt-40" ) )
    apply(batch,1,function(x) pairwiseDE.deseq2(dds,x,savePath = ""))

## RSEM vs polyCat with respect to DE analysis
**EBSeq detected fewer DEs than DESeq2, and the overlaps make a much bigger proportion, which suggests that EBSeq is more stringent.**

    # polycat
    #                  deseq2.sig ebseq.sig both.sig
    # A2-10vsD5-10      12758      9609     9068
    # A2-20vsD5-20       8304      8289     6495
    # A2-30vsD5-30      13333      9741     9247
    # A2-40vsD5-40       6300      4651     3414
    # At-10vsDt-10      11841      8910     8373
    # At-20vsDt-20       7130      7168     5391
    # At-30vsDt-30      12129      8877     8318
    # At-40vsDt-40       5345      4166     2866

    # rsem
    #                  deseq2.sig ebseq.sig both.sig
    # A2.10vsD5.10      12911      9504     8950
    # A2.20vsD5.20       8382      8458     6463
    # A2.30vsD5.30      13429      9832     9288
    # A2.40vsD5.40       6216      4667     3334
    # At.10vsDt.10      11974      9101     8450
    # At.20vsDt.20       7530      7890     5792
    # At.30vsDt.30      12269      9285     8562
    # At.40vsDt.40       5446      4310     2885


Regardless of using RSEM or polycat, DEseq2 shows slightly higher sensitivity, while EBSeq shows much higher specificity. 

**More importantly, RSEM and polyCat results are quite consistent with respect to homoeolog assignment sensitivity/specificity and homoeolog differential expression results (75-88% overlap). **

**RSEM**      |A2D5.sig|ADs.sig|both.sig|both.not.sig|sensitivity|specificity 
--------------|--------|-------|--------|------------|-----------|-----------|:
ebseq.10dpa   |  9504  |  9101 |  7942  |    19905   | 0.8356481 |  0.7180995
deseq2.10dpa  | 12911  | 11974 | 11071  |    15557   | 0.8574859 |  0.6398898
ebseq.20dpa   |  8458  |  7890 |  6720  |    21232   | 0.7945141 |  0.7381192
deseq2.20dpa  |  8382  |  7530 |  6819  |    19205   | 0.8135290 |  0.6658923
ebseq.30dpa   |  9832  |  9285 |  8180  |    19502   | 0.8319772 |  0.7119857
deseq2.30dpa  | 13429  | 12269 | 11411  |    16125   | 0.8497282 |  0.6776919
ebseq.40dpa   |  4667  |  4310 |  3591  |    25803   | 0.7694450 |  0.7925728
deseq2.40dpa  |  6216  |  5446 |  4877  |    21152   | 0.7845882 |  0.6821685

**Polycat**   |A2D5.sig|ADs.sig|both.sig|both.not.sig|sensitivity|specificity 
--------------|--------|-------|--------|------------|-----------|-----------|:
ebseq.10dpa   |  9609  |  8910 |  8240  |    19716   | 0.8575294 |  0.7139857
deseq2.10dpa  | 12758  | 11841 | 11252  |    15501   | 0.8819564 |  0.6335990
ebseq.20dpa   |  8289  |  7168 |  6434  |    21278   | 0.7762094 |  0.7353978
deseq2.20dpa  |  8304  |  7130 |  6755  |    19323   | 0.8134634 |  0.6681766
ebseq.30dpa   |  9741  |  8877 |  8184  |    19507   | 0.8401601 |  0.7098101
deseq2.30dpa  | 13333  | 12129 | 11663  |    16024   | 0.8747469 |  0.6707409
ebseq.40dpa   |  4651  |  4166 |  3665  |    25482   | 0.7880026 |  0.7823284
deseq2.40dpa  |  6300  |  5345 |  5047  |    21045   | 0.8011111 |  0.6805614
------------------------------------------------------------------

Direct comparison of AtvsDt DEs from RSEM and polycat. 
    setwd("~/jfw-lab/Projects/Eflen/seed_for_eflen_paper/")
    comp <- c(10,20,30,40)
    results <-data.frame(rsem.sig=NA, polycat.sig=NA, both.sig =NA, both.not.sig=NA)
    pdf("RSEMvsPolycat.pdf")
    par(mfrow=c(2,2))
    for( i in comp)
    {
        # rsem results
        rsem.ebseq   <- read.table(file=paste0('rsemTests/DE/', grep('ebseq.At',grep(i,list.files("rsemTests/DE"),  value=TRUE),value=TRUE)), header=TRUE)
        rsem.deseq2  <- read.table(file=paste0('rsemTests/DE/', grep('deseq2.At',grep(i,list.files("rsemTests/DE"), value=TRUE),value=TRUE)), header=TRUE)
        #polycat results
        polycat.ebseq   <- read.table(file=paste0('polycatTests/DE/', grep('ebseq.At',grep(i,list.files("polycatTests/DE"),  value=TRUE),value=TRUE)), header=TRUE)
        polycat.deseq2  <- read.table(file=paste0('polycatTests/DE/', grep('deseq2.At',grep(i,list.files("polycatTests/DE"), value=TRUE),value=TRUE)), header=TRUE)
        
        # write result table
        results[paste0("ebseq.",i,"dpa"),]   <- c(table(rsem.ebseq$Status == "DE")["TRUE"], table(polycat.ebseq$Status == "DE")["TRUE"], table(rsem.ebseq$Status == "DE" & polycat.ebseq$Status == "DE")["TRUE"], table(polycat.ebseq$Status == "EE" & rsem.ebseq$Status == "EE")["TRUE"] )
        results[paste0("deseq2.",i,"dpa"),]  <- c(table(rsem.deseq2$padj < 0.05)["TRUE"], table(polycat.deseq2$padj < 0.05)["TRUE"], table(polycat.deseq2$padj < 0.05 & rsem.deseq2$padj < 0.05)["TRUE"], table(rsem.deseq2$padj >= 0.05 & polycat.deseq2$padj >= 0.05)["TRUE"])


        colors <- rep(addTrans("grey",10), length(row.names(rsem.ebseq)))
        # color DE genes: both - red, polycat - blue, rsem-green, esle=grey
        colors[polycat.ebseq$Status == "DE" ] <- addTrans("blue",10)
        colors[rsem.ebseq$Status == "DE" ] <- addTrans("green",10)
        colors[polycat.ebseq$Status == "DE"  & rsem.ebseq$Status == "DE" ] <- addTrans("red",10)
        plot(log2(polycat.ebseq$PostFC), log2(rsem.ebseq$PostFC), type="p", pch=16, col=colors, main=paste0("AtvsDt: ebseq.",i,"dpa"))
        legend("bottomright", inset=.05, title="signicant DE", c("polycat","rsem","both"), fill=c("blue","green","red") )
        
                    
        colors <- rep(addTrans("grey",10), length(row.names(rsem.deseq2)))
        # color DE genes: both - red, polycat - blue, rsem-green, esle=grey
        colors[polycat.deseq2$padj  < 0.05 ] <- addTrans("blue",10)
        colors[rsem.deseq2$padj < 0.05 ] <- addTrans("green",10)
        colors[polycat.deseq2$padj < 0.05  & rsem.deseq2$padj  < 0.05 ] <- addTrans("red",10)
        plot(polycat.deseq2$log2FoldChange , rsem.deseq2$log2FoldChange, type="p", pch=16, col=colors, main=paste0("AtvsDt: deseq2.",i,"dpa"))
        legend("bottomright", inset=.05, title="signicant DE", c("polycat","rsem","both"), fill=c("blue","green","red") )
    }
    dev.off()
    print(results)
    write.table(results,file="DEcomparisons.txt",sep="\t", row.names=TRUE)

