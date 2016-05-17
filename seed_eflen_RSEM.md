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
    D5t<-count[grep(".D",count$gene_id), c(25:35)]
    A2 <-count[grep(".D",count$gene_id), c(2:5,7:13)] + A2t
    D5 <-count[grep(".A",count$gene_id), c(25:35)]    + D5t
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
 
EBSeq is an empirical Bayesian DE analysis tool developed in UW-Madison, can take variance due to read mapping ambiguity into consideration by grouping isoforms with parent gene’s number of isoforms. In addition, it is more robust to outliers. RSEM also includes EBSeq in its folder named ‘EBSeq’.

    library(EBSeq)
    pairwiseDE.ebseq<-function(expr, contrast) {
        # for 2 conditions, e.g.
        # geneMat <- as.matrix( expr[,coldata$sample %in% c("A2.10", "D5.10" )] )
        geneMat <- as.matrix( expr[,gsub(".R.","",colnames(expr)) %in% contrast])
        rownames(geneMat)<-goraiIDs    
        # library size factor
        Sizes = MedianNorm(geneMat)
        
        EBOut = EBTest(Data = geneMat, Conditions =as.factor(gsub("[.]R.","",colnames(geneMat))), sizeFactors=Sizes, maxround=5 )
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
    # A2.20vsD5.20       4139      4584     2826
    # A2.30vsD5.30      13429      9832     9288
    # A2.40vsD5.40       3646      3898     1993
    # At.10vsDt.10      11974      9101     8450
    # At.20vsDt.20       7530      7890     5792
    # At.30vsDt.30      12269      9285     8562
    # At.40vsDt.40       5446      4310     2885

At 10 and 30 dpa, EBSeq detected fewer DEs than DESeq2, and the overlaps make a much bigger proportion, which suggests that EBSeq is more stringent. But 20 and 40 dpa results are hard to make comparisons, with ~50% overlap for both.
----------------------------
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

To detect differential homoeolog expression, do we want to recover as many as possible TRUE DEs (high sensitivity), or we want to maximize the correct number of not DEs (high specificity)?? **DEseq2 shows higher sensitivity, while EBSeq shows higher specificity.**

results       |A2D5.sig|ADs.sig|both.sig|both.not.sig|sensitivity|specificity 
--------------|--------|-------|--------|------------|-----------|-----------|:
ebseq.10dpa   |  9504  |  9101 |  7942  |    19905   | 0.8356481 |  0.7180995
deseq2.10dpa  | 12911  | 11974 | 11071  |    15557   | 0.8574859 |  0.6398898
ebseq.20dpa   |  4584  |  7890 |  3239  |    21608   | 0.7065881 |  0.6620301
deseq2.20dpa  |  4139  |  7530 |  3218  |    19290   | 0.7774825 |  0.5830613
ebseq.30dpa   |  9832  |  9285 |  8180  |    19502   | 0.8319772 |  0.7119857
deseq2.30dpa  | 13429  | 12269 | 11411  |    16125   | 0.8497282 |  0.6776919
ebseq.40dpa   |  3898  |  4310 |  2757  |    25894   | 0.7072858 |  0.7770143
deseq2.40dpa  |  3646  |  5446 |  2908  |    21336   | 0.7975864 |  0.6354350


    
    
------------------------------------------------------------------
     
    # prepare network datasets
    col_names  <- c( "dev10.R1", "dev10.R2", "dev10.R3", "dev20.R1", "dev20.R3", "dev30.R1", "dev30.R2", "dev30.R3", "dev40.R1", "dev40.R2", "dev40.R3")
    row_names<- c(grep(".A",count$gene_id, value=TRUE), grep(".D",count$gene_id, value=TRUE))
    
    prepDatasets<-function(a,d) 
    {
        names(a) <- col_names
        names(d) <- col_names
        ad       <- rbind (a,d)
        row.names(ad) <- row_names
        return(ad)
    }
    A2D5      <-  prepDatasets(A2,D5)
    A2D5.tech <-  prepDatasets(A2t,D5t)
    ADs       <-  prepDatasets(At,Dt)
    # double check content: 74446    12
    # save
    save(A2D5, A2D5.tech, ADs, file = "R-01-networkDatasets.rsem.RData")

# rlog transformation need DESeq2
library(DESeq2)
rlogTransformation<-function(x)
{
    count <- round( x[,-1] ) #have to round ADs.ncorrect
    coldata<-data.frame( sample = names(count), dpa = gsub("dev|[.]R.","", names(count)), rep = gsub(".*[.]","", names(count)) )
    rownames(count) <-x$gene
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    # count.rld <- assay(rld)
    return(rld)
}
networks.rld <- list(A2D5      = rlogTransformation(A2D5),
                     A2D5.tech = rlogTransformation(A2D5.tech),
                     ADs       = rlogTransformation(ADs),
                  ADs.ncorrect = rlogTransformation(ADs.ncorrect)   )
# save list
networks  <- list(A2D5      = A2D5,
                  A2D5.tech = A2D5.tech,
                  ADs       = ADs,
               ADs.ncorrect = ADs.ncorrect   )
# coldata
count <- networks[[1]][,-1]
coldata<-data.frame( sample = names(count), dpa = as.numeric(gsub("dev|[.]R.","", names(count))), rep = gsub(".*[.]","", names(count)) )
coldata
#    sample dpa rep
#  dev10.R1  10  R1
#  dev10.R2  10  R2
#  dev10.R3  10  R3
#  dev20.R1  20  R1
#  dev20.R3  20  R3
#  dev30.R1  30  R1
#  dev30.R2  30  R2
#  dev30.R3  30  R3
#  dev40.R1  40  R1
#  dev40.R2  40  R2
#  dev40.R3  40  R3
#save
save(coldata, networks, networks.rld, file = "R-01-networkDatasetsListed.RData")


# do some plots
# We can also build the PCA plot from scratch using ggplot2. This is done by asking the plotPCA function to return the data used for plotting rather than building the plot. See the ggplot2 documentation for more details on using ggplot.
library(genefilter)
library(ggplot2)
sumPCA<-
function (x, intgroup = "condition", ntop = 500)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
    1, paste, collapse = " : "))
    return(pca)
}

nnames<-names(networks)
pdf("s1.exploratory.plots.pdf")
for(i in 1:4)
{
    # raw counts
    nn <- networks[[i]]
    # boxplot
    boxplot(log2((nn[,-1])), main=paste(nnames[i],"raw"))
    
    # rld
    rld <- networks.rld[[i]]
    nn.rld <- assay(rld)
    # boxplot
    boxplot(log2((nn.rld)), main=paste(nnames[i],"rlog"))
    
    # PCA
    print(plotPCA(rld, intgroup =  "dpa") )
    # Here, we have used the function plotPCA which comes with DESeq2. The two terms specified by intgroup are the interesting groups for labeling the samples; they tell the function to use them to choose colors.
    
    # We can also build the PCA plot from scratch using ggplot2. This is done by asking the plotPCA function to return the data used for plotting rather than building the plot. See the ggplot2 documentation for more details on using ggplot.
    data<-sumPCA(rld, intgroup = c( "dpa"))
    # summary(data)
    ### Importance of components:
    ### PC1     PC2     PC3      PC4      PC5     PC6     PC7     PC8    PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16
    ### Standard deviation     50.4710 21.0492 15.2206 13.39625 11.68574 9.38357 7.38384 6.89884 6.7829 6.11396 5.69427 4.70881 4.59109 4.24778 4.01172 3.87475
    ### Proportion of Variance  0.6256  0.1088  0.0569  0.04408  0.03354 0.02163 0.01339 0.01169 0.0113 0.00918 0.00796 0.00545 0.00518 0.00443 0.00395 0.00369
    ### Cumulative Proportion   0.6256  0.7345  0.7913  0.83543  0.86897 0.89059 0.90398 0.91567 0.9270 0.93615 0.94412 0.94956 0.95474 0.95917 0.96312 0.96681
    # library(ggplot2)
    pp<-qplot(PC1, PC2, main=paste(nnames[i],"rlog", "top 500"),
    # shape=coldata$genome,
    color=as.factor(coldata$dpa), data=as.data.frame(data$x)) +
    xlab(paste("PC1:", summary(data)$importance[2,1]*100, "% variance")) +
    ylab(paste("PC2:", summary(data)$importance[2,2]*100, "% variance"))
    plot(pp)
    # From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.
    
    data<-sumPCA(rld, intgroup = c("dpa"), ntop=37223)
    # summary(data)
    ### Importance of components:
    ### PC1     PC2     PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10     PC11     PC12
    ### Standard deviation     94.5976 63.7282 54.8723 52.64956 40.00305 32.31420 28.16078 26.19804 24.86588 22.63824 19.82023 18.47903
    ### Proportion of Variance  0.3079  0.1398  0.1036  0.09539  0.05507  0.03593  0.02729  0.02362  0.02128  0.01764  0.01352  0.01175
    ### Cumulative Proportion   0.3079  0.4477  0.5513  0.64669  0.70175  0.73768  0.76497  0.78859  0.80987  0.82750  0.84102  0.85277
    pp<-qplot(PC1, PC2, main=paste(nnames[i],"rlog", "top 37223"), color=as.factor(coldata$dpa ), data=as.data.frame(data$x)) +
    xlab(paste("PC1:", summary(data)$importance[2,1]*100, "% variance")) +
    ylab(paste("PC2:", summary(data)$importance[2,2]*100, "% variance"))
    # From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.
    plot(pp)
    
    
}
dev.off()



    
    # false positive, 

    #######################################
    # similiar compared to polycat results
    load('/Volumes/jfw-lab/Projects/Eflen_networks/R-01-countDatasets.RData')
    A2.At.S <-colSums(A2.At[,-1])
    A2.Total.S <-colSums(A2.Total[,-1])
    A2.N.S <-A2.Total.S - A2.At.S
    
    
    D5.Dt.S <-colSums(D5.Dt[,-1])
    D5.Total.S <-colSums(D5.Total[,-1])
    D5.N.S <-D5.Total.S - D5.Dt.S
    
    ADs.At.S <-colSums(ADs.At[,-1])
    ADs.Dt.S <-colSums(ADs.Dt[,-1])
    ADs.Total.S <-colSums(ADs.Total[,-1])
    ADs.N.S <-ADs.Total.S - ADs.At.S- ADs.Dt.S
    
    pAtS <- c(A2.At.S, ADs.At.S, D5.At.S<-rep(0,11) )
    pDtS <- c(A2.Dt.S<-rep(0,11), ADs.Dt.S, D5.Dt.S)
    pNS  <- c(A2.N.S,  ADs.N.S,  D5.N.S)
    
    pdf("~/Downloads/polyCat.pdf")
    tt<-rbind(pAtS,pDtS,pNS)
    colnames(tt)<-names(pNS)
    pdf("~/Downloads/polyCat.pdf")
    barplot(tt, main="Total read counts -PolyCat: At vs Dt vs N",col=c("brown","blue","grey"),legend=c("At","Dt","N"), las=2)
    
    
    par(mfrow=c(2,2))
    plot(log2(as.matrix(A2.Total[,-1])),log2(as.matrix(D5.Total[,-1])), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(ADs.At[,-1])),  log2(as.matrix(ADs.Dt[,-1])),   pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.Total[,-1])),log2(as.matrix(ADs.At[,-1])),   pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(D5.Total[,-1])),log2(as.matrix(ADs.Dt[,-1])),   pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.Total[,-1]/D5.Total[,-1])),log2(as.matrix(ADs.At[,-1]/ADs.Dt[,-1])), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.At[,-1]/D5.Dt[,-1])),log2(as.matrix(ADs.At[,-1]/ADs.Dt[,-1])), pch=16, col=addTrans("black",10))
    plot(log2(as.matrix(A2.Total[,-1]/D5.Total[,-1])),log2(as.matrix(A2.At[,-1]/D5.Dt[,-1])), pch=16, col=addTrans("black",10))
    
    dev.off()
    
    load('/Volumes/jfw-lab/Projects/Eflen_networks/R-01-countDatasets.RData')
----------------------------
