## Analysis by Guanjing Hu from Dec 1st, 2016
## Adapted from Eflen_Networks110516.r
## This R script is for network analysis using 56 G. hisutum RNA-seq samples from Zhang et al. G. hirsutum var. TM1 genome sequence. Both Gsnap-polycat and Bowtie-rsem mapping was conducted to generate read count tables for the polyploid genome.

ssh hugj2006@bigram.ent.iastate.edu
screen -S zhang
cd /home/hugj2006/jfw-lab/Projects/AD1_mapping
module load R
# Module name: R                          Version: 3.3.1
R
# start R analysis



############## Install WGCNA ####################
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
biocLite("impute")
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("WGCNA")
install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") )


sessionInfo()


############### Step 1a. Prepare polycat network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# cd /home/hugj2006/jfw-lab/Projects/AD1_mapping/SRA/networkAnalysis
# R

getwd()
dir <- "/home/hugj2006/jfw-lab/Projects/AD1_mapping/SRA/fastq_raw/htseqCount/"
fileL<-grep( "N|[.]D|A[.]", list.files(dir), value=TRUE)
# 56 samples x 3(A, D, N)  = 168
remove(allcount)
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t")
    names(x) <- c("gene", file)
    if(!exists("allcount")) {allcount = x}
    else {allcount <- merge(allcount, x, by="gene")}
}
row.names(allcount) <- allcount$gene
allcount<-allcount[,-1]
# Prepare read count tables including only Chr genes
allcount13 <- allcount[grep("Gorai.0",row.names(allcount)),]
dim(allcount13)

# check countable reads
check<-cbind( colSums(allcount ), colSums( allcount[grep("Gorai",row.names(allcount)),] ), t(allcount[grep("__",row.names(allcount)),])  )
colnames(check)[1:2]<-c("Mapped","Counted")
check<-as.data.frame(check)
check$percentageCount <- check$Counted/check$Mapped

check$sample<-gsub("[.].*","",row.names(check) )
aggregate(check$Mapped,list(check$Mapped),sum)
aggregate(check$Counted,list(check$Mapped),sum)

write.table(check,file="1.checkmapping_polycat121916.txt", sep="\t")




##################
At  <- allcount13[,grep('*sort.A.txt',names(allcount13))]
Dt  <- allcount13[,grep('*sort.D.txt',names(allcount13))]
N   <- allcount13[,grep('*sort.N.txt',names(allcount13))]
##################
Aportion <- At/( At+ Dt )
summary(Aportion)   # NA means At=0 and (At+Dt)=0
Aportion[is.na(Aportion)]<-0
AtN<- At + N * Aportion
##################
Dportion <-Dt/( At+ Dt)
summary(Dportion)   # NA means At=0 and (At+Dt)=0
Dportion[is.na(Dportion)]<-0
DtN<- Dt + N * Dportion
##################

# double check files
colSums(At)/colSums(AtN)
colSums(Dt)/colSums(DtN)
(colSums(At) + colSums(Dt) + colSums(N) )/ (colSums(AtN) + colSums(DtN) )    # close to 1



# prepare network datasets
info <- read.table("~/jfw-lab/Projects/AD1_mapping/SRR_Acc_List56.txt",skip=5,sep="\t", header=TRUE)
name <- as.character( info$label[match(gsub(".nsort.*","",names(At)), info$sample)] )
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    rownames(a) <- paste0(rownames(a) , "a")
    rownames(d) <- paste0(rownames(d) , "d")
    ad        <- rbind (a,d)
    return(ad)
}


AtDt         <- prepDatasets(At, Dt)
AtDtN        <- prepDatasets(AtN, DtN   )
# double check content: 74446    56


# rlog transformation
# need DESeq2
library(DESeq2)
# coldata
coldata<-data.frame( sample = name )

rlogTransformation<-function(x)
{
    count <- round( x ) #have to round ADs.ncorrect
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    count.rld <- as.data.frame(assay(rld))
    names(count.rld)<-names(count)
    return(count.rld)
}

AtDt.rld= rlogTransformation(AtDt)
AtDtN.rld= rlogTransformation(AtDtN)

# RPKM transformation
# RPKM stands for Reads Per Kilobase of transcript per Million mapped reads.
# RPKM= 10^9 * number of mapped reads in gene
#             ----------------------------------------------------
#             total reads *  exon length
geneLen<-read.table("~/jfw-lab/Projects/Eflen/eflen_recheck/D5.trueANDtheo.lengths",header=TRUE, sep="\t")
table(geneLen$true>=geneLen$theoretical)
quantile(geneLen$true-geneLen$theoretical)
# prepare RPKM tranformation
rpkm<-function(count, length) {
    # match gene order
    gene <-  gsub(".$","", rownames(count)[1:37223])
    # table(match(gene, lenth$gene)==1:37223)
    len <- length[ match(gene, length$gene),2 ]
    len <- c(len, len)
    librarySize <- colSums(count)
    
    # OR x <- t( t(10^9 * count[,-1]/ geneLen) / librarySize )
    x<-sweep(sweep(count*10^9, 2 ,librarySize, FUN="/"), 1, len, FUN ="/" )
    return(x)
}

AtDt.rpkm <- rpkm(AtDt, geneLen[,c("gene","true" )])
AtDtN.rpkm <- rpkm(AtDtN, geneLen[,c("gene","true" )])

table(is.na(as.matrix(AtDt.rpkm)) )
table(is.infinite(as.matrix(AtDt.rpkm)) )
# ~ gene length cannot be zero, so no NA (=0/0) or Inf (e.g. =4/0)



#save
save(info, AtDt, AtDtN, AtDt.rld, AtDtN.rld, AtDt.rpkm, AtDtN.rpkm, file = "R-01-polycatNetworkDatasets.RData")


############### Step 1b. Prepare RSEM network datasets  ###############
## unix command: R CMD BATCH s1.R
################

#### 07/29/16 Run bowtie and RSEM against AtDt pseudo-transcriptome
#### 12/16 rerun

# relevant R analysis as https://github.com/huguanjing/AD1_RNA-seq/blob/master/seed_eflen_RSEM.md
dir <- "/home/hugj2006/jfw-lab/Projects/AD1_mapping/SRA/fastq_trimmed/rsem/AtDt/"
fileL<- grep("genes.results", list.files(dir),value=TRUE)
# 56 samples
remove(count)
remove(rpkm)
for(file in fileL)
{
    sample<-gsub("_.*","",file)
    x <- read.table(paste0(dir,file), sep="\t", header=TRUE)
    cc <- x[,c("gene_id", "expected_count")]
    names(cc)[2] <- sample
    rr  <- x[,c("gene_id", "FPKM")]
    names(rr)[2] <- sample
    if(!exists("count")) {count=cc} else {count <- merge(count, cc, by="gene_id")}
    if(!exists("rpkm")) {rpkm = rr} else {rpkm  <- merge(rpkm,  rr, by="gene_id")}
}


AtDt<-count[,-1]
row.names(AtDt) <- gsub(".D","d", gsub(".A","a", count$gene_id))

AtDt.rpkm<-rpkm[,-1]
row.names(AtDt.rpkm) <- gsub(".D","d", gsub(".A","a", rpkm$gene_id))

names(AtDt) ==names(AtDt.rpkm)

info <- read.table("~/jfw-lab/Projects/AD1_mapping/SRR_Acc_List56.txt",skip=5,sep="\t", header=TRUE)
name <- as.character( info$label[match(names(AtDt), info$sample)] )
names(AtDt) =name
names(AtDt.rpkm) =name

# rlog transformation
# need DESeq2
library(DESeq2)
# coldata
coldata<-data.frame( sample = name )

rlogTransformation<-function(x)
{
    count <- round( x ) #have to round ADs.ncorrect
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    count.rld <- as.data.frame(assay(rld))
    names(count.rld)<-names(count)
    return(count.rld)
}

AtDt.rld= rlogTransformation(AtDt)


#save
save(info, AtDt, AtDt.rld, AtDt.rpkm, file = "R-01-rsemNetworkDatasets.RData")

pdf("s1.checkSamples.rsem.pdf")
boxplot(log2(AtDt), las=2, main = "log2(AtDt)")
boxplot(AtDt.rld, las=2, main = "AtDt.rld")
boxplot(log2(AtDt.rpkm ), las=2, main="log2(AtDt.rpkm)")
dev.off()



############### Step 2. Check outliers  ###############
## unix command: R CMD BATCH s2.R
################
library(WGCNA)
library(ggplot2)
library(flashClust)

checkOutliers <- function(datExpr, title)
{
    # sample network based on squared Euclidean distance
    A=adjacency(datExpr,type="distance")
    # this calculates the whole network connectivity
    k=as.numeric(apply(A,2,sum))-1
    # standardized connectivity
    Z.k=scale(k)
    # Designate samples as outlying
    # if their Z.k value is below the threshold
    thresholdZ.k1 <- -5 # often -2.5
    thresholdZ.k2 <- -2.5 # often -2.5
    # the color vector indicates outlyingness (red)
    outlierColor1=ifelse(Z.k<thresholdZ.k1,"red","black")
    outlierColor2=ifelse(Z.k<thresholdZ.k2,"red","black")
    datColors=data.frame(outlierColor1,outlierColor2)
    # Plot the sample dendrogram and the colors underneath.
    sampleTree = flashClust(as.dist(1-A), method = "average")
    plotDendroAndColors(sampleTree,groupLabels=names(datColors),
    colors=datColors,main=title)
}

polycat <- load("R-01-polycatNetworkDatasets.RData")
# info, AtDt, AtDtN, AtDt.rld, AtDtN.rld, AtDt.rpkm, AtDtN.rpkm
pdf("s2.checkOutlier.polycat.pdf")
boxplot(colSums(AtDt)~info$tissue, las=2, main="library Size - AtDt")
boxplot(colSums(AtDtN)~info$tissue, las=2, main="library Size - AtDtN")
for(i in setdiff(polycat, "info")){
    net.dat <- get(i)
    # remember to log2 for not rld sets
    if(!grepl("rld",i)){net.dat<-log2(net.dat+1)}
    pca = prcomp(t(net.dat))
    dat = as.data.frame(pca$x)
    var = summary(pca)$importance[2,]
    dat$group= info$tissue[match(rownames(dat),info$label)]
    print( ggplot(aes(PC1, PC2, color=group),data=dat) + geom_point() + geom_text(aes(label=rownames(dat)),hjust=0.5, vjust=0) + ggtitle(i) + xlab(paste0("PC1: ",var[1]*100,"% variance")) + ylab(paste0("PC2: ",var[2]*100,"% variance"))
    )
    
    checkOutliers(net.dat,i)
}
dev.off()


rsem <- load("R-01-rsemNetworkDatasets.RData")
# info, AtDt, AtDtN, AtDt.rld, AtDtN.rld, AtDt.rpkm, AtDtN.rpkm
pdf("s2.checkOutlier.rsem.pdf")
boxplot(colSums(AtDt)~info$tissue, las=2, main="library Size - AtDt")
for(i in setdiff(rsem, "info")){
    net.dat <- get(i)
    # remember to log2 for not rld sets
    if(!grepl("rld",i)){net.dat<-log2(net.dat+1)}
    pca = prcomp(t(net.dat))
    dat = as.data.frame(pca$x)
    dat$group= info$tissue[match(rownames(dat),info$label)]
    print( ggplot(aes(PC1, PC2, color=group),data=dat) + geom_point() + geom_text(aes(label=rownames(dat)),hjust=0.5, vjust=0) + ggtitle(i) + xlab(paste0("PC1: ",var[1]*100,"% variance")) + ylab(paste0("PC2: ",var[2]*100,"% variance"))
    )
    
    checkOutliers(net.dat,i)
}
dev.off()



# RSEM total library size indicate an outlier CK0
m<-mean(colSums(AtDt) )
colSums(AtDt)/m

# RSEM - remove CK0
# polyCat - keep CK0
# other samples are kept, while make sure to use bicor for network construction for outliers


############### Step 3. Prep for WGCNA  ###############
## unix command: R CMD BATCH s2.R
################

library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# enableWGCNAThreads(nThreads=6)

nSets=6

# remove outliers
rsem    <- load("R-01-rsemNetworkDatasets.RData") # "info"      "AtDt"      "AtDt.rld"  "AtDt.rpkm"
AtDt.rld <- AtDt.rld[,setdiff(names(AtDt.rld), c("CK0", "stamen"))]
AtDt.rpkm <- AtDt.rpkm[,setdiff(names(AtDt.rpkm), c("CK0"))]
# Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = t(AtDt.rld) )
multiExpr[[2]] = list(data = t(log2(AtDt.rpkm+1)) )
polycat    <- load("R-01-polycatNetworkDatasets.RData") # "info"     "info"       "AtDt"       "AtDtN"    "AtDt.rpkm"  "AtDtN.rpkm"
AtDt.rld <- AtDt.rld[,setdiff(names(AtDt.rld), c( "stamen"))]
AtDtN.rld <- AtDtN.rld[,setdiff(names(AtDtN.rld), c( "stamen"))]
multiExpr[[3]] = list(data = t(AtDt.rld) )
multiExpr[[4]] = list(data = t(log2(AtDt.rpkm+1)))
multiExpr[[5]] = list(data = t(AtDtN.rld))
multiExpr[[6]] = list(data = t(log2(AtDtN.rpkm+1)))
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 6
# $nGenes 74446
# $nSamples 54 55 55 56 55 56
# $structureOK TRUE
shortLabels = c("rsem.rld","rsem.rpkm","polycat.rld","polycat.rpkm","polycat.rld.Npar","polycat.rpkm.Npar")


# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 8616 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Excluding genes from the calculation due to too many missing samples or zero variance.
# Update exprSize
checkSets(multiExpr)
# $nSets 6
# $nGenes 69554
# $nSamples 54 54 55 55 55 55
# $structureOK TRUE
save(multiExpr,nSets, shortLabels, file = "R-03-dataInput.RData")


pdf(file = "s3.SampleClusteringS.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", shortLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;


#######################################
#  Choosing the soft-thresholding power: analysis of network topology

rm(list=ls())
rdatafiles<-c("R-03-dataInput.RData")
corFunction <-  "bicor"

for(file in rdatafiles)
{
     print(file)
     load(file)
     tag<-gsub(".*Input[.]|[.].*","",file)

# Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
# types<-c("unsigned", "signed", "signed hybrid")
    type <- "signed"
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type, blockSize=10000, corFnc= corFunction )[[2]])      }
    collectGarbage()
    
    # Plot the results:
    colors=brewer.pal(5,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
    }
    
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    # sizeGrWindow(8, 6)
    pdf(paste0("s3.ChooseSoftThresholdPower_", tag,".pdf") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    dev.off()
    assign(paste0("powerTables.",tag),powerTables)
}

# blockSize(65830, rectangularBlocks = TRUE) can calculates suitable block size for Biocrunch
# I just use 10000 here


Powers <- unlist(lapply(powerTables, function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.8 & x$slope<0][1]} ) )
#   7 20  7 22  6 22
names(Powers)<-shortLabels
# "rsem.rld"          "rsem.rpkm"         "polycat.rld"
# "polycat.rpkm"      "polycat.rld.Npar"  "polycat.rpkm.Npar"
Powers
# ADs
### polycat_rld: 24 28 24 24  =>28
### polycat_rpkm: 24 24 20 24 24 20  =>24
### rsem_rld: 30 28 26  => 30
### rsem_rpkm: 24,22,28 => 28
save( powerTables,Powers, file = "R-03-chooseSoftThreshold.Rdata")


############### Step X.  Differential coexpression analysis  ###############
## nohup R CMD BATCH s4.networkConstruction.R &
################
library(DiffCorr)

rm(list=ls())
rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

rp <- rbind(c("Set", "X", "Y", "DCpairs" ) )
for(file in rdatafiles)
{
    print(file)
    load(file) # # multiExpr,nSets, setLabels, shortLabels
    tag<-gsub(".*Input[.]|[.].*","",file)
    
    # Analysis of differntial coexpression gene pairs
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2))
    {
        print(paste0("=== Comparing ", shortLabels[ pwset[1,i] ], " vs ", shortLabels[ pwset[2,i] ]))
        X<-t(multiExpr[[ pwset[1,i] ]]$data)
        Y<-t(multiExpr[[ pwset[2,i] ]]$data)
        outfile <- paste0("DC/",tag, ".", shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".res.txt" )
        
        # DiffCorr code cannot handle too big a correlation matrix
        # if( nrow(X) < 40000)
        # comp.2.cc.fdr(output.file=outfile, X, Y, threshold=0.05)
        res <- comp.2.big.cc(output.file=outfile,X,Y, 2000)
        rp <- rbind(rp, c( tag, shortLabels[ pwset[1,i] ], shortLabels[ pwset[2,i] ], nrow(res)   )  )
    }
    
}
DC<-rp
DC<-as.data.frame(rp[-1,])
names(DC)<-rp[1,]

n<-c()
for(file in rdatafiles)
{
    print(file)
    load(file) # # multiExpr,nSets, setLabels, shortLabels
    tag<-gsub(".*Input[.]|[.].*","",file)
    nGenes <- ncol(multiExpr[[1]]$data)
    
    # Analysis of differntial coexpression gene pairs
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2))
    {
        n<- c(n, nGenes)
    }
    
}
DC$nGenes <- n
DC$DCpercent <- as.numeric(as.character(DC$DCpairs))/(DC$nGenes*(DC$nGenes-1)/2)
write.table(DC, file<-"sX.differentialCoexpression.txt", row.names=FALSE, sep="\t")

DC[DC$X=="A2D5"&DC$Y=="ADs",]
#           Set    X   Y  DCpairs   DCpercent nGenes
#     rsem_rpkm A2D5 ADs 11504289 0.005280842  66008
#      rsem_rld A2D5 ADs  8665348 0.004082008  65159
#  polycat_rpkm A2D5 ADs  3277727 0.001589234  64226
#   polycat_rld A2D5 ADs  3059846 0.001481285  64276


##### skipped


############### Step 4.  Network Construction  ###############
## nohup R CMD BATCH s4.networkConstruction.R &
################

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

ptm <- proc.time()

load("R-03-chooseSoftThreshold.Rdata")
Powers

rdatafiles<-c("R-03-dataInput.RData")

# rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.polycat_rpkm.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file)
    # compr<-gsub("R-02-dataInput.|.RData","",file)
    # powerEach = Powers[compr]
    # print(paste0(file,", construct network with b = ",powerEach))
    print(shortLabels)
    print(checkSets(multiExpr)$nGenes)
    
    
    ###### calculate individual TOMs
    print("###### Calculate individual TOMs:")
    iTOMs = blockwiseIndividualTOMs(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "bicor",
    # Adjacency and topology overlap function options
    #power = powerEach,
    power=max(Powers),
    networkType = "signed", TOMType = "signed",
    # Save individual TOMs?
    saveTOMs = TRUE,
    individualTOMFileNames = paste0(compr,".iTOM-%N-block.%b.RData")  )
    
    ###### calculate consensus modules
    print("###### Construct consensus networks:")
    cnet = blockwiseConsensusModules(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "bicor",
    # Adjacency and topology overlap function options
    #power = powerEach,
    power=max(Powers),
    networkType = "signed", TOMType = "signed",
    # load previous TOMs
    individualTOMInfo = iTOMs,
    # Saving the consensus TOM
    saveConsensusTOMs = TRUE,
    consensusTOMFileNames = paste0(compr,".cTOM-Block%b.RData"),
    # Basic tree cut options
    deepSplit = 2,  #default, known to reasonable
    minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
    pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
    # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
    mergeCutHeight = 0.25,
    # others
    reassignThreshold = 0,
    numericLabels = TRUE,
    verbose = 3)
    
    ###### calculate individual modules
    print("###### Construct individual networks:")
    # stupid blockwiseModules only load TOM rdata file with "TOM", not "tomDS"
    tomFiles<-grep(paste0(compr,".iTOM"),list.files(), value=TRUE)
    for(fl in tomFiles)
    {
        load(fl)
        TOM<-tomDS
        save(TOM, file=fl)
    }
    rm(TOM)
    collectGarbage()
    for(i in 1:nSets)
    {
        inet = blockwiseModules(
        # Input data
        multiExpr[[i]]$data,
        # Data checking options
        checkMissingData = TRUE,
        # Options for splitting data into blocks
        blocks =  iTOMs$blocks,
        #randomSeed = 12345,
        #maxBlockSize =  500,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
        # Network construction arguments: correlation options, use bicor instead of default pearson
        corType = "bicor",
        # Adjacency and topology overlap function options
        #power = powerEach,
        power=max(Powers),
        networkType = "signed", TOMType = "signed",
        # load previous TOMs
        loadTOM = TRUE,
        saveTOMFileBase = paste0(compr,".iTOM-",i),
        # Basic tree cut options
        deepSplit = 2,  #default, known to reasonable
        minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
        pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
        # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
        mergeCutHeight = 0.25,
        # others
        reassignThreshold = 0,
        numericLabels = TRUE,
        verbose = 3)
        
        assign(paste0("inet",i), inet)
        
    }
    
    save(list=c( "iTOMs", "cnet", grep("inet.",ls(), value=TRUE)), file = paste0("R-04-buildNetwork.", compr,".RData"))
    collectGarbage()

}

proc.time() - ptm

# With maxBlockSize = 20000, it took roughly a week



############### Step 5.  General network topology analysis  ###############
## nohup R CMD BATCH s5.networkTopology.R &
################

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

plotRefModulePres<-function(mp, file= "")
{
    pdf(file)
    for(set in 2:nSets ){
        # specify the reference and the test networks
        ref=1; test = set
        Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
        Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
        # Look at the observed preservation statistics
        # Obs.PreservationStats
        # Z statistics from the permutation test analysis
        # Z.PreservationStats
        
        # Let us now visualize the data.
        modIDs = rownames(Obs.PreservationStats)
        modColors=labels2colors(order(as.numeric(modIDs) )-1 )
        moduleSize = Obs.PreservationStats$moduleSize
        # we will omit the grey module (background genes)
        # and the gold module (random sample of genes)
        selectModules = !(modColors %in% c("grey", "gold"))
        # Text labels for points
        point.label = modIDs[selectModules]
        # Composite preservation statistics
        medianRank=Obs.PreservationStats$medianRank.pres
        Zsummary=Z.PreservationStats$Zsummary.pres
        
        par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
        # plot medianRank versus module size
        plot(moduleSize[selectModules],medianRank[selectModules],col=1, bg=modColors[selectModules],
        pch = 21,main=paste("medianRank -",shortLabels[ref], "vs",shortLabels[test]),
        cex = 2, ylab ="medianRank",xlab="Module size", log="x")
        labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
        
        # plot Zsummary versus module size
        plot(moduleSize[selectModules],Zsummary[selectModules], col = 1, bg=modColors[selectModules],pch = 21,
        main=paste("Zsummary -",shortLabels[ref], "vs",shortLabels[test]),
        cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
        labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
        # Add threshold lines for Zsummary
        abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
    }
    dev.off()
}

load("R-03-chooseSoftThreshold.Rdata")
load("R-04-buildNetwork.Gh.RData")


Powers

rdatafiles<-c("R-03-dataInput.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file) # "multiExpr"   "nSets"   "shortLabels"
    # compr<-gsub("R-03-dataInput.|.RData","",file)
    compr="Gh"
    # powerEach = Powers[compr]
    powerEach=max(Powers)
    print(paste0(file,", construct network with b = ",powerEach))
    # print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    load(paste0("R-04-buildNetwork.", compr,".RData")) ->Names
    print(Names)  # "iTOMs" "cnet"  "inet1" "inet2" "inet3"  "inet4" "inet5" "inet6"
    
    ###  1.Comparing basic topological networks parameters
    ######################################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        gg<-net$goodGenes    #get the good genes descision
        
        ###  Get basic topological networks parameters
        print(paste("Start to extract network parameter for:",genome))
        # network concept
        adj <- adjacency(subDat, power = powerEach,type = "signed")   # calculate adjacency table
        # fundmentalNetworkConcept taking forever
        # net1C <- fundamentalNetworkConcepts(adj, GS = NULL)           # computes fundamental network concepts
        Size = dim(adj)[1]
        Connectivity = apply(adj, 2, sum)-1
        Density = sum(Connectivity)/(Size * (Size - 1))
        Centralization = Size * (max(Connectivity) - mean(Connectivity))/((Size -
        1) * (Size - 2))
        Heterogeneity = sqrt(Size * sum(Connectivity^2)/sum(Connectivity)^2 -
        1)
        # clustering coef takes forever to computate, skip
        # ClusterCoef = .ClusterCoef.fun(adj)
        fMAR = function(v) sum(v^2)/sum(v)
        MAR = apply(adj, 1, fMAR)
        ScaledConnectivity = Connectivity/max(Connectivity, na.rm = T)
        output = list(Connectivity = Connectivity, ScaledConnectivity = ScaledConnectivity,  MAR = MAR, Density = Density, Centralization = Centralization, Heterogeneity = Heterogeneity)
        param <- unlist(lapply(output,mean))
        
        ### Get module and dupllication related topology
        # how many modules detected??
        param["nModules"]<-length(unique(net$colors))
        # interal and crossing connections
        dt<-grep("d$",rownames(adj))
        size.dt<-length(dt)
        at<-grep("a$",rownames(adj))
        size.at<-length(at)
        param["at.density"] <- (sum(adj[at,at])-size.at)/(size.at*(size.at-1))
        param["dt.density"] <- (sum(adj[dt,dt])-size.dt)/(size.dt*(size.dt-1))
        param["ad.density"] <- mean(adj[at,dt])
        ### write results table
        if(!exists("Rtable"))
        {
            Rtable <- data.frame(param)
            names(Rtable)<-genome
        }else
        {Rtable[,genome] <- param }
        
        
        assign(paste0(genome,"Adj"),adj)
        assign(paste0(genome,"Concepts"),output)
        
    }
    #  Save basic topological networks parameters
    rownames(Rtable)[1:3]<- paste0("mean",rownames(Rtable)[1:3])
    write.table(Rtable, file = paste0("s5.parameters.",compr,".txt"),sep="\t")
    save( list=c( grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
   
    
    ###  2.Correlating node-specific network properties
    ###################################################
    pwset<-combn(nSets,2)
    
    # pdf(paste0("s5.correlatingParameters.",compr,".pdf") )
    # par(mfrow=c(2,2))
    # for( j in 1:3) # "Connectivity"       "ScaledConnectivity" "MAR"
    # {
    #    for(i in 1:ncol(pwset) )
    #    {   xnet <- get(paste0(shortLabels[ pwset[1,i]  ],"Concepts"))
    #        ynet <- get(paste0(shortLabels[ pwset[2,i]  ],"Concepts"))
    #        pp <- names(xnet)[j]
    #        corr <- cor.test(xnet[[j]], ynet[[j]])
    #        plot (xnet[[j]], ynet[[j]], main = paste(pp, ": cor=", round(corr$estimate,2), ", p=", corr$p.value, sep=""), xlab= shortLabels[pwset[1,i]], ylab = shortLabels[pwset[2,i]], type ="p", pch=16, col = rgb(0, 0, 0, 0.2) )
    #    }
    # }
    # dev.off()
    
    df_Connectivity <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$Connectivity, get(paste0(shortLabels[ 2  ],"Concepts"))$Connectivity)
    df_ScaledConnectivity <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$ScaledConnectivity, get(paste0(shortLabels[ 2  ],"Concepts"))$ScaledConnectivity)
    df_MAR <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$MAR, get(paste0(shortLabels[ 2  ],"Concepts"))$MAR)
    df_expr <- cbind(  as.numeric(as.matrix(multiExpr[[1]]$data)),  as.numeric(as.matrix(multiExpr[[2]]$data))  )
    for(i in 3:nSets)
    {
        df_Connectivity <- cbind( df_Connectivity, get(paste0(shortLabels[ i  ],"Concepts"))$Connectivity)
        df_ScaledConnectivity <- cbind( df_ScaledConnectivity , get(paste0(shortLabels[ i  ],"Concepts"))$ScaledConnectivity)
        df_MAR <- cbind( df_MAR, get(paste0(shortLabels[ i  ],"Concepts"))$MAR)
        df_expr <- cbind(  df_expr,  as.numeric(as.matrix(multiExpr[[i]]$data))  )
    }
    colnames(df_Connectivity) <-shortLabels
    colnames(df_ScaledConnectivity) <-shortLabels
    colnames(df_MAR) <-shortLabels
    colnames(df_expr) <-shortLabels
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }
    
    pdf(paste0("s5.correlatingParameters.",compr,".k.pdf") )
    pairs(df_Connectivity, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Connectivity Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".sk.pdf") )
    pairs(df_ScaledConnectivity, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Scaled Connectivity Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".mar.pdf") )
    pairs(df_MAR, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene MAR Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".expr.pdf") )
    pairs(df_expr, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Expression Scatterplot Matrix")
    dev.off()



    ###  3.WGCNA provided methods for comparison
    ############################################
    # get all anova p and put together
    info <- read.table("~/jfw-lab/Projects/AD1_mapping/SRR_Acc_List56.txt",skip=5,sep="\t", header=TRUE)
    for(i in 1:nSets)
    {
        net<-get(paste0("inet",i))
        shortLabels[i]
        
        pdf("s1.sample_dendrogram_and_trait_heatmap.rlog.pdf")
        # calculate the cluster tree using flahsClust or hclust
        sampleTree = flashClust(as.dist(1-A), method = "average")
        # Convert traits to a color representation:
        # where red indicates high values
        traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
        dimnames(traitColors)[[2]]=paste(names(datTraits),"C",sep="")
        datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, traitColors)
        
        # Plot the sample dendrogram and the colors underneath.
        plotDendroAndColors(net$dendrograms,groupLabels=names(datColors),
        colors=datColors,main="Sample dendrogram and trait heatmap")
        dev.off()
    }

    anovaP<-list()
    for(set in 1:nSets)
    {
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        MEs<-net$MEs
        sample <- rownames(multiExpr[[set]]$data)
        tissue<-info$tissue[match(sample,info$label)]
        pval<-apply(MEs,2,function(x){round(anova(aov(x~tissue) )$"Pr(>F)"[1],4)})
        pval<-as.data.frame(pval)
        pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
        pval$numeric<-as.numeric(substring(rownames(pval),3) )
        pval<-pval[order(pval$numeric),]
        pval$symbol[1]<-" "  # ME0 always meaningless
        anovaP[[set]]<-pval
    }
    names(anovaP)<-shortLabels
    
    # Module correspondence test
    pdf(paste0("s5.moduleCorrepondence.",compr,".pdf"),width=10,height=7)
    par(mfrow=c(1,1));
    par(cex = 1.0);
    par(mar=c(8, 10.4, 2.7, 1)+0.3);
    # loop pairwise comparison
    for(i in 1:ncol(pwset))
    {
        print(paste0("Plot for pairwise set ",i))
        colnet<-get(paste0("inet",pwset[1,i]) )
        rownet<-get(paste0("inet",pwset[2,i]) )
        coln<-ncol(colnet$MEs )  # number of MEs
        rown<-ncol(rownet$MEs )
        # color list of MEs in the color of decreasing numbers of memebers
        colModule  <- labels2colors(as.numeric(names(table(colnet$colors)) ))
        rowModule  <- labels2colors(as.numeric(names(table(rownet$colors)) ))
        # colors for each gene
        colColors  <- labels2colors(colnet$colors )
        rowColors  <- labels2colors(rownet$colors )   # colors for each gene
        # anova significance sybol
        colP  <- anovaP[[pwset[1,i]]]$symbol
        rowP  <- anovaP[[pwset[2,i]]]$symbol
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # Execute all pairwaise comparisons
        for (rmod in 1:rown)
        for (cmod in 1:coln)
        {
            rMembers = (rowColors == rowModule[rmod] );
            cMembers = (colColors == colModule[cmod] );
            pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
            CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
        }
        
        # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        rModTotals = apply(CountTbl, 1, sum)
        cModTotals = apply(CountTbl, 2, sum)
        # Use function labeledHeatmap to produce the color-coded table with all the trimmings
        labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
        xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule),
        xSymbols = paste0(shortLabels[pwset[1,i]], ".",colModule, "_", cModTotals, colP),
        ySymbols = paste0(shortLabels[pwset[2,i]], ".",rowModule, "_", rModTotals, rowP),
        textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of modules",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
        
        # plotting only singficant modules
        labeledHeatmap( Matrix = pTable[rowP=="*",colP=="*"], colorLabels = TRUE,
        xLabels = paste(" ", colModule[colP=="*"]), yLabels = paste(" ", rowModule[rowP=="*"]),
        xSymbols = paste0(shortLabels[pwset[1,i]], ".", colModule[colP=="*"], ": ", cModTotals[colP=="*"], colP[colP=="*"]),
        ySymbols = paste0(shortLabels[pwset[2,i]], ".", rowModule[rowP=="*"], ": ", rModTotals[rowP=="*"], rowP[rowP=="*"]),
        textMatrix = CountTbl[rowP=="*",colP=="*"], colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of modules, significant only",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
    }
    dev.off()
    
    # plot module eigengenes
    source("~/jfw-lab/scripts_old/scripts_HGJ/summarySE.r")
    source("~/jfw-lab/scripts_old/scripts_HGJ/multiplot.r")
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
        # Displaying module heatmap and the eigengene
        # sizeGrWindow(8,7);
        MEs<-net$MEs
        plots <- list()  # new empty list
        ss = as.factor(info$tissue[match(rownames(subDat),info$label)])
        Nmodules <- ncol(MEs)
        pdf(paste0("s5.moduleEigengenes.",genome,".pdf"),width=10,height=7)
        for(me in 0:(Nmodules-1)) {
            which.module=paste("ME",me,sep="")
            module.color=labels2colors(me)
            #heatmap
            par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
            plotMat(t(scale(subDat[,net$colors==me ]) ),
            nrgcols=30,rlabels=T,rcols=module.color,
            main=paste(which.module, module.color, sep=": "), cex.main=2)
            #barplot
            par(mar=c(5, 4.2, 0, 0.7))
            barplot(MEs[,which.module], col=module.color, main="", cex.main=1, las=2,
            ylab="eigengene expression",xlab="tisse types", names.arg=paste(rownames(subDat),ss,sep=".") )
            #line, anova
            df<-data.frame(ME=MEs[,which.module], ss, module = which.module )
            fit<-aov(ME~ss,df)
            dfc<-summarySE(df, measurevar="ME", groupvars=c("ss", "module"))
            plots[[me+1]]<- ggplot(dfc, aes(x=ss, y=ME, fill = ss)) +
            geom_bar(stat="identity",position=position_dodge(), color="black", size=0.3) +
            geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.3,position=position_dodge(0.9)) +
            ggtitle(paste(which.module," ",module.color,", P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
            theme_bw() +
            theme(plot.title=element_text( size=11),legend.position = "none", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

            }
        for(page in 1:ceiling(Nmodules/9))
        {
            if(Nmodules>(9*page))
            {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
            else
            {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        }
        plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
        dev.off()# end the big plot
    }

    
    ###  4.Homoeolog focused analysis
    ############################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        inPairs<-intersect( genes[subgenome =="a"], genes[subgenome=="d"] )
        homoeo <- c(length( inPairs ), length(setdiff( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="d"], genes[subgenome=="a"] ) ))
        names(homoeo)<-c( "ADinNet", "AinNet", "DinNet" )
        probes <- colnames(subDat )
        genes <- gsub(".$","",probes)
        subgenome <- gsub(".*00","",probes)
        homoeo <- c(length(intersect( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="d"], genes[subgenome=="a"] ) ))
        names(homoeo)<-c( "ADinNet", "AinNet", "DinNet" )
        colors <-net$colors
        inModuleNo<-aggregate(net$colors, by=list(genes), function(x)length(unique(x)))
        rownames(inModuleNo)<-inModuleNo$Group.1
        inModuleNo<-inModuleNo[inPairs,]
        inSameModule <- length(which(inModuleNo$x==1))
        homoeo["ADinSameModule"]<-inSameModule
        
        dominance<- as.data.frame(unclass(xtabs(~net$colors+subgenome) ) )
        dominance$chisq.p<-apply(dominance,1,function(x)chisq.test(x)$p.value)
        dominance$dominant <- apply(dominance,1,function(x)ifelse(x[3]>0.05,"-",ifelse(x[1]>x[2],"a","d")) )
        homoeo["a"] <- length(which(dominance$dominant=="a"))
        homoeo["b"] <- length(which(dominance$dominant=="b"))
        homoeo["-"] <- length(which(dominance$dominant=="-"))
        if(!exists("Rtableh"))
        { Rtableh<-data.frame(cbind(homoeo))
            names(Rtableh)<-genome
        }else {Rtableh[,genome]<-homoeo}

        assign(paste0(genome,".dominance"),dominance)
    }
    write.table(Rtableh, file = paste0("s5.homoeolog.",compr,".txt"),sep="\t")
    
    # save everything, AGAIN
    save( list=c("anovaP",  grep(".dominance",ls(), value=TRUE), grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
}

proc.time() - ptm



############### Step END.  Extra analysis  ###############
## nohup R CMD BATCH s.networkTopology.R &
################

# Plot grapgh parameters with PCA
library(ggplot2)
pdf("s6.graphParametersPCA.pdf")
f<-read.table(file = "s5.parameters.Gh.txt", header=TRUE, sep="\t" )
pca = prcomp(t(f), scale.=T)
# scale = 0 ensures that arrows are scaled to represent the loadings.
# focus on the extreme ends (top, bottom, left, right), e.g. first principal component corresponds to a measure of left and right arrows.
# For exact measure of a variable in a component, pca$rotation
biplot(pca, scale = 0 , main=i)

dat = as.data.frame(pca$x)
dat$group= rownames(dat)
portion.var <- as.numeric( summary(pca)$importance[2,] )
print(
ggplot(aes(PC1, PC2, color=group),data=dat) +
geom_point() + ggtitle(i) +
xlab(paste0("PC1: ",portion.var[1]*100,"% variance")) +
ylab(paste0("PC2: ",portion.var[2]*100,"% variance"))
)


# so empir corrected profiles are more similiar to uncorrected rpkm than theoretical eflen corrected profiles??
dev.off()
###############bookmark


?
Are duplicated gene pairs more likely to be found in the same or different modules? Using co-residence in the same module to indicate conservation in co-expression patterns, in my opinion, should be more informative and robust than comparing whole network neighborhoods between duplicate genes. And the closeness of co-residence could be quantified within module to measure the level of conservation. Gene pairs found in different modules would suggest diverged co-expression patterns, which may be quantified through inter-modular relationships.

Whether and how are conserved and diverged co-expression patterns related to biological function, gene family, metabolic pathways, cellular compartmentation, etc? These information could be layered on modular structure for testing, and because modules restrict tests to genes with evidence association, it should be more accurate than simply testing, for example, all genes annotated in one pathway while half not co-expressed at all.

With respect to subgenome asymmetry, are duplicated genes enriched in or occupying different sub-network space? This was already tested by searching "dominant modules" in both maize and wheat papers showing different results. When statistical equivalent numbers of maize1 and maize2 were found within a module, subgenome "symmetry" was called, but what if maize1 genes and maize2 genes occupy very different space within the module? This can tested for better resolution.



