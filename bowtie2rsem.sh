#!usr/bin/bash
echo ''
echo ''
echo "Starting Job on "
stampStart=`date`
echo $stampStart 

module load bowtie2/2.2.6
module load rsem/1.2.22

SUFFIX=.cut.fq.gz

for j in $( ls eflen_seed/ | grep 'cut.fq.gz' |sed "/D5-20-R2/d" );	do
   
echo ''
echo "==========Running bowtie2-RSEM for"
echo $j
echo ""
rsem-calculate-expression -p 4 --bowtie2 --time eflen_seed/$j ~/jfw-lab/GenomicResources/pseudogenomes/A2D5 bowtie2RSEM/${j%$SUFFIX} >bowtie2RSEM/${j%$SUFFIX}.log 2>&1

done

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd
