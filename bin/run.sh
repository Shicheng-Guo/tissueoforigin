## Tissue-of-orign
# fastq2cov
# cov2bedgraph
# cov2bw

## checking alignment status 
cd /gpfs/home/guosa/hpc/methylation/PRJNA577192

## alignment with bismarker were conducted in GPFS disk not bigdata (bigdata is slow for IO)

cd /gpfs/home/guosa/hpc/methylation/PRJNA577192
cp /gpfs/home/guosa/PRJNA577192/methyfreq/*.cov.gz ./
cp /gpfs/home/guosa/PRJNA577192/methyfreq/*.bedGraph.gz ./
gunzip *.cov.gz
mv *.cov COV
SRRmerge.pl contig.txt

## copy aligned cov and bedgraph to bigdata for further analysis
cd /gpfs/home/guosa/hpc/methylation/PRJNA577192/COV
for i in `ls *.bedgraph`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/methylation/PRJNA577192/COV >> $i.job
echo awk \'{print \$1,\$2-1,\$3,\$4}\' $i OFS=\"\\t\" \>$i.temp >>$i.job
echo wigToBigWig $i.temp ~/hpc/db/hg19/hg19.chrom.sizes $i.bedgraph.bw >> $i.job
echo bigWigAverageOverBed $i.bedgraph.bw ~/hpc/db/hg19/BUR.GRCH37.bed $i.BUR.tab >> $i.job
echo bigWigAverageOverBed $i.bedgraph.bw ~/hpc/db/hg19/CpGISF.hg19.bed $i.CPGI.tab >> $i.job
qsub $i.job
done
wget https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/medip/tab2matrix.pl -O tab2matrix.pl
perl tab2matrix.pl BUR > BUR.matrix.txt
perl tab2matrix.pl CPGI > CPGISF.matrix.txt

cd ~/hpc/methylation/PRJNA577192/COV
library("gplots")

for(i in 1:8){
data<-read.table("CPGISF.matrix.txt",head=T)
id1<-unlist(lapply(strsplit(colnames(data),".S2019"),function(x) x[1]))
id2<-unlist(lapply(strsplit(colnames(data),"[.]"),function(x) x[4]))
colnames(data)=paste(id1,id2,sep="-")
data<-data[-which(unlist(apply(data,1,function(x) sum(is.na(x))))>i),]
png(paste("heatmap",i,"png",sep="."))
heatmap.2(data.matrix(data), main = "Secretome Profile", trace ="none", na.color="gray", margins = c(5,15), cexRow = 0.7)
dev.off()
}



