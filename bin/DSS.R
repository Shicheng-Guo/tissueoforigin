setwd("/home/local/MFLDCLIN/guosa/hpc/methylation/PRJNA577192/COV")
library("DSS")

read_methBED <- function(filename)
{
   meth = read.table(file=filename,stringsAsFactors=F,skip=1)
   x = meth[,c(1,2,5)]
   x = cbind(x,as.integer(meth[,5]*meth[,4]))
   colnames(x) = c('chr','pos','N','X');
   return(x)
}

read_bismarkcov <- function(filename){
   meth = read.table(file=filename,stringsAsFactors=F,skip=0)
   N= meth[,4]+meth[,5] 
   X= meth[,4]
   chr=meth[,1]
   pos=meth[,2]
   x= data.frame(chr,pos,N,X)
   return(x)
}

############################################################
#Parameters
fdr_cutoff = 0.05

#Input
filenames_1 = c(
            'LYMP-S2019-0003-1.bedgraph','LYMP-S2019-0004-1.bedgraph'
            )
samplenames_1 = c('LYMP-S2019-0003-1','LYMP-S2019-0004-1')
filenames_2 = c(
            'LYMP-S2019-0003-2.bedgraph',
            'LYMP-S2019-0007-2.bedgraph'
            )
samplenames_2 = c('LYMP-S2019-0003-2','LYMP-S2019-0007-2')
############################################################
#Read input files
input = list()
filenames=c(filenames_1,filenames_2)
for(file in filenames)
{
  input = append(input,list(read_bismarkcov(file)))
}
save(input,file="input.RData")

#Calling differentially methylated sites
BSobj <- makeBSseqData(input,c('LYMP-S2019-0003-1','LYMP-S2019-0004-1','LYMP-S2019-0003-2','LYMP-S2019-0007-2') )
save(BSobj,file="BSobj.RData")
dmlTest <- DMLtest(BSobj, group1=c('LYMP-S2019-0003-1','LYMP-S2019-0004-1'), group2=c('LYMP-S2019-0003-2','LYMP-S2019-0007-2'),smoothing = T, smoothing.span = 200)
save(dmlTest,file="dmlTest.region.RData")
write.table(dmlTest,file="output.LYMP.DMC.tsv",row.names=F,quote=F,sep="\t")

#Write output to hard disk
DMS = callDML(dmlTest)
fdr_cutoff=0.05
index = which(DMS$fdr <= fdr_cutoff)
write.table(DMS[index,],file="output.LYMP.DMR.tsv",row.names=F,quote=F,sep="\t")

# Manhattan plot and qqplot
ManhattanmyDMP<-function(myDMP){
  library("qqman")
  library("Haplin")
  
  colnames(myDMP)[match("chr",colnames(myDMP))]<-"CHR"
  colnames(myDMP)[match("pos",colnames(myDMP))]<-"MAPINFO"
  colnames(myDMP)[match("pval",colnames(myDMP))]<-"P.Value"
     
  SNP=rownames(myDMP)
  CHR=gsub("chr","",myDMP$CHR)
   
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=myDMP$MAPINFO
  P=myDMP$P.Value
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=log(0.05/nrow(na.omit(manhattaninput)),10)
  pdf("manhattan.pdf")
  manhattan(na.omit(manhattaninput),col = c("blue4", "orange3"),ylim = c(0,7),lwd=2, suggestiveline=F,genomewideline=F)
  dev.off()
  pdf("qqplot.pdf")
  pQQ(P, nlabs =length(P), conf = 0.95)
  dev.off()
}

