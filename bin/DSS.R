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
DMS = callDML(dmlTest)
#Write output to hard disk
fdr_cutoff=0.05
index = which(DMS$fdr <= fdr_cutoff)
write.table(DMS[index,],file="output.LYMP.DMR.tsv",row.names=F,quote=F,sep="\t")
