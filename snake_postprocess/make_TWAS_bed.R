#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: make_TWAS_bed.R <gene> <trait> <tissue> <data_source> \n", call.=FALSE)
}
library(data.table)
library(tidyverse)

#Need to make bed file:
  #Get chr
  #Get stop/start of each SNP
  #Name: rsid
  #No score
  #No strand

trait=args[1]
tissue=args[2]
gene=args[3]
data_source=args[4]


cmd=paste("awk '$2==\"'",gene,"'\" {print $3}' /project2/nobrega/grace/expression/",trait,"/",tissue,"/results/",data_source,".all.dat.top",sep='')
chrom=system(cmd,intern=TRUE)
rsid_locs<-fread(paste("/project2/nobrega/grace/genos/",data_source,"/plink/chr",chrom,"_rsids.bim",sep=""))
dir.create(paste("/project2/nobrega/grace/expression/",trait,"/",tissue,"/results/posthoc",sep=""), showWarnings = FALSE)

get_rsid_vals<-function(rsid) {
  start<-rsid_locs$V4[rsid_locs$V2 == rsid][1]
  stop<-start+1
  return(c(start,stop,rsid))
}

make_TWAS_bed<-function(gene) {
  chr=paste("chr",chrom,sep="")
  load(paste("/project2/nobrega/grace/expression/",trait,"/",tissue,"/WEIGHTS/",gene,".wgt.RDat",sep=""))
  model=which.max(cv.performance[1,])[[1]][1]
  
  if (model==1 | model==2 ) { #Top1 and BLUP assigns a weight to every variant, so make an arbitrary threshold of top 10% of variants to visualize
    vals<-wgt.matrix[which(wgt.matrix[,model] > 0),] 
    threshold<-quantile(abs(vals[,model]),probs=0.9)[[1]][1] 
    rsids<-names(which(abs(vals[,model]) > threshold))
  } else {
    rsids<-names(which(abs(wgt.matrix[,model]) > 0)) 
  }
  chr<-rep(paste("chr",chrom,sep=""),length(rsids))
  score<-rep(".",length(rsids))
  strand<-rep(".",length(rsids))
  locs<-t(as.data.frame(lapply(rsids,get_rsid_vals)))
  if (nrow(locs)>0) {
    bed<-cbind(chr,locs,score,strand)
    rownames(bed)<-NULL
    write.table(bed,paste("/project2/nobrega/grace/expression/",trait,"/",tissue,"/results/posthoc/beds/",gene,"_TWAS.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

  } 
}

make_TWAS_bed(gene)
