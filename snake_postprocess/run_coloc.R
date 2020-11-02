#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: run_coloc.R <trait> <tissue> <gene> <offset>\n", call.=FALSE)
}

library(coloc)
library(data.table)
library(tidyverse)
library(gridExtra)
source ("~/midway/expression/scripts/coloc/claudia.R")
trait=args[1]
tissue=args[2]
data_source=args[3]
gene=args[4]
window=args[5]
exp_N=args[6]

#Get chrom
if (data_source=="GTEx_v8") {
  cmd=paste("zcat /project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz | grep \"gene_name \\\"",gene,"\\\"\" | awk '/\t'gene'\t/'",sep='')
} else {
  cmd=paste("zcat /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf.gz | grep \"gene_name \\\"",gene,"\\\"\" | awk '/\t'gene'\t/'",sep='')
}
chrom<-strsplit(strsplit(system(cmd,intern=TRUE),'\t')[[1]][1],'r')[[1]][2]

#Grab p values for both datasets
get_coloc_data<-function(dataset1_name,dataset2_name) { ##Names must be filenames in the same directory as wd
  dataset1<-fread(dataset1_name)
  dataset2<-fread(dataset2_name)
  
  ##Make sure both datasets have a p-value or beta
  get_P_or_B<-function(dataset) {
    if ("P" %in% colnames(dataset)) {
      test_type="P"
    } else {
      if ("B" %in% colnames(dataset)) {
        test_type="B"
      } else if ("Z" %in% colnames(dataset)) {
        dataset$P=2*pnorm(-abs(dataset$Z))
        test_type="P"
      } else {
        print("p value or beta cannot be found for one of the datasets, cannot run coloc")
      }
    }
    return(list(data=dataset,test_type=test_type))
  }

  #Find N for both datasets
  find_N<-function(dataset,dataset_name) {
    if ("N" %in% colnames(dataset)) {
      N<-dataset$N
    } else {
      N<-as.numeric(exp_N) #If there isn't an 'N' colum in the data, there needs to be a [filename]_N file from which to read the N
    }
    return(N)
  }

  dset1<-get_P_or_B(dataset1)$data
  dset1_test_type<-get_P_or_B(dataset1)$test_type
  dset2<-get_P_or_B(dataset2)$data
  dset2_test_type<-get_P_or_B(dataset2)$test_type
  dset1_N<-find_N(dataset1,dataset1_name)
  dset2_N<-find_N(dataset2,dataset2_name)
  
  if (dset1_test_type==dset2_test_type) {
    return(list(dset1=dset1,dset2=dset2,dset1_N=dset1_N,dset2_N=dset2_N,test_type=dset1_test_type))
  } else {
    print("One dataset has p-values and one has beta values; this needs to be fixed before coloc can be run")
  }
}

data<-get_coloc_data(paste("coloc/loci/",gene,"_GWAS_sumstats",sep=''),paste("coloc/loci/",gene,"_exp.txt",sep=''))
GWAS<-data$dset1
eQTL<-data$dset2
eQTL_N<-data$dset2_N
test_type<-data$test_type
GWAS<- GWAS %>% distinct(rsid,.keep_all = TRUE)
eQTL<- eQTL %>% distinct(rsid,.keep_all = TRUE)

if (nrow(GWAS)>0 && nrow(eQTL>0)) {
  #Window data to desired range
  if (data_source=="GTEx_v8") {
    cmd=paste("zcat /project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz | grep \"gene_name \\\"",gene,"\\\"\" | awk '/\t'gene'\t/'",sep='')
  } else {
    cmd=paste("zcat /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf.gz | grep \"gene_name \\\"",gene,"\\\"\" | awk '/\t'gene'\t/'",sep='')
  }
  ln<-strsplit(system(cmd,intern=TRUE),'\t')
  if (ln[[1]][7]=="+") {
    TSS=ln[[1]][4]
  } else {
    TSS=ln[[1]][5]
  }
  TSS=as.numeric(TSS)
  genes<-NULL
  
  if (!("genom_distance" %in% names(eQTL))) {
    GWAS$tss_distance<-eQTL$tss_distance
    eQTL$genom_distance<-eQTL$tss_distance+TSS
    GWAS$genom_distance<-eQTL$tss_distance+TSS
  } else {
    eQTL$tss_distance<-eQTL$genom_distance-eQTL$TSS
    GWAS$genom_distance<-eQTL$genom_distance
    GWAS$tss_distance<-eQTL$tss_distance
  }
  GWAS<-GWAS[abs(GWAS$tss_distance)<=window,]
  eQTL<-eQTL[abs(eQTL$tss_distance)<=window,]
  
  MAF<-eQTL$MAF
  
  ##Determine which snp is the likely causal variant in each dataset (we're assuming 1 here, which is probs wrong)
  eQTL_finemap<-finemap.abf(dataset=list(pvalues=eQTL$P, N=eQTL_N, MAF=MAF, snp=eQTL$rsid,type="quant"))
  GWAS_finemap<-finemap.abf(dataset=list(pvalues=GWAS$P, N=GWAS$N, MAF=MAF, snp=GWAS$rsid,type="quant"))
  
  ##Perform coloc.abf test to test whether the causal variant is the same
  if (test_type=="P") {
    coloc <- coloc.abf(dataset1=list(pvalues=GWAS$P, N=GWAS$N, snp=GWAS$rsid, type="quant"),
                       dataset2=list(pvalues=eQTL$P, N=eQTL_N, snp=eQTL$rsid, type="quant"),
                       MAF=MAF)
  } else {
    coloc <- coloc.abf(dataset1=list(beta=GWAS$B, varbeta=GWAS$var, N=GWAS$N, snp=GWAS$rsid, type="quant"),
                       dataset2=list(beta=eQTL$B, varbeta=eQTL$var, N=eQTL_N, snp=eQTL$rsid, type="quant"),
                       MAF=MAF)
  }
  
  #Determine winning model. H4 means likeliest scenario is shared causal variant. H0=no causal variants, H1=only 1st dset has causal variant, 
  #H2=only 2nd dset has causal variant, H3=both have causal variants but they're different
  model<-names(which.max(coloc$summary[2:6]))
  shared_var=FALSE
  if (model=="PP.H4.abf") {
    shared_var=TRUE
  }
  
  ##Plots
  eQTL_plot<-head(eQTL_finemap,-1)
  GWAS_plot<-head(GWAS_finemap,-1)
  
  eQTL_plot$genom_distance<-eQTL$genom_distance
  GWAS_plot$genom_distance<-GWAS$genom_distance
  
  #Get index of most probable causal SNP
  max_eQTL<-eQTL_plot[which.max(eQTL_plot$SNP.PP),]
  max_GWAS<-GWAS_plot[which.max(GWAS_plot$SNP.PP),]
  
  #Plot of -log10 p values of each
  Q<-ggplot(eQTL_plot)+
    geom_point(aes(x=genom_distance,y=SNP.PP),color="turquoise4")+
    ggtitle(paste(gene," eQTL p-values",sep=""))+
    ylab("log10 PP of SNP being causal")+
    xlab(paste("Chromosome ",chrom,sep=''))+
    annotate("point",x=max_eQTL$genom_distance,y=max_eQTL$SNP.PP,color="darkred",size=3)+
    annotate("text",x=max_eQTL$genom_distance,y=max_eQTL$SNP.PP,vjust=-0,hjust=-0.1,label=max_eQTL$snp,size=3)+  
    theme_minimal()
  
  G<-ggplot(GWAS_plot)+
    geom_point(aes(x=genom_distance,y=SNP.PP),color="steelblue4")+
    ggtitle(paste(trait," GWAS summary statistics",sep=""))+
    ylab("log10 PP of SNP being causal")+
    xlab(paste("Chromosome ",chrom,sep=''))+
    annotate("point",x=max_GWAS$genom_distance,y=max_GWAS$SNP.PP,color="darkred",size=3)+
    annotate("text",x=max_GWAS$genom_distance,y=max_GWAS$SNP.PP,vjust=-0,hjust=-0.1,label=max_GWAS$snp,size=3)+  
    theme_minimal()
  
  if (shared_var==TRUE) {
    #Plot of posterior probability of H4
    coloc_plot<-coloc$results
    genom_dists<-eQTL[,c("rsid","genom_distance")]
    coloc_plot<-merge(coloc_plot,genom_dists,by.x="snp",by.y="rsid")
    max_coloc<-coloc_plot[which.max(coloc_plot$SNP.PP.H4),]
    
    C<-ggplot(coloc_plot)+
      geom_point(aes(x=genom_distance,y=SNP.PP.H4),color="darkgreen")+
      ggtitle(paste("Probability that SNP is causal for ",gene," eQTL and GWAS",sep=""))+
      ylab("log10 posterior probability")+
      xlab(paste("Chromosome ",chrom,sep=''))+
      annotate("point",x=max_coloc$genom_distance,y=max_coloc$SNP.PP.H4,color="darkred",size=3)+
      annotate("text",x=max_coloc$genom_distance,y=max_coloc$SNP.PP.H4,vjust=-0,hjust=-0.1,label=max_coloc$snp,size=3)+  
      annotate("text",x=quantile(coloc_plot$genom_distance,0.3),y=0.5,label=paste("PP of shared variant: \n",coloc$summary[model],sep=''))+
      theme_minimal()
    
    pdf(paste("coloc/",gene,"_",trait,"_coloc_plots.pdf",sep=""),width=8,height=10)
    grid.arrange(Q,G,C,nrow=3)
  } else {
    pdf(paste("coloc/",gene,"_",trait,"_coloc_plots.pdf",sep=""),width=8,height=7)
    grid.arrange(Q,G,nrow=2)
  }
  
  
  write.table(eQTL_finemap,paste("coloc/eQTL_",gene,"_coloc_finemap",sep=""),quote=FALSE,row.names=FALSE,sep='\t')
  write.table(GWAS_finemap,paste("coloc/",gene,"_GWAS_coloc_finemap",sep=""),quote=FALSE,row.names=FALSE,sep='\t')
  write.table(coloc$summary,paste("coloc/",gene,"_coloc_summary",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE,sep='\t')
  write.table(coloc$results,paste("coloc/",gene,"_coloc_full",sep=""),quote=FALSE,row.names=FALSE,sep='\t')
  write("done",paste("coloc/",gene,"_coloc_status",sep=""))
} else {
  write("couln't run; empty locus files",paste("coloc/",gene,"_coloc_status",sep=""))
}
