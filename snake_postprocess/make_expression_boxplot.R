#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
arguments <- parse_args(OptionParser(usage = "make_expression_boxplot.R <trait> <tissue> <gene>",option_list=list()),
                        positional_arguments = 4)
opt<-arguments$opt
trait<-arguments$args[1]
tissue<-arguments$args[2]
gene<-arguments$args[3]
data_source<-arguments$args[4]
setwd(paste("/project2/nobrega/grace/expression/",trait,"/",tissue,"/results/posthoc",sep=""))
if (data_source=="GTEx_v7") {
  chrom=system(paste("zcat ../../",tissue,"_",data_source,"_exp.txt.gz | awk '{ if ($3 == \"",gene,"\") print $1}'",sep=''),intern=TRUE)
} else if (data_source=="CMC") {
  cmd=paste("zcat /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf.gz | grep \"gene_name \\\"",gene,"\\\"\" | awk '/\t'gene'\t/'",sep='')
  info=system(cmd,intern=TRUE)
  chrom=strsplit(strsplit(info,'\t')[[1]][1],'r')[[1]][2]
} else {
  cmd=paste("zcat /project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz | grep \"gene_name \\\"",gene,"\\\"\" | awk '/\t'gene'\t/'",sep='')
  info=system(cmd,intern=TRUE)
  chrom=strsplit(strsplit(info,'\t')[[1]][1],'r')[[1]][2]
}
pony_colors<-fread("/project2/nobrega/grace/pony_palette")

#Count number of subjects per allele
count_genos<-function(gene_dat,rsid) {
  plink_rsids<-fread(paste(gene,".map",sep=''))
  line<-which(plink_rsids$V2==rsid)
  plink_genos<-fread(paste(gene,".ped",sep=''))
  IDs<-plink_genos[,1]
  plink_genos<-plink_genos[,7:ncol(plink_genos)]
  genos<-as.data.frame(plink_genos)[,line]
  genos<-gsub(" ", "/", genos)
  genos<-cbind(IDs,genos)
  allele_counts<-matrix(nrow=2,ncol=0)
  for (c in c("A","C","T","G")) {
    num<-sum(str_count(genos$genos,c))
    allele_counts<-cbind(allele_counts,as.data.frame(c(c,num),stringsAsFactors = FALSE))
  }
  allele_counts<-as.data.frame(t(allele_counts),stringsAsFactors = FALSE) %>% arrange(.,V2)
  write.table(t(allele_counts[3:4,]),file=paste("reports/",gene,"_report",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  geno_exp<-merge(gene_dat,genos,by.x="P1",by.y="V1")
  geno_exp$Expression<-as.numeric(as.character(geno_exp$Expression))
  return(geno_exp)
}

#Get rsid: coloc'ed rsid if in coloc, top finemapped eQTL rsid if not in coloc
coloc_genes<-fread("coloc/colocalizing_genes_rsids",header=FALSE)
if (gene %in% coloc_genes$V1) {
  rsid_dat<-fread(paste("coloc/",gene,"_coloc_full",sep=''))
  rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP.H4)]
} else if (length(which(list.files("coloc/")==paste("eQTL_",gene,"_coloc_finemap",sep='')))>0) {
  rsid_dat<-fread(paste("coloc/eQTL_",gene,"_coloc_finemap",sep=''))
  rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP)] 
} else {
  cmd=paste("awk -F '\t' '$2==\"",gene,"\" {print $9}' ../",data_source,".all.dat.top",sep='')
  rsid=system(cmd,intern=TRUE)
}
if (rsid=='null') {
  cmd=paste("awk -F '\t' '$2==\"",gene,"\" {print $9}' ../",data_source,".all.dat.top",sep='')
  rsid=system(cmd,intern=TRUE)
}

#Make expression phenotpe file
if (data_source=="GTEx_v8") {
  tpms<-fread(cmd="zcat ~/midway/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | tail -n +3")
  samples<-scan(paste("~/midway/expression/",trait,"/",tissue,"/",tissue,"_sample_IDs",sep=''),what="character",sep="\n")
  colnames<-c("Description",samples)
  tpms<-tpms[,..colnames]
  tpms<-tpms[tpms$Description==gene]
  tpms$Description<-NULL
  subj_IDs<-paste("GTEX-",sapply(strsplit(names(tpms),'-'),'[[',2),sep='')
  names(tpms)<-subj_IDs
  write.table(tpms,paste(gene,"_pheno",sep=''),row.names = FALSE,quote=FALSE)
  gc()
} else if (data_source=="CMC") {
  cmd=paste("zcat ../../",tissue,"_",data_source,"_tpm.txt.gz | head -1 > ",gene,"_pheno",sep='')
  system(cmd)
  cmd=paste("zcat ../../",tissue,"_",data_source,"_tpm.txt.gz | awk '{ if ($3 == \"",gene,"\") print $0}' >> ",gene,"_pheno",sep='')
  system(cmd)
} else {
  cmd=paste("zcat ../../",tissue,"_",data_source,"_exp_1.txt.gz | head -1 > ",gene,"_pheno",sep='')
  system(cmd)
  cmd=paste("zcat ../../",tissue,"_",data_source,"_exp_",chrom,".txt.gz | awk '{ if ($3 == \"",gene,"\") print $0}' >> ",gene,"_pheno",sep='')
  system(cmd)
}
gene_dat<-fread(paste(gene,"_pheno",sep=''))
gene_dat<-as.data.frame(t(gene_dat),stringsAsFactors = FALSE)
rownames<-rownames(gene_dat)[4:length(rownames(gene_dat))]
gene_dat<-as.data.frame(gene_dat[4:nrow(gene_dat),])
colnames(gene_dat)<-"Expression"
if (data_source=="GTEx_v7" || data_source=="GTEx_v8") {
  IDs<-paste("GTEX-",sapply(strsplit(rownames,'-'),'[[',2),sep='')
} else {
  IDs<-rownames
}
gene_dat$P1<-IDs
gene_dat$P2<-IDs
gene_dat<-gene_dat %>% select(P1,P2,everything())
write.table(gene_dat,paste(gene,".pheno",sep=''),row.names = FALSE,col.names = FALSE,quote=FALSE)

#Add rsid genotypes to expression file from plink
cmd=paste("awk '{ if ($2 == \"",rsid,"\") print $4}' /project2/nobrega/grace/genos/",data_source,"/plink/*bim",sep="")
rsid_loc=as.numeric(system(cmd,intern=TRUE))
start=rsid_loc-5
stop=rsid_loc+5
cmd<-paste("plink --bfile /project2/nobrega/grace/genos/",data_source,"/plink/chr",chrom,"_rsids --pheno ",gene,".pheno --make-bed --out ",gene," --keep ",gene,".pheno --chr ",chrom,
           " --from-bp ",start," --to-bp ",stop,sep="")
system(cmd)
cmd<-paste("plink --bfile ",gene," --recode tab --out ",gene,sep='')
system(cmd)

geno_exp<-count_genos(gene_dat,rsid)
geno_exp<-geno_exp[geno_exp$genos!="0/0",]
#Make boxplot of gene expression per genotype
pdf(paste("boxplots/",gene,"_expression_boxplot.pdf",sep=''),width=4,height=4)
ggplot(geno_exp,aes(x=genos,y=Expression,fill=genos))+
  geom_boxplot(width=0.5,outlier.shape=NA)+
  geom_jitter(shape=16,size=1.5,position=position_jitter(0.1))+
  ylab("Expression (tpm)")+
  xlab(paste(rsid," genotype",sep=''))+
  ggtitle(gene)+
  scale_fill_manual(values=c(rgb(pony_colors[8,1:3]),rgb(pony_colors[10,1:3]),rgb(pony_colors[13,1:3])))+
  guides(fill=FALSE)+
  theme_minimal()+
  theme(axis.title.x=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=15))
dev.off()

system(paste("rm ",gene,"_pheno ",gene,".*",sep=''))
