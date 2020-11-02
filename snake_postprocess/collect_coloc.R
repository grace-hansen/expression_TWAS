#/usr/bin/Rscript
library(tidyverse)
library(data.table)
library(optparse)
arguments <- parse_args(OptionParser(usage = "collect_coloc.R <trait> <tissue> <data_source>",option_list=list()),
                        positional_arguments = 3)
opt=arguments$opt
trait=arguments$args[1]
tissue=arguments$args[2]
data_source=arguments$args[3]

setwd(paste("/project2/nobrega/grace/expression/",trait,"/",tissue,"/results/posthoc/coloc/",sep=''))
files<-list.files(pattern="*summary")

genes<-character()
rsids<-character()

for ( i in 1:length(files)) {
  dat<-fread(files[i])
  hyp<-which.max(dat$V2[2:6])
  if ( hyp == 5) {
    gene<-strsplit(files[i],"_")[[1]][1]
    rsid_dat<-fread(paste(gene,"_coloc_full",sep=''))
    rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP.H4)]
    
    genes<-c(genes,gene)
    rsids<-c(rsids,rsid)
  }
}

#Order results by TWAS p-value
out<-as.data.frame(cbind(genes,rsids),stringsAsFactors = FALSE)
res<-fread(paste("/project2/nobrega/grace/expression/",trait,"/",tissue,"/results/",data_source,".all.dat.top",sep=''))
out<-inner_join(res,out,by=c("ID"="genes"))
out<-out %>% select(ID,rsids)

write.table(out,"colocalizing_genes_rsids",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
