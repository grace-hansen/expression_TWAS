#!/usr/bin/python
import pandas as pd
import subprocess
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("trait", help="trait (name of directory in which tissue is located")
parser.add_argument("tissue", help="tissue")
parser.add_argument("genelist", help="genelist, in which each entry is in the format gene\trsid")
parser.add_argument("prefix", help="prefix for output file")
parser.add_argument("data_source", help="transcriptome data source, e.g. GTEx, CMC")
args = parser.parse_args()

##############Input: text file containing list of gene-rsid pairs, which represent the colocalizing variants for each TWAS result gene
    #For each gene-variant pair, get snpeff annotation and chromHMM annotation, and put into output file
trait=args.trait
tissue=args.tissue
genelist=args.genelist
data_source=args.data_source

########### get location annotation by parsing gtf #####################

def snpeff_annotation(rsid,chrom):
    snpeff=open("/project2/nobrega/grace/genos/snpeff/chr%s_rsid_anno.vcf"%(chrom),'r')
    #Get line in snpeff vcf corresponding to rsid
    rsid_verdict=''
    for ln in snpeff:
        if ln[0]!='#':
            chrom,pos,snpeff_rsid,ref,alt,score,PASS,anno=ln.strip().split('\t')[0:8]
            if snpeff_rsid==rsid:
                rsid_verdict=anno
    return(rsid_verdict)

########### get chromHMM annotation #####################
def chromHMM_annotation(rsid):
    rsid_verdict=''
    if data_source== "GTEx_v8":
        chromHMM=open("/project2/nobrega/grace/expression/%s/%s/epigenomeroadmap_%s_chromHMM_hg38.bed"%(trait,tissue,tissue),'r')
    else:
        chromHMM=open("/project2/nobrega/grace/expression/%s/%s/epigenomeroadmap_%s_chromHMM_hg19.bed"%(trait,tissue,tissue),'r')
    for ln in chromHMM:
        chr,start,stop,anno=ln.strip().split('\t')[0:4]
        if chr=='chr'+chrom:
            if rsid >= start and rsid < stop:
                rsid_verdict=anno
    return(rsid_verdict)


#####################

out=open("/project2/nobrega/grace/expression/%s/%s/results/posthoc/%s"%(trait,tissue,args.prefix),'w')
genelist=open(genelist,'r')
for ln in genelist:
    if ln != '':
        gene,rsid=ln.strip().split('\t')
        cmd="""awk -F '\t' ' {if ($2 == "%s") print $0}' ~/midway/genos/%s/plink/*bim"""%(rsid,data_source)
        std_out=subprocess.check_output(cmd,shell=True).rstrip()
        rsid_pos=std_out.split('\t')[3]
        chrom=std_out.split('\t')[0]
        if chrom[0]=='c':
            chrom=chrom.split('r')[1]
        snpeff_anno=snpeff_annotation(rsid,chrom)
        chromHMM_anno=chromHMM_annotation(rsid_pos)
        out.write('\t'.join([chrom,gene,rsid,rsid_pos,chromHMM_anno,snpeff_anno])+'\n')
out.close()

if not os.path.isdir("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/beds"%(trait,tissue)):
    os.mkdir("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/beds"%(trait,tissue))
cmd="""cat /project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/coloc_annotation.txt | while read line; do 
    gene=$(echo $line | cut -f2 -d' '); chr=$(echo $line | cut -f1 -d' '); 
    rsid=$(echo $line | cut -f3 -d' '); loc=$(echo $line | cut -f4 -d' '); 
    loc1=$(($loc +1)); echo -e "chr$chr\t$loc\t$loc1\t$rsid\t.\t." > /project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/beds/${gene}_${rsid}_coloc.bed; 
    bgzip /project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/beds/${gene}_${rsid}_coloc.bed; 
    tabix -p bed /project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/beds/${gene}_${rsid}_coloc.bed.gz; done"""%(trait,tissue,trait,tissue,trait,tissue,trait,tissue)
subprocess.call(cmd,shell=True)