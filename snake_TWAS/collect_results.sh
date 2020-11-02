#! /bin/sh
######### Author: Grace Hansen #########
#This script collects results from TWAS.
#This script is part of the snakemake TWAS pipeline.
trait=$1
tissue=$2
data_source=$3

cd /project2/nobrega/grace/expression/$trait/$tissue/results

#Collect significant genes
head -1 ${data_source}.1.dat.top > ${data_source}.all.dat.top
for i in $(seq 1 22); do
    sed '1d' ${data_source}.$i.dat.top >> temp
done
sort -g -k19,19 temp >> ${data_source}.all.dat.top
rm temp

#Collect all tested genes
head -1 ${data_source}.1.dat > ${data_source}.all.dat
for i in $(seq 1 22); do
    sed '1d' ${data_source}.$i.dat >> temp
done
sort -g -k19,19 temp >> ${data_source}.all.dat
rm temp

#Make list of significant genes
mkdir posthoc
sed '1d' ${data_source}.all.dat.top | cut -f2  > posthoc/sig_genes