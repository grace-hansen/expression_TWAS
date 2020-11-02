import os
import argparse
parser = argparse.ArgumentParser()
#Mandatory argument: which chromosome
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("data_source", help="source of data (CMC, GTEx)")
parser.add_argument("out", help="Output directory")
parser.add_argument("sumstats", help="GWAS summary stats")
parser.add_argument("LD", help="chromosome to analyze, in format 'chr1'")
args = parser.parse_args()

chromnum = args.chrom.split('r')[-1:][0]

######### Author: Grace Hansen #########
#This script is run as part of the Snakemake pipeline to run TWAS.
#This script collects TWAS data and preforms Bonferroni correction and plots loci after conditional testing.

##############TWAS postprocessing#################
numweights = sum(1 for f in os.listdir("%s/WEIGHTS"%(args.tissue)) if f.endswith(".wgt.RDat"))
cmd = "cat %s%s.%s.dat | awk 'NR == 1 || $NF < 0.05/%s' > %s%s.%s.dat.top" % (args.out, args.data_source, chromnum, numweights, args.out, args.data_source, chromnum)
os.system(cmd)

cmd = """Rscript /project2/nobrega/grace/expression/scripts/fusion_twas-master/FUSION.post_process.R \
--sumstats %s \
--tissue %s \
--input %s%s.%s.dat.top \
--out %s%s.%s.dat.top.analysis \
--ref_ld_chr %s \
--chr %s \
--plot \
--plot_individual \
--locus_win 100000""" % (args.sumstats, args.tissue, args.out, args.data_source, chromnum, args.out, args.data_source, chromnum, args.LD, chromnum)
os.system(cmd)