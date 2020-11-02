import os, glob, subprocess, argparse
parser = argparse.ArgumentParser()
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("data_source", help="source of data (CMC, GTEx)")
parser.add_argument("out", help="Output directory")
parser.add_argument("sumstats", help="GWAS summary stats")
parser.add_argument("LD", help="chromosome to analyze, in format 'chr1'")
args = parser.parse_args()

######### Author: Grace Hansen #########
#This script takes weights in a chromosome-specific WEIGHTS_* directory and runs TWAS on them.
#This script is part of the snakemake TWAS pipeline.

##############Run TWAS#################

tissue=args.tissue
data_source=args.data_source

chromnum = args.chrom.split('r')[-1:][0]

#Move all WEIGHTS_* contents to WEIGHTS
if not os.path.isdir("%s/WEIGHTS"%(args.tissue)):
	os.mkdir("%s/WEIGHTS"%(args.tissue))
cmd = "mv %s/WEIGHTS_%s/* %s/WEIGHTS && rm -r %s/WEIGHTS_%s" % (args.tissue,chromnum,args.tissue,args.tissue,chromnum)
try:
	os.system(cmd)
except:
	print "Weights for chrom %s have already been moved"%chromnum

#make WEIGHTS*.pos file
posfile="%s/WEIGHTS/GTEx.v7.%s.pos" % (args.tissue,chromnum)
pos = open(posfile,'w')
pos.write("WGT\tID\tCHR\tP0\tP1\n")

#Get TSS from gene_exp file
weights=glob.glob("%s/WEIGHTS/*RDat"%(args.tissue))
for ln in weights:
	line=ln.strip().split('/')[2]
	gene = ln.strip().split('.wgt.RDat')[0].split('/')[2]
	if data_source=="GTEx_v8":
		cmd="""zcat /project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz | grep "gene_name "\\""%s""\\"" | awk '/\t'gene'\t/'"""%(gene)
	else:
		cmd="""zcat /project2/nobrega/grace//gencode.v19.annotation_hg19.gtf.gz | grep "gene_name "\\""%s""\\"" | awk '/\t'gene'\t/'"""%(gene)
	info = subprocess.check_output(cmd,shell=True)
	chrom,source,type,s,e = info.strip().split()[0:5]
	posline=("%s\t%s\t%s\t%s\t%s\n") % (line, gene, chromnum, s, e)
	pos.write(posline)
pos.close()

cmd = """Rscript /project2/nobrega/grace/expression/scripts/fusion_twas-master/FUSION.assoc_test.R \
--sumstats %s \
--weights ./%s/WEIGHTS/GTEx.v7.%s.pos \
--weights_dir ./%s/WEIGHTS/ \
--ref_ld_chr %s \
--chr %s \
--out %s%s.%s.dat""" % (args.sumstats, args.tissue, chromnum, args.tissue, args.LD, chromnum, args.out, args.data_source, chromnum)
os.system(cmd)