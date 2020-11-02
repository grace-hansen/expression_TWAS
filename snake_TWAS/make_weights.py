import os, glob, argparse, commands, gzip, shutil
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("trait", help="tissue being analyzed")
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("data_source", help="origin of expression data (e.g. GTEx, CMC")
parser.add_argument("offset", help="size of offset in scientific notation, e.g. 5e4")
parser.add_argument("covars", help="Covariates to include when making TWAS weights")
parser.add_argument("temp", help="Location of temp directory for snakefile flag")
args = parser.parse_args()

######### Author: Grace Hansen #########
# This  script takes genotype and phenotype(here, gene expression in a [tissue]_[source]_exp.txt format) data and prepares plinks for TWAS.
# It then takes these plinks and  computes TWAS weights.
#This script is part of the snakemake TWAS pipeline. 

###################################################################
chromnum = args.chrom.split('r')[-1:][0]
offset = float(args.offset)
trait=args.trait
tissue=args.tissue
data_source=args.data_source
covars=args.covars

#Make list of all variants in GWAS summary statistics to filter by
if not os.path.exists("GWAS_rsids_%s"%(chromnum)):
	cmd="cut -f2 /project2/nobrega/grace/genos/%s/plink/chr%s_rsids.bim"%(args.data_source,chromnum)
	status,rsids=commands.getstatusoutput(cmd)
	rsids=list(rsids.split('\n'))
	GWAS=pd.read_csv("%s_GWAS_sumstats"%(trait),sep='\t')
	with open("GWAS_rsids_%s"%(chromnum),'w') as file:
		GWAS_rsids=GWAS["SNP"].tolist()
		intersect=list(set(rsids).intersection(GWAS_rsids))
		intersect='\n'.join(intersect)
		file.write(intersect+'\n')
	del GWAS
GWAS_rsid_file="GWAS_rsids_%s"%(chromnum)

subject_IDs="%s/%s_subject_IDs"%(tissue,tissue)

gene_exp_mat=gzip.open('%s/%s_%s_exp_%s.txt.gz'%(tissue,tissue,data_source,chromnum),'r')
for ln in gene_exp_mat:
	if ln[0] == "c":
		pass
	else:
		c, coord, gene = ln.split()[:3]
		print type(c)
		print type(chromnum)
		print c,chromnum
		if c == str(chromnum):
			vals = ln.split()[3:]
			#Make phenotype file of subject IDs + gene expression to filter plink by
			with open("%s_gene_exp"%(gene),'w') as file:
				file.write('\n'.join(vals)+'\n')
			gene_vals="%s_gene_exp"%(gene)
			cmd="paste %s %s > %s_IDs_gene_exp"%(subject_IDs,gene_vals,gene)
			os.system(cmd)
			IDs_gene_exp_file="%s_IDs_gene_exp"%(gene)
			#Define window within which to grab variants
			p0=int(coord)-offset
			if p0 < 0:
				p0 = 0
			p1=int(coord)+offset
			#Create plinks to run compute_weights on
			cmd="plink --bfile /project2/nobrega/grace/genos/%s/plink/chr%s_rsids \
			--pheno %s --make-bed --out %s/%s --keep %s \
			--chr %s --from-bp %s --to-bp %s --extract %s"%(data_source,chromnum,IDs_gene_exp_file,tissue,gene,IDs_gene_exp_file,chromnum,p0,p1,GWAS_rsid_file)
			print cmd
			os.system(cmd)
			os.remove("%s_gene_exp"%(gene))
			os.remove("%s_IDs_gene_exp"%(gene))
			#Run compute_weights, delete plink and temp files after
			tmp = "%s/tmp_%s" % (tissue,chromnum)
			if not os.path.isdir(tmp):
				os.mkdir(tmp)
			weightsdir = "%s/WEIGHTS_%s" % (tissue,chromnum)
			if not os.path.isdir(weightsdir):
				os.mkdir(weightsdir)
			if covars == "None":
				cmd = """Rscript /project2/nobrega/grace/expression/scripts/fusion_twas-master/FUSION.compute_weights.R \
				--bfile %s/%s \
				--tmp %s/tmp_%s/ \
				--save_hsq \
				--out %s/WEIGHTS_%s/%s \
				--PATH_gcta /project2/nobrega/grace/expression/scripts/gcta_1.91.4beta/gcta64 \
				--models top1,blup,lasso,enet""" % (tissue, gene, tissue, chromnum, tissue, chromnum, gene)
			else:
				cmd = """Rscript /project2/nobrega/grace/expression/scripts/fusion_twas-master/FUSION.compute_weights.R \
				--bfile %s/%s \
				--covar %s \
				--tmp %s/tmp_%s/ \
				--save_hsq \
				--out %s/WEIGHTS_%s/%s \
				--PATH_gcta /project2/nobrega/grace/expression/scripts/gcta_1.91.4beta/gcta64 \
				--models top1,blup,lasso,enet""" % (tissue, gene, covars, tissue, chromnum, tissue, chromnum, gene)
			print(cmd)
			os.system(cmd)
			for p in glob.glob("%s/%s.*"%(tissue,gene)):
				os.remove(p)
			os.system("cat %s/WEIGHTS_%s/%s.hsq >> %s/WEIGHTS_%s/chr_%s.hsq && rm %s/WEIGHTS_%s/%s.hsq"%(args.tissue, chromnum, gene, args.tissue, chromnum, chromnum, args.tissue, chromnum, gene))

#Remove temp directories if empty
shutil.rmtree(tmp,ignore_errors=True)

#Create temp file indicating TWAS weights have been made
temp="%smake_weights_%s.txt" % (args.temp,chromnum)
open(temp,'a').close()