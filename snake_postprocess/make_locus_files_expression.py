import gzip, os, argparse, subprocess
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("trait", help="trait (name of directory in which tissue is located")
parser.add_argument("tissue", help="tissue (must be the prefix for expression file, e.g. 'adipose_subcutaneous.allpairs.txt.gz'")
parser.add_argument("data_source", help="data source, e.g. GTEx_v8 or CMC")
parser.add_argument("gene", help="gene for which to extract information.")
parser.add_argument("GWAS_sumstats", help="full path to sumstats that will be used for coloc. Will only include SNPs in GWAS sumstats. rsid must be first column")
parser.add_argument("lookup", help="lookup file matching GTEx variant IDs to rsids")
parser.add_argument("LD", help="LD reference panel being used (for ref/alt allele matching")

##############Grab gene rsids#################
args = parser.parse_args()
data_source=args.data_source
gene=args.gene
trait=args.trait
tissue=args.tissue
GWAS_sumstats=args.GWAS_sumstats
LD=args.LD

def make_gene_file(gene,trait,tissue):
    #Given the GTEx file for all tested eQTL-gene relationships in a specific tissue, extract only the ones for your gene of interest,
        # replace the ENSG with gene name, replace the GTEx variant_id with rsid, and put into a new file.
    print "writing eQTL entries for %s..."%(gene)
    out=gzip.open("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/%s_eQTL_%s_rsid.txt.gz"%(trait,tissue,tissue,gene),'wt')
    eQTL_dat=gzip.open("/project2/nobrega/grace/expression/%s/%s/%s.allpairs.txt.gz"%(trait,tissue,tissue),'rt')
    header=eQTL_dat.readline()

    if data_source=="GTEx_v7":
        out.write(header)
        #Make dictionary of rsids: gene IDs
        lookup=pd.read_csv(args.lookup,sep='\t',compression='gzip')
        rsid_lookup=dict(zip(lookup.variant_id,lookup.rs_id_dbSNP147_GRCh37p13))
        del lookup

        #Make dictionary of ENSG : gene name
        coding=pd.read_csv("/project2/nobrega/grace/expression/coding_genes.bed",sep='\t',names=['chrom','start','stop','gene','score','strand','ENSG'])
        gene_lookup=dict(zip(coding.gene,coding.ENSG))
        gene_ENSG=gene_lookup.get(gene,'.')
        print("%s: %s"%(gene,gene_ENSG))
        del coding

        for ln in eQTL_dat:
            ENSG, var, tss, ma_samples, ma_count, maf, pval, slope, slope_se = ln.split()[:9]
            if ENSG==gene_ENSG:
                rsid=rsid_lookup.get(var,'.')
                if rsid != '.':
                    begin = "\t".join([args.gene, rsid, TSS, ma_samples, ma_count, maf, pval, slope, slope_se]) #Defining TSS the way I do, with gencode ( as vs. the way GTEx is defining it)
                    newln = begin + "\n"
                    out.write(newln)

    if data_source=="GTEx_v8":
        out.write(header)
        #Make dictionary of rsids: gene IDs
        lookup=pd.read_csv(args.lookup,sep='\t',compression='gzip')
        rsid_lookup=dict(zip(lookup.variant_id,lookup.rs_id_dbSNP151_GRCh38p7))
        del lookup

        #Make dictionary of ENSG : gene name
        gene_lookup={}
        with gzip.open("/project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz",'r') as gtf:
            for ln in gtf:
                if ln[0]!='#':
                    chrom,havana,element,start,stop,dot,strand,dot2,info=ln.strip().split('\t')
                    type=info.split('"')[3]
                    gtf_gene=info.split('"')[5]
                    if gtf_gene==gene:
                        gene_ENSG=info.split('"')[1]
                        break
        print("%s: %s"%(gene,gene_ENSG))

        for ln in eQTL_dat:
            ENSG, var, tss, ma_samples, ma_count, maf, pval, slope, slope_se = ln.split()[:9]
            if ENSG==gene_ENSG:
                rsid=rsid_lookup.get(var,'.')
                if rsid != '.':
                    begin = "\t".join([gene, rsid, tss, ma_samples, ma_count, maf, pval, slope, slope_se]) #Defining TSS the way I do, with gencode ( as vs. the way GTEx is defining it)
                    newln = begin + "\n"
                    out.write(newln)


    if data_source=="CMC":
        out.write("gene\trsid\teQTL_ref\tTSS\tgenom_distance\tMAF\tP\tslope\n")
        ### Get TSS
        cmd="""zcat /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf.gz | grep "gene_name "\\""%s""\\"" | awk '/\t'gene'\t/'"""%(gene)
        info = subprocess.check_output(cmd,shell=True)
        strand=info.strip().split('\t')[6]
        if (strand == "+"):
            TSS=info.strip().split('\t')[3]
        else:
            TSS=info.strip().split('\t')[4]
        chrom=info.strip().split()[0].split('r')[1]

        #Make dictionary of MAFs
        cmd="plink --recode --freq --out plink_%s --bfile /project2/nobrega/grace/genos/CMC/plink/chr%s_rsids"%(gene,chrom)
        subprocess.check_call(cmd,shell=True)
        freqs=pd.read_csv("plink_%s.frq"%gene,delim_whitespace=True)
        MAF_dict=dict(zip(freqs.SNP,freqs.MAF))
        cmd="rm plink_%s.*"%gene
        subprocess.check_call(cmd,shell=True)
        del freqs

        #Make dictionary of locations
        locs=pd.read_csv("/project2/nobrega/grace/genos/CMC/plink/chr%s_rsids.bim"%chrom,delim_whitespace=True,header=None)
        locs_dict=dict(zip(locs.iloc[:,1],locs.iloc[:,3]))
        del locs

        for ln in eQTL_dat:
            eQTL_gene, rsid, tss_distance, pval_nominal, slope = ln.strip().split(' ')
            if eQTL_gene==gene:
                cmd="""awk -F '\t' '{if ($2 == "%s") print $6 }' /project2/nobrega/grace/genos/CMC/plink/chr%s_rsids.bim """%(rsid,chrom)
                ref = subprocess.check_output(cmd,shell=True).rstrip()
                MAF=MAF_dict.get(rsid,'.')
                if MAF != '.':
                    genom_dist=locs_dict.get(rsid,'.')
                    if genom_dist != '.':
                        newln = "\t".join([gene, rsid, ref, str(TSS), str(genom_dist), str(MAF), str(pval_nominal), str(slope)]) + "\n"
                        out.write(newln)

    out.close()

def merge_GWAS_eQTL(GWAS_sumstats,gene,trait,tissue):
    #Inner join the file made by make_gene_file and a supplied GWAS summary statistics file, so that the same rsids are in both. 
        #Save out these files with a *coloc suffix

    #Load GWAS sumstats and eQTL sumstats
    GWAS=pd.read_csv(GWAS_sumstats,sep='\t')
    GWAS.rename(columns={ GWAS.columns[0]: "rsid"},inplace=True)
    GWAS_cols=list(GWAS.columns)
    eQTL=pd.read_csv("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/%s_eQTL_%s_rsid.txt.gz"%(trait,tissue,tissue,gene),sep='\t',compression='gzip')
    if data_source=="GTEx_v8"   or data_source=="GTEx_v7":
        eQTL.rename(columns={ "variant_id": "rsid", "maf": "MAF", "pval_nominal": "P", 'slope' : 'B', 'slope_se': 'SE'},inplace=True)
    
    eQTL_cols=list(eQTL.columns)

    merged=pd.merge(GWAS,eQTL,on='rsid')

    ##Flip alleles
    print "Flipping alleles so eQTL, GWAS, and LD panel match..."
    flips=[]
    if data_source == "GTEx_v8" or data_source == "GTEx_v7":
        lookup=pd.read_csv(args.lookup,sep='\t',compression="gzip")
        if data_source=="GTEx_v8":
            eQTL_refs=pd.Series(lookup.ref.values,index=lookup.rs_id_dbSNP151_GRCh38p7).to_dict()
        else:
            eQTL_refs=pd.Series(lookup.ref.values,index=lookup.rs_id_dbSNP147_GRCh37p13).to_dict()
        del lookup
        for row in merged.itertuples():
            new_row=list(row)
            rsid=row.rsid
            cmd="""awk -F ' ' '{ if ($2 == "%s") print $6}' %s*.bim"""%(rsid,args.LD)
            LD_ref=subprocess.check_output(cmd,shell=True).decode("utf-8").rstrip()
            if LD_ref != "":
                eQTL_ref=eQTL_refs.get(rsid,'.')
                cmd="""awk -F '\t' '{ if ($1 == "%s") print $3}' %s """%(rsid,GWAS_sumstats)
                GWAS_ref=subprocess.check_output(cmd,shell=True).decode("utf-8").rstrip()
                if eQTL_ref != LD_ref:
                    B=float(new_row[12])
                    new_row[12]=-B
                if GWAS_ref != LD_ref:
                    A1, A2, Z = new_row[2],new_row[3],float(new_row[4])
                    new_row[2]=A2
                    new_row[3]=A1
                    new_row[4]=-Z
                new_row=new_row[1:]
                flips.append(new_row)
    if data_source == "CMC":
        for row in merged.itertuples():
            new_row=list(row)
            rsid=row.rsid
            eQTL_ref=row.eQTL_ref
            cmd="""awk -F '\t' '{ if ($1 == "%s") print $3}' %s """%(rsid,GWAS_sumstats)
            GWAS_ref=subprocess.check_output(cmd,shell=True).rstrip()
            if GWAS_ref != eQTL_ref:
                A1, A2, Z = new_row[2],new_row[3],float(new_row[4])
                new_row[2]=A2
                new_row[3]=A1
                new_row[4]=-Z
            new_row=new_row[1:]
            flips.append(new_row)

    merged_flipped=pd.DataFrame(flips,columns=merged.columns)
    print "Alleles flipped!"

    #Keep only GWAS sumstats where rsid is in eQTL dataset
    GWAS_out=merged_flipped.ix[:,GWAS_cols]
    GWAS_out.to_csv("coloc/loci/%s_GWAS_sumstats"%(gene),sep='\t',index=False)

    #Keep only eQTL sumstats where rsid is in GWAS dataset
    eQTL_out=merged_flipped.ix[:,eQTL_cols]
    eQTL_out.to_csv("coloc/loci/%s_exp.txt"%(gene),sep='\t',index=False)

    os.remove("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/%s_eQTL_%s_rsid.txt.gz"%(trait,tissue,tissue,args.gene))

try:
    os.makedirs("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/"%(trait,tissue))
except:
    pass

if not os.path.isfile("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/%s_exp.txt"%(
    args.trait,tissue,gene)) and not os.path.isfile("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/%s_GWAS_sumstats"%(
    trait,tissue,gene)):
    if not os.path.isfile("/project2/nobrega/grace/expression/%s/%s/results/posthoc/coloc/loci/%s_eQTL_%s_rsid.txt.gz"%(trait,tissue,tissue,gene)):
        make_gene_file(gene,trait,tissue)
        merge_GWAS_eQTL(GWAS_sumstats,gene,trait,tissue)
    else: 
        merge_GWAS_eQTL(GWAS_sumstats,gene,trait,tissue)
