import os, commands, gzip, argparse, re
import pandas as pd
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("chrom", help="chromosome to analyze, in format '1'")
parser.add_argument("temp", help="Location of temp directory for snakefile flag")
parser.add_argument("data_source", help="Location of temp directory for snakefile flag")
parser.add_argument("--exclude", help="List of subject IDs to exclude (e.g. non-Europeans",default="None")
parser.add_argument("--search_terms", help="search terms to identify tissue in E-MTAB-5214. Max of 2",
    action="append",nargs="+")
args = parser.parse_args()

############## Author: Grace Hansen ###############
# For a given tissue, this script identifies samples for which there are both GTEx genotype data and gene expression data.
# It then writes out a subject ID list of these subjects, a sample ID list of the corresponding samples, and a matrix of all the expression values for these samples.
# These files will be used to make plinks and compute TWAS weights.#This script is part of a snakemake TWAS pipeline.

#######################################################################

###### Function to grab TSS of primary isoform from gencode gtf ##########

chromnum = args.chrom.split('r')[-1:][0]
data_source=args.data_source
tissue=args.tissue
temp=args.temp
exclude=args.exclude
search_terms=args.search_terms

#Grab every subject ID with a genotype from the *vcf with bcftools
cmd="cut -f1  -d' ' /project2/nobrega/grace/genos/%s/plink/chr%s_rsids.fam"%(data_source,chromnum)
status, output = commands.getstatusoutput(cmd)
geno_IDs=list(output.split('\n'))

#If exclude flag, remove ones that aren't white :/
if exclude != "None":
    with open(exclude,'r') as excl_file:
        exclude_IDs=excl_file.readlines()
        exclude_IDs = ' '.join(exclude_IDs).replace('\n','').split()
    geno_IDs=list(set(geno_IDs)-set(exclude_IDs))


#Open up output files
if data_source == "GTEx_v7":
    write_IDs=False
    if not os.path.exists("%s/%s_subject_IDs"%(args.tissue,args.tissue)) or not os.path.exists("%s/%s_sample_IDs"%(args.tissue,args.tissue)):
        write_IDs=True
        subject_IDs=open("%s/%s_subject_IDs"%(args.tissue,args.tissue),"w")
        sample_IDs=open("%s/%s_sample_IDs"%(args.tissue,args.tissue),"w")

    tissue_exp=open("%stemp_%s"%(args.temp,chromnum),"w")
    with open("/project2/nobrega/grace/expression/GTEx_exp.gct") as GTEx_exp:
        header=GTEx_exp.readlines()[1]
    tissue_exp.write(header)

    #For every subject with a genotype, find sample IDs from the tissue of interest
    for subject_ID in geno_IDs:
        terms=args.search_terms[0]
        if len(terms)==1 and terms != "None":
            term1=terms[0]
            cmd="grep %s /project2/nobrega/GTEx/E-MTAB-5214.sdrf.txt | grep %s | cut -f1 | sort -u"%(subject_ID,term1)
        else:
            term1,term2=terms
            cmd="grep %s /project2/nobrega/GTEx/E-MTAB-5214.sdrf.txt | grep %s | grep %s | cut -f1 | sort -u"%(subject_ID,term1,term2)
        status, output = commands.getstatusoutput(cmd)
        subject_sample_IDs=list(output.split('\n'))

        #For each sample ID from the tissue of interest for a given subject, look for gene expression data. If found, write sample and subject IDs to file, write gene expression data to file. 
        if len(subject_sample_IDs)==1:
            if not len(subject_sample_IDs[0])==0:
                cmd="grep %s /project2/nobrega/grace/expression/GTEx_exp.gct"%(subject_sample_IDs[0])
                status, output = commands.getstatusoutput(cmd)
                if len(output) != 0:
                    if write_IDs==True:
                        subject_IDs.write(subject_ID+'\t'+subject_ID+'\n')
                        sample_IDs.write(subject_sample_IDs[0]+'\n')
                    tissue_exp.write(output+'\n')
        else: 
            for subject_sample_ID in subject_sample_IDs:
                cmd="grep %s /project2/nobrega/grace/expression/GTEx_exp.gct"%(subject_sample_ID)
                status, output = commands.getstatusoutput(cmd)
                if len(output) != 0:
                    if write_IDs==True:
                        subject_IDs.write(subject_ID+'\t'+subject_ID+'\n')
                        sample_IDs.write(subject_sample_ID+'\n')
                    tissue_exp.write(output+'\n')
                    break
    if write_IDs==True:
        subject_IDs.close()
        sample_IDs.close()
    tissue_exp.close()

    #Transpose the output so genes are rows and GTEx IDs are columns
    transposed=pd.read_csv("%stemp_%s"%(args.temp,chromnum),sep='\t',index_col=0).T
    out=gzip.open("%s/%s_GTEx_v7_exp_%s.txt.gz"%(args.tissue,args.tissue,chromnum),'w')
    col_vals='\t'.join(str(v) for v in transposed.columns.values)
    out.write('chr'+'\t'+'coord'+'\t'+'gene'+'\t'+col_vals+'\n')

    #For each gene, add chr and TSS to start of its row
    for ln in transposed.itertuples():
        gene=ln[0]
        line='\t'.join(str(i) for i in ln)
        cmd="awk '/\t'%s'\t/' /project2/nobrega/grace/expression/coding_genes.bed"%(gene)
        status, output = commands.getstatusoutput(cmd)
        if len(output) != 0:
            chrom, s, e, gname, dot, strand, ENSG = output.split('\t')
            c=chrom.split('r')[1]
            if strand == "+":
                coord=s
            else:
                coord=e
            row=c+'\t'+str(coord)+'\t'+line+'\n'
            out.write(row)

    out.close()
    os.remove("%stemp_%s"%(args.temp,chromnum))

if data_source == "GTEx_v8":
    term=search_terms[0]
    term=term[0]
    write_IDs=False
    if not os.path.exists("%s/%s_subject_IDs"%(tissue,tissue)):
        write_IDs=True
        subject_IDs=open("%s/%s_subject_IDs"%(tissue,tissue),"w")

    tissue_exp=open("%stemp_%s"%(temp,chromnum),"w")

    with gzip.open ("/project2/nobrega/grace/expression/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz"%term) as expression:
        expression_colnames=expression.readline().strip().split()

    #Grab column with ENSG gene names, put it into first row of output
    cmd="""zcat /project2/nobrega/grace/expression/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz | awk 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i }} { print $(f["gene_id"]) }'"""%(term)
    status, output = commands.getstatusoutput(cmd)
    ENSG_ids=list(output.split('\n'))
    tissue_exp.write('\t'.join(ENSG_ids)+'\n')

    print "Grabbing %s expression data from GTEx v8..."%tissue
    for subject_ID in geno_IDs:
        #Write subject IDs into text file
        if subject_ID in expression_colnames:
            if write_IDs==True:
                subject_IDs.write(subject_ID+'\t'+subject_ID+'\n')
            col=int(expression_colnames.index(subject_ID))+1
            cmd="""zcat /project2/nobrega/grace/expression/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz | awk 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i }} { print $(f["%s"]) }'"""%(term,subject_ID)
            status, output = commands.getstatusoutput(cmd)
            subj_exp_dat=list(output.split('\n'))
            tissue_exp.write('\t'.join(subj_exp_dat)+'\n')

    if write_IDs==True:
        subject_IDs.close()
    tissue_exp.close()

    print "Combining expression data with TSS and gene name data..."
    #Get TSS and gene name from GTF, combine it with expression info
    gtf=pd.read_csv("/project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz",sep='\t',compression='gzip',header=None,skiprows=range(0,5))
    gtf.columns=['chr','source','type','start','stop','dot','strand','zeroes','info']
    gtf['ENSG']=[f.split('"')[1] for f in [e[0] for e in gtf['info'].str.split(';')]]
    gtf=gtf[gtf['type']=="gene"]
    gtf=gtf[gtf['info'].str.contains('protein_coding')]

    #Import gene expression
    gene_exp=pd.read_csv("%stemp_%s"%(temp,chromnum),sep='\t',index_col=0).T
    gene_exp['ENSG']=gene_exp.index

    #Get TSS, considering strand
    gtf['coord'] = np.where(gtf['strand'] == '+', gtf['start'],gtf['stop'])

    #get chr in '1' format
    gtf['chrom']=[e[1] for e in gtf['chr'].str.split('r',n=1)]
    gtf=gtf[gtf['chrom']!='M']
    gtf=gtf[gtf['chrom']==str(chromnum)]

    #get gene name from info
    info=gtf['info'].str.split(';')
    genes=list()
    for line in info:
        ind=[line.index(i) for i in line if 'gene_name' in i]
        if len(ind)==0:
            gene_name="NA"
        else:
            gene_name=line[ind[0]].split('"')[1]
        genes.append(gene_name)

    gtf['gene']=genes
    gtf = gtf[gtf.gene != "NA"]

    #Take only the columns of interest
    gtf=gtf[['chrom','coord','gene','ENSG']]
    merged=pd.merge(gtf,gene_exp,on="ENSG")
    merged=merged.drop(columns=['ENSG'])
    print "Writing output..."
    merged.to_csv("%s/%s_GTEx_v8_exp_%s.txt.gz"%(tissue,tissue,chromnum),sep='\t',index=False,compression='gzip')

    os.remove("%stemp_%s"%(args.temp,chromnum))

if data_source == "CMC":
    raw=pd.read_csv("%s/%s_%s_exp.txt.gz"%(tissue,tissue,data_source),sep=' ',index_col=0)
    raw=raw[list(set(geno_IDs) & set(raw.columns))]
    nsubjs=raw.shape[1]

    #Get chrom and coord for each gene in gene list
    genes=pd.read_csv("/project2/nobrega/grace/expression/coding_genes_hg19.bed",sep='\t',header=None,index_col=3)
    genes_info=pd.DataFrame(index=genes.index, columns=["chr","coord","gene"])
    genes_info['chr']=[re.sub("chr","",c) for c in list(genes.iloc[:,0])]
    coords=[0]*genes.shape[0]

    for i in range(0,genes.shape[0]):
        row=list(genes.iloc[i,:])
        strand=row[4]
        if strand=='+':
            coords[i]=row[1]
        if strand=='-':
            coords[i]=row[2]
            
    genes_info['coord']=coords
    genes_info['gene']=genes_info.index
    raw_coords=genes_info.merge(raw,left_index=True,right_index=True)

    exp=raw_coords.copy()

    #Select chromosome of interest
    exp=exp.loc[exp['chr'] == str(chromnum)]
    exp.to_csv("%s/%s_CMC_exp_%s.txt.gz"%(tissue,tissue,chromnum),sep='\t',index=False,compression='gzip')

    #Write subject IDs if they don't exist
    if not os.path.exists("%s/%s_subject_IDs"%(tissue,tissue)):
        IDlist=list(exp.columns.values)[3:]
        subj_IDs = pd.DataFrame(index=range(0,len(IDlist)),columns=['A','B'])
        subj_IDs['A']=IDlist
        subj_IDs['B']=IDlist
        subj_IDs.to_csv("%s/%s_subject_IDs"%(tissue,tissue),sep='\t',index=False,header=False)