#!/usr/bin/env python
# -*- coding: ASCII -*-

#####This python script applyies to allele frequency estimator for pool seq data from Lynch et al. 2014 GBE.
#####Please note it is slow, and is only designed for obtain allele frequencies at specific SNPs in smallish numbers (tens to hundreds) of individuals. I do not recommend using this as a general variant caller.
#####Pysam and Numpy must be installed in the version of python used.
#####Bam files must be indexed.
#####To guard against mis-mappings and CNVs, the program only outputs sites with coverage between a third and twice the mean at the set of SNPs considered
#####The SNP file must have the tab seperated fields in the following order: chromosome, position (one-based), reference allele, alternate allele
#####If you want to change things like mapping and base quality threshold, edit the python code under the section "Input arguments"
#####Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)

#####usage is ./LynchPool_called.py <bamfile> <fileoutname> <target_SNP_file> <reference_genome>


###import libraries
import string
import numpy as np
import pysam
import gzip
import math
import copy
from sys import argv
import time

###Input arguments
BAMin=argv[1] #filename of bam with poolseq data
filenameout=argv[2] #creates a vcf
SNPfile=argv[3] #must have the tab seperated fields in the following order chromosome, position (one-based), reference allele, alternate allele
#ref_file=argv[4]
MQ_t=20 #mapping quality threshold
BQ_t=20 #base_qualitythreshold

###converts phred score to probability
def phred2prob(x):
    return 10.0**(-x/10.0)

###converts probability to phred score
def prob2phred(x):
    return -10*math.log10(x)

###extract base of reads for a give position in a bam.
def extract_bam_SNP_base_only(samfile,chromo,pos,BQ,MQ):
    var_list=[]
    for pileupcolumn in samfile.pileup(chromo,pos-1,pos,truncate=True,stepper='all'):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ):
                        var_list.append(pileupread.alignment.query_sequence[pileupread.query_position])
    return var_list

#####open reference file
##ref=pysam.FastaFile(ref_file)
##chromos=ref.references

###Read SNPlist
file=open(SNPfile,'r')
SNPs=file.read()
SNPs=string.split(SNPs,'\n')
if SNPs[-1]=='':
    del(SNPs[-1])

nb_SNPs=len(SNPs)


samfile = pysam.AlignmentFile(BAMin, "rb")
samp_name=samfile.header['RG'][0]['SM']

all_counts=np.zeros((nb_SNPs,3),dtype='float32')  #major_allele, minor_allele, other_allele
Mm_allele=np.zeros((nb_SNPs),dtype='int32')   #0=ref is major allele, 1=alt is major allele. In a tie, ref is chosen as major

###Iterate through SNPs
for g in range(len(SNPs)):
    k=string.split(SNPs[g])
    chromo=k[0]
    k[1]=int(k[1])
    SNPs[g]=k
    ref_all=k[2]
    alt_all=k[3]
    
    read_list=extract_bam_SNP_base_only(samfile,k[0],k[1],BQ_t,MQ_t)
    #all_counts[g]=[read_list.count('A'),read_list.count('C'),read_list.count('G'),read_list.count('T')]
    if read_list.count(ref_all)>=read_list.count(alt_all):
        all_counts[g]=[read_list.count(ref_all),read_list.count(alt_all),len(read_list)-read_list.count(ref_all)-read_list.count(alt_all)]
    else:
        all_counts[g]=[read_list.count(alt_all),read_list.count(ref_all),len(read_list)-read_list.count(ref_all)-read_list.count(alt_all)]
        Mm_allele[g]=1

    if g%100==0:
        print g,k



###calulate depth stats
depth=np.sum(all_counts,axis=1)
min_DP=round(np.mean(depth)/3)
max_DP=round(np.mean(depth)*2)


###first iteration through Lynch
#calculate error
p_e=all_counts[:,2]/depth
E=3*p_e/2

p_m=all_counts[:,0]/np.sum(all_counts[:,0:2],axis=1)

p_hat=p_m*(1.0-(2.0*E/3.0))-(E/3.0)/1-(4*E/3)

###second iteration
a=p_hat>0.9
b=depth<min_DP
c=depth>max_DP

exclude=a+b+c
exclude_swap=exclude == False

E_mean=np.mean(E[exclude_swap])

for g in range(len(E)):
    if p_hat[g]>0.9:
        E[g]=E_mean

p_hat2=p_m*(1.0-(2.0*E/3.0))-(E/3.0)/1-(4*E/3)



outfile=open(filenameout,'w')

header='chrom\tpos\tref\talt\talt_AF\talt_AF_correct\tref:alt:other_dp\tDP\tDP_ok?\n'
outfile.write(header)

for g in range(len(SNPs)):
    out=SNPs[g][0]+'\t'+str(SNPs[g][1])+'\t'+SNPs[g][2]+'\t'+SNPs[g][3]+'\t'

    if Mm_allele[g]==0:
        out=out+str(1-round(p_hat[g],3))+'\t'+str(1-round(p_hat2[g],3))+'\t'+str(int(all_counts[g][0]))+':'+str(int(all_counts[g][1]))
    else:
        out=out+str(round(p_hat[g],3))+'\t'+str(round(p_hat2[g],3))+'\t'+str(int(all_counts[g][1]))+':'+str(int(all_counts[g][0])) 

    out=out+':'+str(int(all_counts[g][2]))+'\t'+str(int(depth[g]))+'\t'

    if min_DP<=depth[g]<=max_DP:
        out=out+'Y\n'
    else:
        out=out+'N\n'

    outfile.write(out)

outfile.close()
