#!/usr/bin/env python3.4

'''
#-------------------------------------------------------------------------------
#Author:Yu Wang(wang_yu@nwsuaf.edu.cn)
#Time:  2016/1/20
#Version: 1.0
#-------------------------------------------------------------------------------
'''
#######################modle####
import sys
import re
import getopt
import os
import gc
def usage():
    print('''Useage: python script.py [option] [parameter]
    -s/--snp_file        input the snp file
    -b/--bam_file        input the bam/sam file
    -g/--gtf_file        input the genes annotation files 
    -o/--output          the output results file
    -r/reads_info        output the med file reads info 
    -m/med_file          output the med file reads number
    -h/--help            show possible options''')
#######################default#####
####################### argv ######
opts, args = getopt.getopt(sys.argv[1:], "hs:b:g:o:r:m:",["help","snp_file=","bam_file=","gtf_file=","output=","reads_info=","med_file="])
for op, value in opts:
    if op == "-s" or op == "--snp_file":
        snp_file = value
    elif op == "-b" or op == "--bam_file":
        bam_file = value
    elif op == "-g" or op == "--gtf_file":
        gtf_file = value    
    elif op == "-o" or op == "--output":
        output = value
    elif op == "-r" or op == "--reads_info":
        reads_info = value
    elif op == "-m" or op == "--med_file":
        med_file = value
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) < 9:
    usage()
    sys.exit(1)
f1=open(snp_file)
f2=os.popen('samtools view '+bam_file)
f3=open(gtf_file)
f4=open(output,'w')
f5=open(reads_info,'w')
f6=open(med_file,'w')

#load snp dictionary########################
'''
1       612     T       C
1       638     A       C
1       681     G       C
1       1596    T       C
'''
#load sam file ########################
'''
HWI-D00524:32:C4D07ACXX:6:1101:2199:1984        83      2       39256853        255     100M    =       39256777        -176    GTTTCTAAGGTCAAGCTGGTCTAGATCCATACATTTCCTCCAAGGCAAAATGTACCCATGACTTTTTTGTG
HWI-D00524:32:C4D07ACXX:6:1101:2199:1984        163     2       39256777        255     100M    =       39256853        176     AGCGGCCCTAGGAGTCACCTCAGCACCTTCTCCCCACCCTCCGACGCCACCTCCCTCTGGGGCTTGCTTGC
'''
#load gtf file ########################
'''
1       ensembl gene    11193   15975   .       +       .       gene_id "ENSECAG00000012421"; gene_version "1"; gene_name "SYCE1"; gene_source "ensembl"; gene_biotype "protein_coding";
1       ensembl transcript      11193   15975   .       +       .       gene_id "ENSECAG00000012421"; gene_version "1"; transcript_id "ENSECAT00000013004"; transcript_version "1"; gene_name "SYCE1";
'''
ref_dict={}
alt_dict={}
genome_dict={}
for snp in f1:
    snp=snp.split()
    index=snp[0]+'-'+snp[1]
    ref_dict[index]=snp[2]
    alt_dict[index]=snp[3]
    genome_dict[snp[0]]=genome_dict.get(snp[0],0)+1
###############define ########################

def get_region():
    '''
    reads_region_ref: the reads in which chromosome
    reads_region_start: the reads's start in chromosome
    reads_region_end: the reads's end in chromosome
    '''
    reads_region_ref[reads[0]]=reads[2]
    if re.search('-',reads[8]):
        reads_region_end[reads[0]]=int(reads[3])+reads_length
    else:
        reads_region_start[reads[0]]=int(reads[3])

def decide_ref_or_alt_sub(i,sub_reads_info):
    snp_index=reads[2]+'-'+str(int(reads[3])+i)
    if snp_index in ref_dict.keys():
        if sub_reads_info==ref_dict[snp_index]:
            reads_ref[reads[0]]=reads_ref.get(reads[0],0)+1
        elif sub_reads_info==alt_dict[snp_index]:
            reads_alt[reads[0]]=reads_alt.get(reads[0],0)+1
        else:
            pass
    else:
        pass
    
def decide_ref_or_alt():
    '''
    in this function, we match 7 different modle.
    for example:
        125M
        10M5N115M
        10M5N5M10N110M
        107M2D18M
        112M1I12M
        6M1I63M2I55M
        31M1D54M1D40M
        
    '''
    if re.search('^(\d+)M$',reads[5]):  ###125M###
        for i in range(reads_length):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)N(\d+)M$',reads[5]): ###10M5N115M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$',reads[5]):  ###10M5N5M10N110M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3])+int(pos_number[4])):
            sub_reads_info=reads[9][i-int(pos_number[1])-int(pos_number[3])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)D(\d+)M$',reads[5]):     ###107M2D18M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)I(\d+)M$',reads[5]): ###112M1I12M##
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0]),int(pos_number[0])+int(pos_number[2])):
            sub_reads_info=reads[9][i+int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M$',reads[5]): ###31M1D54M1D40M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3])+int(pos_number[4])):
            sub_reads_info=reads[9][i-int(pos_number[1])-int(pos_number[3])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M$',reads[5]): ###31M1I54M1I40M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0]),int(pos_number[0])+int(pos_number[2])):
            sub_reads_info=reads[9][i+int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[2]),int(pos_number[0])+int(pos_number[2])+int(pos_number[4])):
            sub_reads_info=reads[9][i+int(pos_number[1])+int(pos_number[3])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)S(\d+)M$',reads[5]):    ###32S118M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0]),reads_length):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)S$',reads[5]):  ####136M14S##
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)S(\d+)M(\d+)S$',reads[5]):###7S127M16S###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0]),int(pos_number[0])+int(pos_number[1])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    else:
        pass                       
################## define dictionary #############
reads_ref={}
reads_alt={}
reads_name={}
reads_region_ref={}
reads_region_start={}
reads_region_end={}
readsdic={}
sitedicalt={}
sitedicref={}
gene_ref={}
gene_alt={}
genedic={}
genedic1={}
geneinfo={}
genename={}
genelength={}
ASEgenes={}
##################################################
for reads in f2:
    reads=reads.split()
    try:
        if int(reads[4]) >= 50:
            if reads[2] in genome_dict:
                reads_length=len(reads[9])
                decide_ref_or_alt()
                get_region()
            else:
                pass
        else:
            pass
    except:
        pass  
#################### put all keys into reads_ref ###################
del ref_dict
gc.collect()
del alt_dict
gc.collect()
for key,value in reads_alt.items():
    if key in reads_ref.keys():
        pass
    else:
        reads_ref[key]=0
############################## genarate distinguish species reads dictionary ################
for key,value in reads_ref.items():
    try:
        readsdic[key]=str(key)+'\t'\
                +str(reads_region_ref[key])+'\t'\
                +str(reads_region_start.get(key,(reads_region_end.get(key,0)-125)))+'\t'\
                +str(reads_region_end.get(key,(reads_region_start.get(key,0)+125)))+'\t'\
                +str(reads_ref.get(key,0))+'\t'\
                +str(reads_alt.get(key,0))
        f5.write(str(key)+'\t'\
                +str(reads_region_ref[key])+'\t'\
                +str(reads_region_start.get(key,(reads_region_end.get(key,0)-125)))+'\t'\
                +str(reads_region_end.get(key,(reads_region_start.get(key,0)+125)))+'\t'\
                +str(reads_ref.get(key,0))+'\t'\
                +str(reads_alt.get(key,0))+'\n')
    except KeyError:
        pass
f5.close()       
###################### calculate distingguish species reads number ###############################        
species_info_reads_number=int(len(readsdic))
del reads_alt
gc.collect()
del reads_ref
gc.collect()
############################## gene select reads count ################

for key,value in readsdic.items():
    reads=value.split()
    if reads[4]=='0':
        sitedicalt[reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]]=''
    elif reads[5]=='0':
        sitedicref[reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]]=''
    else:
        pass

for gene in f3:
    gene=gene.split()
    try:
        if gene[2]=='gene':
            geneinfo[gene[0]+'\t'+gene[3]+'\t'+gene[4]]=gene[9].strip('"";')+'\t'+gene[13].strip('"";')
            gene_alt[gene[0]+'\t'+gene[3]+'\t'+gene[4]]=0
            gene_ref[gene[0]+'\t'+gene[3]+'\t'+gene[4]]=0
            a=[]
            b=[]
            a=[gene[3],gene[4]]
            genename[gene[3]+'\t'+gene[4]]=[]
            if gene[0]+'\t'+str(round(int(gene[3])/1000000)) in genedic:
                genedic[gene[0]+'\t'+str(round(int(gene[3])/1000000))].append(a)
            else:
                genedic[gene[0]+'\t'+str(round(int(gene[3])/1000000))]=[]
                genedic[gene[0]+'\t'+str(round(int(gene[3])/1000000))].append(a)
            if gene[0]+'\t'+str(round(int(gene[4])/1000000)) in genedic1:
                genedic1[gene[0]+'\t'+str(round(int(gene[4])/1000000))].append(a)
            else:
                genedic1[gene[0]+'\t'+str(round(int(gene[4])/1000000))]=[]
                genedic1[gene[0]+'\t'+str(round(int(gene[4])/1000000))].append(a)
        elif gene[2]=='exon':
            length=int(gene[4])-int(gene[3])
            a=gene[9].strip('"";')
            genelength[a]=genelength.get(a,0)+length
        else:
            pass
    except IndexError:
        pass

for key,value in sitedicref.items():
    key=key.split()
    try:
        a=genedic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=genedic1[key[1]+'\t'+str(round(int(key[3])/1000000))]        
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass 
    except KeyError:
        pass
for key,value in sitedicalt.items():
    key=key.split()
    try:
        a=genedic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=genedic1[key[1]+'\t'+str(round(int(key[3])/1000000))]
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass 
    except KeyError:
        pass
del sitedicref
gc.collect()
del sitedicalt
gc.collect()
for key,value in gene_alt.items():
    if key in gene_ref.keys():
        pass
    else:
        gene_ref[key]=0
###################### calculate distingguish species reads number ###############################
for key,value in gene_ref.items():
    ASEgenes[key]=str(key)+'\t'+geneinfo[key]+'\t'\
            +str(gene_ref.get(key,0))+'\t'\
            +str(gene_alt.get(key,0))+'\t'\
            +str(gene_ref.get(key,0)+gene_alt.get(key,0))
    f6.write(str(key)+'\t'+geneinfo[key]+'\t'\
            +str(gene_ref.get(key,0))+'\t'\
            +str(gene_alt.get(key,0))+'\t'\
            +str(gene_ref.get(key,0)+gene_alt.get(key,0))+'\n')
f6.close()
del gene_ref
gc.collect()
del gene_alt
gc.collect()
for key,value in ASEgenes.items():
    reads=value.split()
    ref=int(reads[5])*1000000000/(int(species_info_reads_number)*float(int(genelength[reads[3]])))
    alt=int(reads[6])*1000000000/(int(species_info_reads_number)*float(int(genelength[reads[3]])))
    total=int(reads[7])*1000000000/(int(species_info_reads_number)*float(int(genelength[reads[3]])))
    f4.write(reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]+'\t'+reads[4]+'\t'+str(round(ref,3))+'\t'+str(round(alt,3))+'\t'+str(round(total,3))+'\n')
        
f1.close()
f2.close()
f3.close()
f4.close()


