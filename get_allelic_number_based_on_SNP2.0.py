#!/usr/bin/env python3.4

'''
#-------------------------------------------------------------------------------
#Author:WangYu (wang_yu@nwsuaf.edu.cn)
#Time:  2016/3/10
#Version: 2.0
#useage: get allelic numbers based on SNP from sam file
#-------------------------------------------------------------------------------
'''
import os
import sys
import re
import getopt
def usage():
    print('''Useage: python script.py [option] [parameter]
    -s/--snp_file        input the snp file
    -b/--bam_file        input the bam/sam file
    -o/--output          the output results file
    -h/--help            show possible options''')
#######################default

opts, args = getopt.getopt(sys.argv[1:], "hs:b:o:",["help","snp_file=","bam_file=","output="])
for op, value in opts:
    if op == "-s" or op == "--snp_file":
        snp_file = value
    elif op == "-b" or op == "--bam_file":
        bam_file = value
    elif op == "-o" or op == "--output":
        output = value
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) < 7:
    usage()
    sys.exit(1)

f1=open(snp_file)
f2=os.popen('samtools view '+bam_file)
f3=open(output,'w')

#load snp dictionary########################
'''
1       612     T       C
1       638     A       C
1       681     G       C
1       1596    T       C
'''

ref_dict={}
alt_dict={}
genome_dict={}
allele_ref={}
allele_alt={}
allele_other={}
snpinfo={}
for snp in f1:
    snp=snp.split()
    index=snp[0]+'-'+snp[1]
    ref_dict[index]=snp[2]
    alt_dict[index]=snp[3]
    snpinfo[index]=snp[2]+'\t'+snp[3]
    genome_dict[snp[0]]=genome_dict.get(snp[0],0)+1
f1.close()
########################################

def decide_ref_or_alt_sub(i,sub_reads_info):
    snp_index=reads[2]+'-'+str(int(reads[3])+i)
    if snp_index in ref_dict.keys():
        if sub_reads_info==ref_dict[snp_index]:
            allele_ref[snp_index]=allele_ref.get(snp_index,0)+1
        elif sub_reads_info==alt_dict[snp_index]:
            allele_alt[snp_index]=allele_alt.get(snp_index,0)+1
        else:
            allele_other[snp_index]=allele_other.get(snp_index,0)+1
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
    else:
        pass                       
################## define dictionary #############
reads_name={}
reads_region_ref={}
reads_region_start={}
reads_region_end={}
##################################################


for reads in f2:
    reads=reads.split()
    try:
        if int(reads[4]) >= 50:
            if reads[2] in genome_dict:
                reads_length=len(reads[9])
                decide_ref_or_alt()
            else:
                pass
        else:
            pass
    except:
        pass  
#f3.write('reads_name\tchr\tstart\tend\tref\talt\n')


############################## output what we want ################
f4=open(snp_file)
for snp in f4:
    snp=snp.split()
    try:
        f3.write(snp[0]+'\t'+snp[1]+'\t'\
                +snpinfo[snp[0]+'-'+snp[1]]+'\t'\
                +str(allele_ref.get(snp[0]+'-'+snp[1],0))+'\t'\
                +str(allele_alt.get(snp[0]+'-'+snp[1],0))+'\t'\
                +str(allele_other.get(snp[0]+'-'+snp[1],0))+'\n')
    except KeyError:
        pass

###################### close file ###############################
f4.close()
f2.close()
f3.close()
