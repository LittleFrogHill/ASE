import sys
import os
if len(sys.argv) != 5:
    print('Useage:' 'python_script.py''<input_sam_file>' '<input_reads_species_info_file>''<output_file_ref>' '<output_file_alt>')
    sys.exit(1)
f1=open(sys.argv[1])
f2=open(sys.argv[2])

'''
####load sam f1###
HWI-7001457:246:C79W4ANXX:2:1101:14832:2309     339     27      37368206        1       125M    =       37368047        -284    GCGAGGGCGGTGAGTTTTATGTAAGTAGGTATGGTTATTTCTGGGACGGTTGTGGGGTAGATATTGTTGGAGATGAAGAATCCGGCAAAAATGCTGCCAATT
HWI-7001457:246:C79W4ANXX:2:1101:14832:2309     419     27      37368047        1       125M    =       37368206        284     CTTATTGATAGGTTAGCGAGCGGTGGGAGTCGGTGTATAATTGTTGGGTAGTATCCTAGGAGGTTGGAGAATTTGAATACGCTGGTGGAGTGTTCTAGTTTT
HWI-7001457:246:C79W4ANXX:2:1101:14769:2263     77      *       0       0       *       *       0       0       TAGGAGGAAAATGTTTCCCAAAAGACTCTAATTTTCACTTGATATATTTCCAAGATGAAAATCTTCTAAATATCAAAATCTTTCTTAAACCACACACTCATTGTATGACTTTGTACAG
HWI-7001457:246:C79W4ANXX:2:1101:14769:2263     141     *       0       0       *       *       0       0       CGGCAGTGGCAGCGGGAGCAGCTGCGCCGCCAGCGTGGAGCGGGAGGCGTGTGGGTGTGGAGAGAAGGGAGCGCCGGACCGGAACGGATGTGTACATCTCATGCTTGATGTGAAGAAA
HWI-7001457:246:C79W4ANXX:2:1101:14793:2291     99      MT      13152   1       125M    =       13271   244     TCTTAATTGGCAGCATTTTTGCCGGATTCTTCATCTCCAACAATATCTACCCCACAACCGTCCCAGAAATAACCATACCTACTTACATAAAACTCACCGCCCTCGCAGTAACCATCCT
#### load reads count file f2###
HWI-7001457:246:C79W4ANXX:2:2216:1528:44977     26      35498819        35499050        4       0
HWI-7001457:246:C79W4ANXX:2:2309:2519:89849     5       9257700 9257928 1       0
HWI-7001457:246:C79W4ANXX:2:2114:3383:52075     12      33034824        33035075        1       0
HWI-7001457:246:C79W4ANXX:2:1103:10384:70864    1       92340920        92341159        0       2
HWI-7001457:246:C79W4ANXX:2:1203:8835:98972     2       56472223        56472508        0       1
HWI-7001457:246:C79W4ANXX:2:2304:15605:9808     30      28465488        28467111        2       0
HWI-7001457:246:C79W4ANXX:2:1215:14218:79456    22      34800418        34800777        2       0
'''
ref={}
alt={}
reads_ref={}
reads_alt={}

for readscount in f2:
    readscount=readscount.split()
    if int(readscount[5])==0:
        ref[readscount[0]]=''
    elif int(readscount[4])==0:
        alt[readscount[0]]=''
    else:
        pass
os.system("samtools view -H "+sys.argv[1]+">"+sys.argv[3])
os.system("samtools view -H "+sys.argv[1]+">"+sys.argv[4])
for reads in f1:
    reads=reads.split()
    try:
        if reads[0] in ref:
            info='\t'.join(reads)
            reads_ref[reads[0]]=reads_ref.get(reads[0],'')+'\n'+info
        elif reads[0] in alt:
            info='\t'.join(reads)
            reads_alt[reads[0]]=reads_alt.get(reads[0],'')+'\n'+info
        else:
            pass
    except IndexError:
        pass
f3=open(sys.argv[3],'a')
f4=open(sys.argv[4],'a')        
for key,value in reads_ref.items():
    f3.write(reads_ref[key].strip('^\n')+'\n')
for key,value in reads_alt.items():
    f4.write(reads_alt[key].strip('^\n')+'\n')
f1.close()
f2.close()
f3.close()
f4.close()
