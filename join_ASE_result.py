#!/usr/bin/env python3.4
import sys
import getopt
def usage():
    print('''Useage: python script.py [option] [parameter]
    -l/--ASE_list            input the ASE result file list
    -o/--results             the  results
    -c/--classes               the class gene or lncRNA (defout:gene)
    -h/--help                show possible options''')
#######################default
classes='gene'
opts, args = getopt.getopt(sys.argv[1:], "hl:c:o:",["help","ASE_list=","classes=","results="])
for op, value in opts:
    if op == "-l" or op == "--ASE_list":
        ASE_list = value
    elif op == "-c" or op == "--classes":
        classes = str(value)
    elif op == "-o" or op == "--results":
        results = value
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) < 5:
    usage()
    sys.exit(1)


f1=open(ASE_list)
f2=open(results,'w')

''' lncRNA ASE result #################
XLOC_022690+lncRNA      5       75966672        75967600        3       374     0.0     0.581   0.581
XLOC_000812+lncRNA      1       88669153        88675534        2       1303    6.837   3.502   10.339
XLOC_026382+lncRNA      8       51411184        51412130        1       945     1.609   0.23    1.839
XLOC_012673+lncRNA      2       58794184        58795896        1       1711    13.715  10.159  23.874
XLOC_017456+lncRNA      26      22936374        22937356        1       981     1.772   3.544   5.316

/home/yjiang/liming/yak/cattle_ref/00.shell/yak_total_data.sam/dzo.ASE_lncRNA_results/2.map.sam.test

load gene ASE result ################
19      29670744        29876650        ENSBTAG00000019107      GAS7    3.91    3.724   7.634
25      21080469        21108140        ENSBTAG00000010163      SCNN1G  0.702   0.2     0.902
3       57390507        57391002        ENSBTAG00000047364      ensembl 0.0     0.0     0.0
20      66662746        66662852        ENSBTAG00000043275      U6      0.0     0.0     0.0
11      93011815        93029730        ENSBTAG00000004295      NDUFA8  53.978  56.718  110.696

'''
if classes=='gene':
    f1=open(ASE_list)
    f2=open(results,'w')
    title='chr\tstart\tend\tgeneID\tname'
    name=''
    dict={}
    for file in f1:
        file=file.strip()
        a=file.split('/')
        b=a[-1].split('_')
        c=b[0]
        name+=c+'-ref'+'\t'+c+'-alt'+'\t'+c+'-sum'+'\t'
        file=open(file)    
        for reads in file:
            reads=reads.split()
            dict[reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]+'\t'+reads[4]]=dict.get(reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]+'\t'+reads[4],'')+'\t'+reads[5]+'\t'+reads[6]+'\t'+reads[7]
    f2.write(title+'\t'+name+'\n')
    for key,value in dict.items():
        f2.write(key+value+'\n')
    f1.close()
    f2.close()
elif classes=='lncRNA':
    f1=open(ASE_list)
    f2=open(results,'w')
    title='name\tchr\tstart\tend\texon\tlength'
    name=''
    dict={}
    for file in f1:
        file=file.strip()
        a=file.split('/')
        b=a[-1].split('.')
        c=b[0]
        name+=c+'-ref'+'\t'+c+'-alt'+'\t'+c+'-sum'+'\t'
        file=open(file)    
        for reads in file:
            reads=reads.split()
            dict[reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]+'\t'+reads[4]+'\t'+reads[5]]=dict.get(reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]+'\t'+reads[4]+'\t'+reads[5],'')+'\t'+reads[6]+'\t'+reads[7]+'\t'+reads[8]
    f2.write(title+'\t'+name+'\n')
    for key,value in dict.items():
        f2.write(key+value+'\n')
    f1.close()
    f2.close()
else:
    pass
