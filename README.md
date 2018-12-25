# Pipeline of analysing ASE genes of hybrid animals

## 1.trim.data Reads 质控
  resequencing data
    java -Xmx30g -jar /home/yjiang/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 SRR1868837_1.fastq SRR1868837_2.fastq SRR1868837_1.clean.fq.gz SRR1868837_1.unpaired.fq.gz SRR1868837_2.clean.fq.gz SRR1868837_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 >SRR1868837.log &
  RNA-seq data
  
  
  resequencing data
    bwa mem -M -t 10 -R "@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}" /home/yjiang/database/yeast/yeast ${i}_1.clean.fq.gz ${i}_2.clean.fq.gz > ${i}.sam

  
    
