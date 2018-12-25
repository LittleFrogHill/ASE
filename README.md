# Pipeline of analysing ASE genes of hybrid animals

## 1.trim.data Reads 质控
  resequencing data
  
    java -Xmx30g -jar /home/yjiang/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 SRR1868837_1.fastq SRR1868837_2.fastq SRR1868837_1.clean.fq.gz SRR1868837_1.unpaired.fq.gz SRR1868837_2.clean.fq.gz SRR1868837_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 >SRR1868837.log &
  
  RNA-seq
  
    java -Xmx30g -jar /home/yjiang/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 SRR1868837_1.fastq SRR1868837_2.fastq SRR1868837_1.clean.fq.gz SRR1868837_1.unpaired.fq.gz SRR1868837_2.clean.fq.gz SRR1868837_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 >SRR1868837.log &
    
## 2.resequencing data mapping

  bwa mapping

    bwa mem -M -t 10 -R "@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}" /home/yjiang/database/yeast/yeast ${i}_1.clean.fq.gz ${i}_2.clean.fq.gz > ${i}.sam

## 3.SNP calling

  sam2bam
  
    samtools view -bS SRR1868848.sam >SRR1868848.bam
    
  sort bam
  
    java -Xmx4g -jar /home/yjiang/liming/bin/picard-tools-1.119/SortSam.jar INPUT=SRR1868848.bam OUTPUT=SRR1868848.sort.bam SORT_ORDER=coordinate
    
  MarkDuplicates
  
    nohup java -Xmx4g -jar /home/yjiang/liming/bin/picard-tools-1.56/MarkDuplicates.jar REMOVE_DUPLICATES=false MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=800 INPUT=SRR1868845.sort.bam OUTPUT=SRR1868845.sort.dedup.bam METRICS_FILE=SRR1868845.sort.dedup.bam.matrix &
    
  Index bam
  
    samtools index Clean_FCHMTJ2CCXX_L_WH1603002611.sort.dedup.bam ./index/Clean_FCHMTJ2CCXX_L_WH1603002611.bai

  Realigner1

    java -Xmx16g -jar /home/yjiang/liming/bin/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -R /home/yjiang/database/yeast/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa -T RealignerTargetCreator -I ${i}.sort.dedup.bam -o ${i}.sort.dedup.bam.intervals
    
  Realigner2

    java -Xmx16g -jar /home/yjiang/liming/bin/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -R /home/yjiang/database/yeast/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa -T IndelRealigner  -targetIntervals ${i}.sort.dedup.bam.intervals -I ${i}.sort.dedup.bam -o ${i}.sort.dedup.realn.bam

  GATK_call_snp
  
    java -jar /home/yjiang/liming/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
        -allowPotentiallyMisencodedQuals \
        -R /home/yjiang/database/yeast/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868837.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868838.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868839.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868843.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868844.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868845.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868846.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868847.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868848.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868849.sort.dedup.realn.bam \
        -I /home/yjiang/liming/download_data/yeast/yeast.bam/SRR1868852.sort.dedup.realn.bam \
        -T UnifiedGenotyper \
        -o gatk.raw.vcf \
        -stand_call_conf 30.0 \
        -stand_emit_conf 0 \
        -glm SNP \
        -rf BadCigar
        
  Filter
  
    nohup java -jar /home/yjiang/liming/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -R /home/yjiang/database/yeast/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa -T VariantFiltration --filterExpression " QD < 20.0 || ReadPosRankSum < -8.0 || FS > 10.0 || QUAL < 30.0" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant ./gatk.raw.vcf -o yeast.filtered.raw.vcf

  Grep Filter  file

    grep -v "Filter" fruitfly.filtered1.raw.vcf > fruitfly1.final.raw.vcf &

    perl -lane 'next if /^#/; print $_ if length($F[4]) == 1' yeast.final.raw.vcf > yeast.final.raw.snp

    les yeast.final.raw.snp |awk '{print $1"\t"$2"\t"$4"\t"$5}' >yeast.final.snp
    
    
## 4.RNA_seq data mapping

  Hisat2 command：
  
    #!/bin/sh
      for i in AFB_L4_I042 \
        AHM-1_L2_I005 \
        ANS-1_L2_I006 \
        BFB_L4_I041 \
        BHM1_L6_I026 \
        BNS-1_L2_I007 \
        CFB_L4_I044 \
        CHM_L2_I008 \
        CNS_L6_I027 \
        DFB_L4_I045 \
        DHM_L2_I010 \
        DNS_L2_I011 \
        EFB_L4_I046 \
        EHM_L2_I021 \
        ENS_L2_I012
        
        do
        
        mkdir /home/yjiang/liming/mule_hinny/analysis/mule.mapping/$i/
        
        hisat --pen-noncansplice 1000000 \
                --no-mixed \
                --no-discordant \
                --un-conc /home/yjiang/liming/mule_hinny/analysis/mule.mapping/$i/ \
                -p 6 \
                --known-splicesite-infile /home/yjiang/database/horse/horse_splice_sites \
                -x /home/yjiang/database/horse/horse \
                -1 /home/yjiang/liming/mule_hinny/rna_seq/$i\.R1.clean.fastq.gz \
                -2 le's
                -S /home/yjiang/liming/mule_hinny/analysis/mule.mapping/$i/$i.map.sam
        
        done
