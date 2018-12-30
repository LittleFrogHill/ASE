# Pipeline of analysing ASE genes of hybrid animals

## 1.Clean data
  resequencing data
  
    java -Xmx30g -jar /home/yjiang/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 SRR1868837_1.fastq SRR1868837_2.fastq SRR1868837_1.clean.fq.gz SRR1868837_1.unpaired.fq.gz SRR1868837_2.clean.fq.gz SRR1868837_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 >SRR1868837.log &
  
  RNA-seq
  
    java -Xmx30g -jar /home/yjiang/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 SRR1868837_1.fastq SRR1868837_2.fastq SRR1868837_1.clean.fq.gz SRR1868837_1.unpaired.fq.gz SRR1868837_2.clean.fq.gz SRR1868837_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 >SRR1868837.log &
    
## 2.Resequencing data mapping

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
        
  extract the unmapped reads
  
    #!/bin/sh
      for i in ENS_L2_I012 \
      EHM_L2_I021 \
      EFB_L4_I046 \
      DNS_L2_I011 \
      DHM_L2_I010 \
      DFB_L4_I045 \
      CNS_L6_I027 \
      CHM_L2_I008 \
      CFB_L4_I044 \
      BNS-1_L2_I007 \
      BHN1_L6_I026 \
      BFB_L4_I041 \
      ANS-1_L2_I006 \
      AHM-1_L2_I005 \
      AFB_L4_I042 \
      
      do
      
      samtools view -f 0x4 \
              /home/yjiang/liming/mule_hinny/analysis/mule.mapping/$i/$i.map.sam \
              > /home/yjiang/liming/mule_hinny/analysis/mule_sam_fastq/$i.unmap.sam
      
      done
      
umapped.sam2fastq

    #!/bin/sh

    java -Xmx4g -jar /home/yjiang/liming/bin/picard-tools-1.119/SamToFastq.jar \
        INPUT=/home/yjiang/liming/mule_hinny/analysis/mule_sam_fastq/BHN1_L6_I026.unmap.sam \
        FASTQ=BHN1_L6_I026_1.fastq \
        SECOND_END_FASTQ=BHN1_L6_I026_2.fastq
        
Tophat mapping
      
    for i in 1 \
    2 \
    4 \
    5 \
    6 \
    9 \
    10 \
    11 \
    12 \
    13 \
    AFB_L4_I042 \
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
    
    mkdir /home/yjiang/liming/mule_hinny/analysis/tophat/$i/
      
    cd /home/yjiang/liming/mule_hinny/analysis/tophat
    
    tophat --read-mismatches 10 \
            --read-edit-dist 10 \
            --output-dir /home/yjiang/liming/mule_hinny/analysis/tophat/$i \
            --GTF /home/yjiang/database/horse/Equus_caballus.EquCab2.79.gtf \
            --no-convert-bam \
            /home/yjiang/database/horse/bowtie2.horse2 \
            /home/yjiang/liming/mule_hinny/analysis/mule_sam_fastq/Last_fastq/$i\_1.fastq \
            /home/yjiang/liming/mule_hinny/analysis/mule_sam_fastq/Last_fastq/$i\_2.fastq
    done

## 5.Divergent sites 
  
    python get_allelic_number_based_on_SNP.py

## 6.Generate the speices info

    python generate_reads_info_from_sam_4.0.py

## 7.Distingguish the bam

    python distingguish_species_sam.py

## 8.Call ASE genes

    python ASEgenesCaller.py
  
## 9.Merge all the samples

    python join_ASE_result.py
  
  
## 10.Quantitate the expression by Stringtie

    stringtie -p 8 -G ~/genome.data/cattle/Cattle.GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff -o /stor9000/apps/users/NWSUAF/2015060145/StringTie/cattle/stringtie.gtf/${i}.gtf -l ${i} /stor9000/apps/users/NWSUAF/2015060145/total_bam/cattle/genome/hisat_merge_tophat.bam/sort.bam/${i}.total.sort.bam
  
    stringtie --merge -p 8 -G ~/genome.data/cattle/Cattle.GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff -o stringtie_merged.gtf cattle.mergelist.txt
  
    stringtie -e -B -p 8 -G ../stringtie_merged.gtf -o ~/StringTie/cattle/stringtie2ballgown/${i}/${i}.gtf /stor9000/apps/users/NWSUAF/2015060145/distinguish_species_bam/cattle/merge/sort.bam/${i}.sort.bam
  
## 11.Calculate the ASE genes by Ballgown in R

      library(ballgown)
      
      library(RSkittleBrewer)
      
      library(genefilter)
      
      library(dplyr)
      
      library(devtools)
      
      bg=ballgown(dataDir = "ballgown.MF", samplePattern = "ASE")
      
      pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0),8))
      
      bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
      
      results_transcripts = stattest(bg_filt,feature="transcript",covariate="group",getFC=TRUE, meas="FPKM")
      
      results_genes = stattest(bg_filt, feature="gene",covariate="group",getFC=TRUE,meas="FPKM")
      
      results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
      
      results_transcripts = arrange(results_transcripts,pval)
      
      results_genes = arrange(results_genes,pval)
      
      transcripts.sig<-subset(results_transcripts,results_transcripts$qval<0.05)
  
      gene.sig<-subset(results_genes,results_genes$qval<0.05)
    
      write.table(file="horse.all.ASE.gene.result",results_genes,quote=F,sep="\t",row.names=F,col.names=T)
    
      write.table(file="horse.all.ASE.gene.sig.result",gene.sig,quote=F,sep="\t",row.names=F,col.names=T)
    
      write.table(file="horse.all.ASE.transcripts.result",results_transcripts,quote=F,sep="\t",row.names=F,col.names=T)
    
      write.table(file="horse.all.ASE.transcripts.sig.result",transcripts.sig,quote=F,sep="\t",row.names=F,col.names=T)
    
      whole_tx_table = texpr(bg, 'all')

## 12.Alternative splcing deceting by rMATs

    python ~/rMATS.4.0.1/rMATS-turbo-Linux-UCS2/rmats.py --b1 alt --b2 ref --gtf /stor9000/apps/users/NWSUAF/2012015204/splicing/ensembl_gtf/horse.NCBI.annotation.gtf --od result -t paired --nthread 8 --readLength 125

## 13.Plot the alternative splcing events

    rmats2sashimiplot --b1 bam1,bam2,bam3 --b2 bam3,bam4,bam5 -e event.txt --l1 donkey --l2 horse --exon_s 1 --intron_s 5 -o out -t SE




