# Pipeline of analysing ASE genes of hybrid animals

##1.trim.data Reads 质控

  BWA mapping
    bwa mem -M -t 10 -R "@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}" \
    /home/yjiang/database/yeast/yeast ${i}_1.clean.fq.gz ${i}_2.clean.fq.gz > ${i}.sam
