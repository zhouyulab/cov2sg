#!/bin/bash
# Follow mk_trs_ex.py generated trsread.fa
mkdir -p wivSTAR
STAR --runMode genomeGenerate --genomeDir wivSTAR --genomeFastaFiles WIV04.fasta

mkdir -p trsmapped
STAR --runThreadN 1 \
    --genomeDir wivSTAR \
    --readFilesIn trsread.fa \
    --outFileNamePrefix trsmapped/ \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.07 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 4

samtools sort -O BAM -o trsread.bam trsmapped/Aligned.out.sam
samtools index trsread.bam

