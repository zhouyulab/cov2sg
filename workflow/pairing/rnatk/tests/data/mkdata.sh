#!/bin/bash

bedtools bedpetobam -i dupbcpe.bedpe -g tk.chromInfo > tmp.bam
samtools sort -O BAM tmp.bam > dupbcpe.bam
samtools index dupbcpe.bam
samtools flagstat dupbcpe.bam
rm -f tmp.bam

