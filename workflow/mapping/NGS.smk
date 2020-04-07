import os
include: "MetaInfo.smk"

RNA_BASE = os.path.join("analysis", "RNA_seq")
rule RNA_cutadapt:
    input:
        fastq_R1 = os.path.join("data", "RNA_seq", "CoV_{rep}.R1.fastq.gz"),
        fastq_R2 = os.path.join("data", "RNA_seq", "CoV_{rep}.R2.fastq.gz"),
    output:
        trimed_fq_R1 = os.path.join(RNA_BASE, "trim", "CoV_{rep}.trim.R1.fastq.gz"),
        trimed_fq_R2 = os.path.join(RNA_BASE, "trim", "CoV_{rep}.trim.R2.fastq.gz"),
    params:
        job_name = "RNA_cutadapt",
        walltime = "480:00:00",
        nodes = "nodes=1:ppn=1",
        log = "log",
        queue = "wo_dell01",
        adaptor_R1 = "AGATCGGAAGAGC",
        adaptor_R2 = "AGATCGGAAGAGC",
        min_read_len = 20,
    shell:
        "cutadapt -e 0.25 --match-read-wildcards -a {params.adaptor_R1} -A {params.adaptor_R2} -m {params.min_read_len} -o {output.trimed_fq_R1} -p {output.trimed_fq_R2} {input.fastq_R1} {input.fastq_R2}; fastqc {output.trimed_fq_R1} {output.trimed_fq_R2}"


rule RNA_rm_rRNA:
    input:
        fq_R1 = rules.RNA_cutadapt.output.trimed_fq_R1,
        fq_R2 = rules.RNA_cutadapt.output.trimed_fq_R2,
    output:
        unmap_R1 = os.path.join(RNA_BASE, "rm_rRNA", "CoV_{rep}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_BASE, "rm_rRNA", "CoV_{rep}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "rm_rRNA", "CoV_{rep}"),
        indx = rRNA_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --genomeDir {params.indx} \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --outFileNamePrefix {params.mapping_dir}/ \
    --outReadsUnmapped Fastx \
    --readFilesCommand gunzip -c
gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
rm -r {params.mapping_dir}
        """

rule RNA_STAR:
    input:
        fq_R1 = rules.RNA_rm_rRNA.output.unmap_R1,
        fq_R2 = rules.RNA_rm_rRNA.output.unmap_R2,
        gtf = REF_GTF,
    output:
        bam = os.path.join(RNA_BASE, "STAR", "CoV_{rep}.sort.bam"),
        unmap_R1 = os.path.join(RNA_BASE, "STAR", "CoV_{rep}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_BASE, "STAR", "CoV_{rep}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "STAR", "CoV_{rep}"),
        indx = chkSab_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --quantMode TranscriptomeSAM \
    --genomeDir {params.indx} \
    --sjdbGTFfile {input.gtf} \
    --sjdbScore 1 \
    --outFileNamePrefix {params.mapping_dir}/ \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outReadsUnmapped Fastx \
    --readFilesCommand gunzip -c
samtools sort -@ {threads} -o {output.bam} {params.mapping_dir}/Aligned.out.sam
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam
gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
        """

rule RNA_STAR_WIV04:
    input:
        fq_R1 = rules.RNA_STAR.output.unmap_R1,
        fq_R2 = rules.RNA_STAR.output.unmap_R2,
    output:
        bam = os.path.join(RNA_BASE, "STAR_WIV04_rm_host", "CoV_{rep}.sort.bam"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "STAR_WIV04_rm_host", "CoV_{rep}"),
        indx = WIV04_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --genomeDir {params.indx} \
    --outFileNamePrefix {params.mapping_dir}/ \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --outFilterMultimapNmax 1 \
    --alignSJoverhangMin 8 \
    --outSJfilterOverhangMin 8 8 8 8 \
    --outSJfilterCountUniqueMin 3 3 3 3 \
    --outSJfilterCountTotalMin 3 3 3 3 \
    --outSJfilterDistToOtherSJmin 0 0 0 0 \
    --scoreGap -4 \
    --scoreGapNoncan -4 \
    --scoreGapATAC -4 \
    --alignIntronMax 30000 \
    --alignMatesGapMax 30000 \
    --alignSJstitchMismatchNmax -1 -1 -1 -1 \
    --readFilesCommand gunzip -c
samtools sort -@ {threads} -o {output.bam} {params.mapping_dir}/Aligned.out.sam
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam
        """

rule RNA_seq:
    input:
        expand(rules.RNA_STAR_WIV04.output, rep=REPs),

