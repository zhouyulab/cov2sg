import os
include: "MetaInfo.smk"

RNA_NARRY_KIM = "analysis/RNA_seq_NarryKim"
rule RNA_cutadapt_NarryKim:
    input:
        fastq_R1 = "data/NarryKimRNAseq/VERO_SARS_CoV_2.R1.fastq.gz",
        fastq_R2 = "data/NarryKimRNAseq/VERO_SARS_CoV_2.R1.fastq.gz",
    output:
        trimed_fq_R1 = os.path.join(RNA_NARRY_KIM, "trim", "NarryKim.R1.fastq.gz"),
        trimed_fq_R2 = os.path.join(RNA_NARRY_KIM, "trim", "NarryKim.R2.fastq.gz"),
    params:
        adaptor_R1 = "AGATCGGAAGAGC",
        adaptor_R2 = "AGATCGGAAGAGC",
        min_read_len = 20,
    shell:
        "cutadapt --cores 0 -e 0.25 --match-read-wildcards -a {params.adaptor_R1} -A {params.adaptor_R2} -m {params.min_read_len} -o {output.trimed_fq_R1} -p {output.trimed_fq_R2} {input.fastq_R1} {input.fastq_R2}; fastqc {output.trimed_fq_R1} {output.trimed_fq_R2}"

rule RNA_rm_rRNA_NarryKim:
    input:
        fq_R1 = rules.RNA_cutadapt_NarryKim.output.trimed_fq_R1,
        fq_R2 = rules.RNA_cutadapt_NarryKim.output.trimed_fq_R2,
    output:
        unmap_R1 = os.path.join(RNA_NARRY_KIM, "rm_rRNA", "NarryKim.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_NARRY_KIM, "rm_rRNA", "NarryKim.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_NARRY_KIM, "rm_rRNA", "NarryKim"),
        indx = rRNA_INDX,
    threads: 64
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

rule RNA_STAR_NarryKim:
    input:
        fq_R1 = rules.RNA_rm_rRNA_NarryKim.output.unmap_R1,
        fq_R2 = rules.RNA_rm_rRNA_NarryKim.output.unmap_R2,
        gtf = REF_GTF,
    output:
        bam = os.path.join(RNA_NARRY_KIM, "STAR", "NarryKim.sort.bam"),
        unmap_R1 = os.path.join(RNA_NARRY_KIM, "STAR", "NarryKim.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_NARRY_KIM, "STAR", "NarryKim.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_NARRY_KIM, "STAR", "NarryKim"),
        indx = chkSab_INDX,
    threads: 64
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

rule RNA_STAR_WIV04_NarryKim:
    input:
        fq_R1 = rules.RNA_STAR_NarryKim.output.unmap_R1,
        fq_R2 = rules.RNA_STAR_NarryKim.output.unmap_R2,
    output:
        bam = os.path.join(RNA_NARRY_KIM, "STAR_WIV04_rm_host", "NarryKim.sort.bam"),
    params:
        mapping_dir = os.path.join(RNA_NARRY_KIM, "STAR_WIV04_rm_host", "NarryKim"),
        indx = WIV04_INDX,
    threads: 64
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
gzip -d -c {input.fq_R1} > {params.mapping_dir}/tmp.fastq
gzip -d -c {input.fq_R2} >> {params.mapping_dir}/tmp.fastq
STAR --runThreadN {threads} \
    --genomeDir {params.indx} \
    --outFileNamePrefix {params.mapping_dir}/ \
    --readFilesIn {params.mapping_dir}/tmp.fastq \
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
    --alignSJstitchMismatchNmax -1 -1 -1 -1
samtools sort -@ {threads} -o {output.bam} {params.mapping_dir}/Aligned.out.sam
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam
rm {params.mapping_dir}/tmp.fastq
        """

rule WIV04_bam2bw_NarryKim:
    input:
        bam = rules.RNA_STAR_WIV04_NarryKim.output.bam,
        size = WIV04_CHROM_SIZE,
    output:
        flag = touch(os.path.join(RNA_NARRY_KIM, "WIV04_bw", "NarryKim.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(RNA_NARRY_KIM, "WIV04_bw"),
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/NarryKim
#wigToBigWig -clip {params.work_dir}/NarryKim.Forward.wig {input.size} {params.work_dir}/NarryKim.Forward.bw
#wigToBigWig -clip {params.work_dir}/NarryKim.Reverse.wig {input.size} {params.work_dir}/NarryKim.Reverse.bw
#rm {params.work_dir}/NarryKim.Forward.wig {params.work_dir}/NarryKim.Reverse.wig
        """

rule cnt_JS_NarryKim:
    input:
        bam = rules.RNA_STAR_WIV04_NarryKim.output.bam,
        genome = WIV04_GENOME_FA,
    output:
        js = os.path.join(RNA_NARRY_KIM, "STAR_WIV04_JS", "NarryKim_js.bed"),
    params:
        rnatk_ja = "python scripts/subgenome/subgenome.py",
        outprefix = os.path.join(RNA_NARRY_KIM, "STAR_WIV04_JS", "NarryKim")
    shell:
        """
source activate py35
{params.rnatk_ja} --window 5 --mingap 30 --minevi 1 --armsize 20 --genome {input.genome} -o {params.outprefix} -v {input.bam}
        """

rule RNA_seq_NarryKim:
    input:
        rules.cnt_JS_NarryKim.output,
        rules.WIV04_bam2bw_NarryKim.output,
