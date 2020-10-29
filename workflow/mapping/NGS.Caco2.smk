import os
include: "MetaInfo.smk"

RNA_TP_BASE = os.path.join("analysis", "RNA_seq_caco2")
RNA_TP_SAMPLE = {"4h": ["rep1", "rep2"], "12h": ["rep1", "rep2"], "24h": ["rep1", "rep2"]}


rule RNA_cutadapt_TP:
    input:
        fastq_R1 = os.path.join("data", "RNA_seq_caco2", "Caco2.{tp}.{part}.R1.fastq.gz"),
        fastq_R2 = os.path.join("data", "RNA_seq_caco2", "Caco2.{tp}.{part}.R2.fastq.gz"),
    output:
        trimed_fq_R1 = os.path.join(RNA_TP_BASE, "trim", "WIV04.{tp}.{part}.trim.R1.fastq.gz"),
        trimed_fq_R2 = os.path.join(RNA_TP_BASE, "trim", "WIV04.{tp}.{part}.trim.R2.fastq.gz"),
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

rule RNA_rm_rRNA_TP:
    input:
        fq_R1 = rules.RNA_cutadapt_TP.output.trimed_fq_R1,
        fq_R2 = rules.RNA_cutadapt_TP.output.trimed_fq_R2,
    output:
        unmap_R1 = os.path.join(RNA_TP_BASE, "rm_rRNA", "WIV04.{tp}.{part}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_TP_BASE, "rm_rRNA", "WIV04.{tp}.{part}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_TP_BASE, "rm_rRNA", "WIV04.{tp}.{part}"),
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

rule RNA_STAR_TP:
    input:
        fq_R1 = rules.RNA_rm_rRNA_TP.output.unmap_R1,
        fq_R2 = rules.RNA_rm_rRNA_TP.output.unmap_R2,
        gtf = hg38_REF_GTF,
    output:
        bam = os.path.join(RNA_TP_BASE, "STAR", "WIV04.{tp}.{part}.sort.bam"),
        unmap_R1 = os.path.join(RNA_TP_BASE, "STAR", "WIV04.{tp}.{part}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_TP_BASE, "STAR", "WIV04.{tp}.{part}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_TP_BASE, "STAR", "WIV04.{tp}.{part}"),
        indx = hg38_INDX,
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


rule RNA_STAR_WIV04_TP:
    input:
        fq_R1 = rules.RNA_STAR_TP.output.unmap_R1,
        fq_R2 = rules.RNA_STAR_TP.output.unmap_R2,
    output:
        bam = os.path.join(RNA_TP_BASE, "STAR_WIV04_rm_host", "WIV04.{tp}.{part}.sort.bam"),
    params:
        mapping_dir = os.path.join(RNA_TP_BASE, "STAR_WIV04_rm_host", "WIV04.{tp}.{part}"),
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

rule merge_WIV04_bam_TP:
    input:
        bams = lambda wildcards: [rules.RNA_STAR_WIV04_TP.output.bam.format(part=part, tp=wildcards.tp) for part in RNA_TP_SAMPLE[wildcards.tp]]
    output:
        merged = os.path.join(RNA_TP_BASE, "STAR_WIV04_merge", "WIV04.{tp}.sort.bam"),
    shell:
        """
samtools merge -@ 28 -c -p {output.merged}.tmp.bam {input.bams}
samtools sort -@ 28 -o {output.merged} {output.merged}.tmp.bam
samtools index {output.merged}
rm {output.merged}.tmp.bam
        """

rule WIV04_bam2bw_TP:
    input:
        bam = rules.merge_WIV04_bam_TP.output.merged,
        size = WIV04_CHROM_SIZE,
    output:
        flag = touch(os.path.join(RNA_TP_BASE, "WIV04_bw", "{tp}.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(RNA_TP_BASE, "WIV04_bw"),
        rule = "1+-,1-+,2++,2--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/{wildcards.tp} -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/{wildcards.tp}.Forward.wig {input.size} {params.work_dir}/{wildcards.tp}.Forward.bw
wigToBigWig -clip {params.work_dir}/{wildcards.tp}.Reverse.wig {input.size} {params.work_dir}/{wildcards.tp}.Reverse.bw
rm {params.work_dir}/{wildcards.tp}.Forward.wig {params.work_dir}/{wildcards.tp}.Reverse.wig
        """

rule cnt_JS_TP:
    input:
        bam = rules.merge_WIV04_bam_TP.output.merged,
        genome = WIV04_GENOME_FA,
    output:
        js = os.path.join(RNA_TP_BASE, "STAR_WIV04_JS", "{tp}_js.bed"),
        jd = os.path.join(RNA_TP_BASE, "STAR_WIV04_JS", "{tp}_jd.csv"),
    params:
        rnatk_ja = "python scripts/subgenome/subgenome.py",
        outprefix = os.path.join(RNA_TP_BASE, "STAR_WIV04_JS", "{tp}")
    shell:
        """
source activate py35
{params.rnatk_ja} --window 5 --mingap 30 --minevi 1 --armsize 20 --genome {input.genome} -o {params.outprefix} -v {input.bam}
        """

rule jspseq:
    input:
        js = rules.cnt_JS_TP.output.jd,
        genome = WIV04_GENOME_FA,
    output:
        jscsv = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.flank.csv"),
    threads: 1
    shell:
        """
        source activate rnatk
        rnatk_jsflank --mingap 30 --numevi 1 --strand '+' -f 'jd' -d {wildcards.direction} --length 20 {input.js} {input.genome} {output.jscsv}
        """

rule jspair2fa:
    input:
        jscsv = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.flank.csv"),
    output:
        jsfa = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}_pair.fa"),
    threads: 1
    shell:
        """
        source activate rnatk
        python scripts/sgcov/workflow/pairing/scripts/jdcsv2fa.py {input.jscsv} > {output.jsfa}
        """

rule rnahybrid:
    input:
        jscsv = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.flank.csv"),
    output:
        hyb = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.hyb.csv"),
    threads: 1
    shell:
        """
        python scripts/sgcov/workflow/pairing/scripts/rnahyb.py {input.jscsv} {output.hyb}
        """

rule hyb2cofold:
    input:
        hyb = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.hyb.csv"),
        fa = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}_pair.fa"),
    output:
        cofold = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.cofold.csv"),
    threads: 1
    shell:
        """
        source activate rnatk
        rnatk_hyb2cofold {input.hyb} {input.fa} {output.cofold}
        """

rule fixendpair:
    input:
        cofold = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}.cofold.csv"),
    output:
        cofold = os.path.join(RNA_TP_BASE, "STAR_WIV04_pairing", "{tp}.{direction}_hybcons.csv"),
    threads: 1
    shell:
        """
        source activate rnatk
        rnatk_fixcofold -d {wildcards.direction} {input.cofold} {output.cofold}
        """

rule RNA_seq_pair:
    input:
        expand(rules.fixendpair.output, tp=RNA_TP_SAMPLE.keys(), direction=["lft", "rgt"])

rule RNA_seq_TP:
    input:
        expand(rules.cnt_JS_TP.output, tp=RNA_TP_SAMPLE.keys()),
        expand(rules.WIV04_bam2bw_TP.output, tp=RNA_TP_SAMPLE.keys()),
