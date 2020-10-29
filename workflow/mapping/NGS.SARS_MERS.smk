import os
include: "MetaInfo.smk"

RNA_OTHER_BASE = os.path.join("analysis", "RNA_seq_other")
RNA_OTHER_SAMPLE = {
    "SARS": ["rep1", "rep2"], 
    "MERS": ["rep1", "rep2"]
}


rule RNA_cutadapt_other:
    input:
        fastq_R1 = os.path.join("data", "other_data", "{species}", "{species}.{rep}.R1.fastq.gz"),
        fastq_R2 = os.path.join("data", "other_data", "{species}", "{species}.{rep}.R2.fastq.gz"),
    output:
        trimed_fq_R1 = os.path.join(RNA_OTHER_BASE, "trim", "{species}", "{species}.{rep}.trim.R1.fastq.gz"),
        trimed_fq_R2 = os.path.join(RNA_OTHER_BASE, "trim", "{species}", "{species}.{rep}.trim.R2.fastq.gz"),
    params:
        adaptor_R1 = "AGATCGGAAGAGC",
        adaptor_R2 = "AGATCGGAAGAGC",
        min_read_len = 20,
    shell:
        "cutadapt --cores 20 -e 0.25 --match-read-wildcards -a {params.adaptor_R1} -A {params.adaptor_R2} -m {params.min_read_len} -o {output.trimed_fq_R1} -p {output.trimed_fq_R2} {input.fastq_R1} {input.fastq_R2}; fastqc {output.trimed_fq_R1} {output.trimed_fq_R2} &"

rule RNA_rm_rRNA_other:
    input:
        fq_R1 = rules.RNA_cutadapt_other.output.trimed_fq_R1,
        fq_R2 = rules.RNA_cutadapt_other.output.trimed_fq_R2,
    output:
        unmap_R1 = os.path.join(RNA_OTHER_BASE, "rm_rRNA", "{species}", "{species}.{rep}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_OTHER_BASE, "rm_rRNA", "{species}", "{species}.{rep}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_OTHER_BASE, "rm_rRNA", "{species}", "{species}.{rep}"),
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

rule RNA_STAR_other:
    input:
        fq_R1 = rules.RNA_rm_rRNA_other.output.unmap_R1,
        fq_R2 = rules.RNA_rm_rRNA_other.output.unmap_R2,
        gtf = hg38_REF_GTF,
    output:
        bam = os.path.join(RNA_OTHER_BASE, "STAR", "{species}", "{species}.{rep}.sort.bam"),
        unmap_R1 = os.path.join(RNA_OTHER_BASE, "STAR", "{species}", "{species}.{rep}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_OTHER_BASE, "STAR", "{species}", "{species}.{rep}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_OTHER_BASE, "STAR", "{species}", "{species}.{rep}"),
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

rule RNA_STAR_WIV04_other:
    input:
        fq_R1 = rules.RNA_STAR_other.output.unmap_R1,
        fq_R2 = rules.RNA_STAR_other.output.unmap_R2,
    output:
        bam = os.path.join(RNA_OTHER_BASE, "rm_host", "{species}", "{species}.{rep}.sort.bam"),
    params:
        mapping_dir = os.path.join(RNA_OTHER_BASE, "rm_host", "{species}", "{species}.{rep}"),
        indx = "data/index/STAR/{species}",
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

rule WIV04_bam2bw_other:
    input:
        bam = rules.RNA_STAR_WIV04_other.output.bam,
        size = "data/genome/{species}.sizes",
    output:
        flag = touch(os.path.join(RNA_OTHER_BASE, "virus_bw", "{species}", "{species}.{rep}.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(RNA_OTHER_BASE, "virus_bw", "{species}"),
        rule = "1+-,1-+,2++,2--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/{wildcards.species}.{wildcards.rep} -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/{wildcards.species}.{wildcards.rep}.Forward.wig {input.size} {params.work_dir}/{wildcards.species}.{wildcards.rep}.Forward.bw
rm {params.work_dir}/{wildcards.species}.{wildcards.rep}.*.wig
        """

rule cnt_JS_other:
    input:
        bam = rules.RNA_STAR_WIV04_other.output.bam,
        genome = "data/genome/{species}.fasta",
    output:
        js = os.path.join(RNA_OTHER_BASE, "virus_JS", "{species}", "{species}.{rep}_js.bed"),
        jd = os.path.join(RNA_OTHER_BASE, "virus_JS", "{species}", "{species}.{rep}_jd.csv"),
    params:
        rnatk_ja = "python scripts/subgenome/subgenome.py",
        outprefix = os.path.join(RNA_OTHER_BASE, "virus_JS", "{species}", "{species}.{rep}")
    shell:
        """
source activate py35
{params.rnatk_ja} --window 5 --mingap 30 --minevi 1 --armsize 20 --genome {input.genome} -o {params.outprefix} -v {input.bam}
        """

rule jspseq:
    input:
        js = rules.cnt_JS_other.output.jd,
        genome = "data/genome/{species}.fasta",
    output:
        jscsv = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.flank.csv"),
    threads: 1
    shell:
        """
        source activate rnatk
        rnatk_jsflank --mingap 30 --numevi 1 --strand '+' -f 'jd' -d {wildcards.direction} --length 20 {input.js} {input.genome} {output.jscsv}
        """

rule jspair2fa:
    input:
        jscsv = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.flank.csv"),
    output:
        jsfa = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}_pair.fa"),
    threads: 1
    shell:
        """
        source activate rnatk
        python scripts/sgcov/workflow/pairing/scripts/jdcsv2fa.py {input.jscsv} > {output.jsfa}
        """

rule rnahybrid:
    input:
        jscsv = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.flank.csv"),
    output:
        hyb = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.hyb.csv"),
    threads: 1
    shell:
        """
        python scripts/sgcov/workflow/pairing/scripts/rnahyb.py {input.jscsv} {output.hyb}
        """

rule hyb2cofold:
    input:
        hyb = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.hyb.csv"),
        fa = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}_pair.fa"),
    output:
        cofold = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.cofold.csv"),
    threads: 1
    shell:
        """
        source activate rnatk
        rnatk_hyb2cofold {input.hyb} {input.fa} {output.cofold}
        """

rule fixendpair:
    input:
        cofold = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}.cofold.csv"),
    output:
        cofold = os.path.join(RNA_OTHER_BASE, "virus_pairing", "{species}", "{species}.{rep}.{direction}_hybcons.csv"),
    threads: 1
    shell:
        """
        source activate rnatk
        rnatk_fixcofold -d {wildcards.direction} {input.cofold} {output.cofold}
        """

rule RNA_seq_other:
    input:
        [rules.WIV04_bam2bw_other.output.flag.format(species=species, rep=rep) for species in RNA_OTHER_SAMPLE.keys() for rep in RNA_OTHER_SAMPLE[species]],
        [rules.fixendpair.output.cofold.format(species=species, rep=rep, direction=direction) for species in RNA_OTHER_SAMPLE.keys() for rep in RNA_OTHER_SAMPLE[species] for direction in ["lft", "rgt"]],
