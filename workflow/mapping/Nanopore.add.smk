import os
include: "Nanopore.smk"

SUBGENOME_BASE = os.path.join("analysis", "subgenome")
SUBGENOME_OTHER_DATA = os.path.join(SUBGENOME_BASE, "other_data")

rule mapping_other_nanopore_fastq:
    input:
        fastq = os.path.join("data", "NanoporeFastq", "{project}", "{sample}.fastq"),
        merge_fa = os.path.join("data", "genome", "chkSab_WIV04.merge.fa"),
    output:
        bam = temp(os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "minimap2_map", "{project}", "{sample}.bam")),
        sort_bam = os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "minimap2_map", "{project}", "{sample}.merge.sort.bam"),
        virus_bam = os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "minimap2_map", "{project}", "{sample}.virus.sort.bam"),
        merge_flagstat = os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "minimap2_map", "{project}", "{sample}.merge_flagstat.txt"),
        virus_flagstat = os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "minimap2_map", "{project}", "{sample}.virus_flagstat.txt"),
    threads: 20,
    shell:
        """
        minimap2 -ax splice -un -k14 --no-end-flt --secondary=no \
                -t {threads} {input.merge_fa} {input.fastq} | samtools view -@ {threads} -bh > {output.bam}
        samtools sort -@ {threads} {output.bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        samtools view -bh {output.sort_bam} MN996528 > {output.virus_bam}
        samtools index {output.virus_bam}
        samtools flagstat {output.sort_bam} > {output.merge_flagstat}
        samtools flagstat {output.virus_bam} > {output.virus_flagstat}
        """

rule other_nanopore_filter_host_bam:
    input:
        bam = rules.mapping_other_nanopore_fastq.output.sort_bam,
    output:
        host_bam = os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "host", "bam", "{project}", "{sample}.host.bam"),
    params:
        filter_host_bam = "python scripts/subgenome/filter_host_bam.py",
    shell:
        """
source activate py35
{params.filter_host_bam} -i {input.bam} -o {output.host_bam} --virus-chrom MN996528
samtools index {output.host_bam}
        """

rule other_nanopore_host_bam2bed:
    input:
        bam = rules.other_nanopore_filter_host_bam.output.host_bam,
    output:
        bed =  os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "host", "bed", "{project}", "{sample}.host.bed"),
    shell:
        "bedtools bamtobed -bed12 -i {input.bam} | sort -k 1,1 -k 2,2n > {output.bed}"

rule other_nanopore_bam2bed:
    input:
        bam = rules.mapping_other_nanopore_fastq.output.virus_bam,
    output: 
        bed =  os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "data", "{project}", "WIV04.{sample}.bed"),
    shell:
        "bedtools bamtobed -bed12 -i {input.bam} | sort -k 1,1 -k 2,2n > {output.bed}"

rule other_nanopore_gap:
    input:
        bed = rules.other_nanopore_bam2bed.output.bed,
    output: 
        gap = os.path.join(SUBGENOME_OTHER_DATA, "nanopore", "data", "{project}", "WIV04.{sample}.gap.bed"),
    shell:
        """
source activate py35
iso2intron -i {input.bed} -o {output.gap}
        """

rule WIV04_other_nanopore_bam2bw:
    input:
        bam = rules.mapping_other_nanopore_fastq.output.virus_bam,
        size = WIV04_CHROM_SIZE,
    output:
        flag = touch(os.path.join(SUBGENOME_OTHER_DATA, "bw", "{project}", "Nanopore.{sample}.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(SUBGENOME_OTHER_DATA, "bw", "{project}"),
        rule = "++,--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/Nanopore.{wildcards.sample} -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/Nanopore.{wildcards.sample}.Forward.wig {input.size} {params.work_dir}/Nanopore.{wildcards.sample}.Forward.bw
wigToBigWig -clip {params.work_dir}/Nanopore.{wildcards.sample}.Reverse.wig {input.size} {params.work_dir}/Nanopore.{wildcards.sample}.Reverse.bw
rm {params.work_dir}/Nanopore.{wildcards.sample}.Forward.wig {params.work_dir}/Nanopore.{wildcards.sample}.Reverse.wig
        """

rule subgenome_other:
    input:
        expand(rules.other_nanopore_gap.output, zip, project=["NarryKim", "TP", "TP", "TP", "TP"], sample=["VeroInf24h", "Vero4h", "Vero6h", "Vero12h", "Vero24h"]),
        expand(rules.WIV04_other_nanopore_bam2bw.output, zip, project=["NarryKim", "TP", "TP", "TP", "TP"], sample=["VeroInf24h", "Vero4h", "Vero6h", "Vero12h", "Vero24h"]),
        expand(rules.other_nanopore_host_bam2bed.output, zip, project=["NarryKim", "TP", "TP", "TP", "TP"], sample=["VeroInf24h", "Vero4h", "Vero6h", "Vero12h", "Vero24h"]),


SUBGENOME_CACO2_DATA = os.path.join(SUBGENOME_BASE, "Caco2")
CACO2_NANOPORE_TP = ["4h", "12h", "24h"]
CACO2_NANOPORE_REP = ["rep1"]
rule mapping_caco2_nanopore_fastq:
    input:
        fastq = os.path.join("data", "Nanopore_caco2", "Caco2.{tp}.{rep}.fastq.gz"),
        merge_fa = os.path.join("data", "genome", "hg38_WIV04.merge.fa"),
    output:
        bam = temp(os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "minimap2_map", "Caco2.{tp}.{rep}.bam")),
        sort_bam = os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "minimap2_map", "Caco2.{tp}.{rep}.merge.sort.bam"),
        virus_bam = os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "minimap2_map", "Caco2.{tp}.{rep}.virus.sort.bam"),
        merge_flagstat = os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "minimap2_map", "Caco2.{tp}.{rep}.merge_flagstat.txt"),
        virus_flagstat = os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "minimap2_map", "Caco2.{tp}.{rep}.virus_flagstat.txt"),

    threads: 20,
    shell:
        """
        minimap2 -ax splice -un -k14 --no-end-flt --secondary=no \
                -t {threads} {input.merge_fa} {input.fastq} | samtools view -@ {threads} -bh > {output.bam}
        samtools sort -@ {threads} {output.bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        samtools view -bh {output.sort_bam} MN996528 > {output.virus_bam}
        samtools index {output.virus_bam}
        samtools flagstat {output.sort_bam} > {output.merge_flagstat}
        samtools flagstat {output.virus_bam} > {output.virus_flagstat}
        """

rule caco2_filter_host_bam:
    input:
        bam = rules.mapping_caco2_nanopore_fastq.output.sort_bam,
    output:
        host_bam = os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "host", "bam", "Caco2.{tp}.{rep}.host.sort.bam"),
    params:
        filter_host_bam = "python scripts/subgenome/filter_host_bam.py",
    shell:
        """
source activate py35
{params.filter_host_bam} -i {input.bam} -o {output.host_bam} --virus-chrom MN996528
samtools index {output.host_bam}
        """
        
rule caco2_nanopore_host_bam2bed:
    input:
        bam = rules.caco2_filter_host_bam.output.host_bam,
    output:
        bed =  os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "host", "bed", "Caco2.{tp}.{rep}.host.bed"),
    shell:
        "bedtools bamtobed -bed12 -i {input.bam} | sort -k 1,1 -k 2,2n > {output.bed}"


rule caco2_nanopore_bam2bed:
    input:
        bam = rules.mapping_caco2_nanopore_fastq.output.virus_bam,
    output:
        bed =  os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "data", "Caco2.{tp}.{rep}.bed"),
    shell:
        "bedtools bamtobed -bed12 -i {input.bam} | sort -k 1,1 -k 2,2n > {output.bed}"

rule caco2_nanopore_gap:
    input:
        bed = rules.caco2_nanopore_bam2bed.output.bed,
    output:
        gap = os.path.join(SUBGENOME_CACO2_DATA, "nanopore", "data", "Caco2.{tp}.{rep}.gap.bed"),
    shell:
        """
source activate py35
iso2intron -i {input.bed} -o {output.gap}
        """

rule caco2_nanopore_bam2bw:
    input:
        bam = rules.mapping_caco2_nanopore_fastq.output.virus_bam,
        size = os.path.join("data", "genome", "WIV04.sizes"),
    output:
        flag = touch(os.path.join(SUBGENOME_CACO2_DATA, "Nanopore_bw", "Caco2.{tp}.{rep}.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(SUBGENOME_CACO2_DATA, "Nanopore_bw"),
        rule = "++,--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep} -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep}.Forward.wig {input.size} {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep}.Forward.bw
#wigToBigWig -clip {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep}.Reverse.wig {input.size} {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep}.Reverse.bw
rm {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep}.Forward.wig {params.work_dir}/Caco2.{wildcards.tp}.{wildcards.rep}.Reverse.wig
        """

rule caco2_nanopore:
    input:
        expand(rules.caco2_nanopore_gap.output, tp=CACO2_NANOPORE_TP, rep=CACO2_NANOPORE_REP),
        expand(rules.caco2_nanopore_bam2bw.output, tp=CACO2_NANOPORE_TP, rep=CACO2_NANOPORE_REP),
        expand(rules.caco2_nanopore_host_bam2bed.output, tp=CACO2_NANOPORE_TP, rep=CACO2_NANOPORE_REP),

SUBGENOME_HCoV_DATA = os.path.join(SUBGENOME_BASE, "HCoV")

rule mapping_HCoV_nanopore_fastq:
    input:
        fastq = os.path.join("data", "HCoV_Nanopore", "HCoV.{sample}.fastq.gz"),
        fa = os.path.join("data", "genome", "HCoV.{sample}.fa"),
    output:
        bam = os.path.join(SUBGENOME_HCoV_DATA, "nanopore", "minimap2_map", "{sample}.bam"),
        flagstat = os.path.join(SUBGENOME_HCoV_DATA, "nanopore", "minimap2_map", "{sample}.flagstat.txt"),
    threads: 20,
    shell:
        """
        minimap2 -ax splice -un -k14 --no-end-flt --secondary=no \
                -t {threads} {input.fa} {input.fastq} | samtools view -@ {threads} -bh > {output.bam}.tmp.bam
        samtools sort -@ {threads} {output.bam}.tmp.bam > {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        rm {output.bam}.tmp.bam
        """

rule HCoV_nanopore_bam2bed:
    input:
        bam = rules.mapping_HCoV_nanopore_fastq.output.bam,
    output: 
        bed =  os.path.join(SUBGENOME_HCoV_DATA, "nanopore", "data", "HCoV.{sample}.bed"),
    shell:
        "bedtools bamtobed -bed12 -i {input.bam} | sort -k 1,1 -k 2,2n > {output.bed}"

rule HCoV_nanopore_gap:
    input:
        bed = rules.HCoV_nanopore_bam2bed.output.bed,
    output: 
        gap = os.path.join(SUBGENOME_HCoV_DATA, "nanopore", "data", "HCoV.{sample}.gap.bed"),
    shell:
        """
source activate py35
iso2intron -i {input.bed} -o {output.gap}
        """

rule HCoV_nanopore_RNAhybrid:
    input:
        intron = rules.HCoV_nanopore_gap.output.gap,
        genome = os.path.join("data", "genome", "HCoV.{sample}.fa"),
    output:
        pair = os.path.join(SUBGENOME_HCoV_DATA, "pairing", "intron.{sample}.tsv"),
    params:
        intron = os.path.join(SUBGENOME_HCoV_DATA, "pairing", "intron.{sample}.bed"),
        RNAhybrid_js = "python scripts/subgenome/RNAhybrid_js.py",
    shell:
        """
source activate py35
{params.RNAhybrid_js} -i {input.intron} -g {input.genome} -o {output.pair}
        """

rule HCoV_nanopore_RNAhybrid_tmp:
    input:
        expand(rules.HCoV_nanopore_RNAhybrid.output, sample=["WT", "SL2"])

rule HCoV_nanopore_bam2bw:
    input:
        bam = rules.mapping_HCoV_nanopore_fastq.output.bam,
        size = os.path.join("data", "genome", "HCoV.{sample}.sizes"),
    output:
        flag = touch(os.path.join(SUBGENOME_HCoV_DATA, "bw", "Nanopore.{sample}.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(SUBGENOME_HCoV_DATA, "bw"),
        rule = "++,--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/Nanopore.{wildcards.sample} -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/Nanopore.{wildcards.sample}.Forward.wig {input.size} {params.work_dir}/Nanopore.{wildcards.sample}.Forward.bw
#wigToBigWig -clip {params.work_dir}/Nanopore.{wildcards.sample}.Reverse.wig {input.size} {params.work_dir}/Nanopore.{wildcards.sample}.Reverse.bw
rm {params.work_dir}/Nanopore.{wildcards.sample}.Forward.wig {params.work_dir}/Nanopore.{wildcards.sample}.Reverse.wig
        """


rule HCoV:
    input:
        expand(rules.HCoV_nanopore_bam2bw.output, sample=["WT", "SL2"]),
        expand(rules.HCoV_nanopore_gap.output, sample=["WT", "SL2"]),

