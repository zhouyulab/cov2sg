import os

SUBGENOME_BASE = os.path.join("analysis", "subgenome")

rule WIV04_nanopore_bam2bed:
    input:
        bam = os.path.join("data", "Nanopore", "WIV04", "CoV.{rep}.sort.bam"),
    output:
        bed =  os.path.join(SUBGENOME_BASE, "nanopore", "data", "WIV04.{rep}.bed"),
    shell:
        "bedtools bamtobed -bed12 -i {input.bam} | sort -k 1,1 -k 2,2n > {output.bed}"

rule WIV04_nanopore_gap:
    input:
        bed = rules.WIV04_nanopore_bam2bed.output.bed,
    output:
        gap = os.path.join(SUBGENOME_BASE, "nanopore", "data", "WIV04.{rep}.gap.bed"),
    shell:
        """
source activate py35
iso2intron -i {input.bed} -o {output.gap}
        """

rule WIV04_NGS_bam2bw:
    input:
        bam = rules.build_WIV04_NGS_pool.output.merged,
        size = WIV04_CHROM_SIZE,
    output:
        flag = touch(os.path.join(SUBGENOME_BASE, "bw", "NGS.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(SUBGENOME_BASE, "bw"),
        rule = "1+-,1-+,2++,2--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/NGS -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/NGS.Forward.wig {input.size} {params.work_dir}/NGS.Forward.bw
wigToBigWig -clip {params.work_dir}/NGS.Reverse.wig {input.size} {params.work_dir}/NGS.Reverse.bw
rm {params.work_dir}/NGS.Forward.wig {params.work_dir}/NGS.Reverse.wig
        """

rule merge_WIV04_nanopore_bam:
    input:
        bams = expand(os.path.join("data", "Nanopore", "WIV04", "CoV.{rep}.sort.bam"), rep=["rep1", "rep2"])
    output:
        merged = os.path.join(SUBGENOME_BASE, "WIV04.merged.bam"),
    threads: 32
    shell:
        """
samtools merge -@ {threads} -c -p {output.merged}.tmp.bam {input.bams}
samtools sort -@ {threads} -o {output.merged} {output.merged}.tmp.bam
samtools index {output.merged}
        """

rule WIV04_nanopore_bam2bw:
    input:
        bam = rules.merge_WIV04_nanopore_bam.output.merged,
        size = WIV04_CHROM_SIZE,
    output:
        flag = touch(os.path.join(SUBGENOME_BASE, "bw", "Nanopore.flag")),
    params:
        total_sum = 1000000,
        work_dir = os.path.join(SUBGENOME_BASE, "bw"),
        rule = "++,--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/Nanopore -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/Nanopore.Forward.wig {input.size} {params.work_dir}/Nanopore.Forward.bw
wigToBigWig -clip {params.work_dir}/Nanopore.Reverse.wig {input.size} {params.work_dir}/Nanopore.Reverse.bw
rm {params.work_dir}/Nanopore.Forward.wig {params.work_dir}/Nanopore.Reverse.wig
        """

rule subgenome:
    input:
        rules.WIV04_NGS_bam2bw.output,
        rules.WIV04_nanopore_bam2bw.output,
        expand(rules.WIV04_nanopore_gap.output, rep=["rep1", "rep2"]),

