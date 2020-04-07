import os
include: "MetaInfo.smk"

DATA = os.path.join("data")
ROOT = os.path.join("results")
GENOME = os.path.join("genome")

rule all:
    input:
        expand("results/ja/{sample}_js.bed", sample=SAMPLEs),


rule ja:
    input:
        fbam = os.path.join(DATA, "NGS", "{sample}.sort.bam"),
        fbambai = os.path.join(DATA, "NGS", "{sample}.sort.bam.bai"),
        genome = os.path.join(GENOME, "WIV04.fasta"),
    output:
        js = os.path.join(ROOT, "ja", "{sample}_js.bed"),
        jg = os.path.join(ROOT, "ja", "{sample}_jg.bed")
    params:
        outprefix = os.path.join(ROOT, "ja", "{sample}"),
    threads: 1
    shell:
        """
        rnatk_ja --window 5 --mingap 30 --minevi 1 --armsize 20 --genome {input.genome} -o {params.outprefix} -v {input.fbam}
        """

