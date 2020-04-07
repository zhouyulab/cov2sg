import os
include: "MetaInfo.smk"

DATA = os.path.join("data")
ROOT = os.path.join("results")
GENOME = os.path.join("genome")

SAMPLEs = ['rand', ]

rule all:
    input:
        expand("results/ja/{sample}_jd.csv", sample=SAMPLEs),


rule ja_rand:
    input:
    output:
        jd = os.path.join(ROOT, "ja", "{sample}_jd.csv"),
    params:
        genomesize = 29891 - 21, # remove #polyA in tail
        outprefix = os.path.join(ROOT, "ja", "{sample}"),
    threads: 1
    shell:
        """
        rnatk_jsrand {params.genomesize} MN996528 {output.jd} -n 10000 --mingap 30 --armsize 20
        """

