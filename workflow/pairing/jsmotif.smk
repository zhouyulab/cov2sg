import os
include: "MetaInfo.smk"

DATA = os.path.join("data")
ROOT = os.path.join("results")
GENOME = os.path.join("genome")

SAMPLEs = ["CoV_merged", 'rand', ]

rule all:
    input:
        "results/jsmotif/TRScore.bed",
        expand(os.path.join(ROOT, "jsmotif", "{sample}_20.bed"), sample=SAMPLEs),
        expand(os.path.join(ROOT, "jsmotif", "{sample}_cnt.csv"), sample=SAMPLEs),

rule ja_rand:
    input:
        genome = os.path.join(GENOME, "WIV04.fasta"),
    output:
        motif = "results/jsmotif/TRScore.bed",
    threads: 1
    shell:
        """
        python scripts/find_cmotif.py {input.genome} {output.motif}
        """

rule jspair2fa:
    input:
        jdcsv = os.path.join(ROOT, "ja", "{sample}_jd.csv"),
    output:
        bed = os.path.join(ROOT, "jsmotif", "{sample}_20.bed"),
    threads: 1
    shell:
        """
        python scripts/js2bed.py {input.jdcsv} 20 {output.bed}
        """

rule cofold:
    input:
        motif = "results/jsmotif/TRScore.bed",
        bed = os.path.join(ROOT, "jsmotif", "{sample}_20.bed")
    output:
        cnt = os.path.join(ROOT, "jsmotif", "{sample}_cnt.csv"),
    threads: 1
    shell:
        """
        #grep -e '+' {input.motif} > {input.motif}.fwd.bed
        bedtools intersect -a {input.bed} -b {input.motif} -s -F 1.0 -c > {output.cnt}
        """
