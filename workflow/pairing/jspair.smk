import os
include: "MetaInfo.smk"

DATA = os.path.join("data")
ROOT = os.path.join("results")
GENOME = os.path.join("genome")
DIRECTIONs = ('lft', 'rgt')

SAMPLEs = SAMPLEs + ['rand']

rule all:
    input:
        expand(os.path.join(ROOT, "jspair", "{sample}_{direction}.csv"), sample=SAMPLEs, direction=DIRECTIONs),
        expand(os.path.join(ROOT, "jspair", "{sample}_{direction}_pair.fa"), sample=SAMPLEs, direction=DIRECTIONs),


rule jspseq:
    input:
        js = os.path.join(ROOT, "ja", "{sample}_jd.csv"),
        genome = os.path.join(GENOME, "WIV04.fasta"),
    output:
        jscsv = os.path.join(ROOT, "jspair", "{sample}_{direction}.csv"),
    threads: 1
    shell:
        """
        rnatk_jsflank --mingap 30 --numevi 1 --strand '+' -f 'jd' -d {wildcards.direction} --length 20 {input.js} {input.genome} {output.jscsv}
        """

rule jspair2fa:
    input:
        jscsv = os.path.join(ROOT, "jspair", "{sample}_{direction}.csv"),
    output:
        jsfa = os.path.join(ROOT, "jspair", "{sample}_{direction}_pair.fa"),
    threads: 1
    shell:
        """
        python scripts/jdcsv2fa.py {input.jscsv} > {output.jsfa}
        """

