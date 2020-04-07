import os
include: "MetaInfo.smk"

DATA = os.path.join("data")
ROOT = os.path.join("results")
GENOME = os.path.join("genome")

DIRECTIONs = ('lft', 'rgt')
SAMPLEs = ["CoV_merged", 'rand']

rule all:
    input:
        expand(os.path.join(ROOT, "jshyb", "{sample}_{direction}.csv"), sample=SAMPLEs, direction=DIRECTIONs),
        expand(os.path.join(ROOT, "jshyb", "cofold_{sample}_{direction}.csv"), sample=SAMPLEs, direction=DIRECTIONs),
        expand(os.path.join(ROOT, "jshyb", "{sample}_{direction}_hybcons.csv"), sample=SAMPLEs, direction=DIRECTIONs),


rule rnahybrid:
    input:
        jscsv = os.path.join(ROOT, "jspair", "{sample}_{direction}.csv"),
    output:
        hyb = os.path.join(ROOT, "jshyb", "{sample}_{direction}.csv"),
    threads: 1
    shell:
        """
        python scripts/rnahyb.py {input.jscsv} {output.hyb}
        """

rule hyb2cofold:
    input:
        hyb = os.path.join(ROOT, "jshyb", "{sample}_{direction}.csv"),
        fa = os.path.join(ROOT, "jspair", "{sample}_{direction}_pair.fa"),
    output:
        cofold = os.path.join(ROOT, "jshyb", "cofold_{sample}_{direction}.csv"),
    threads: 1
    shell:
        """
        rnatk_hyb2cofold {input.hyb} {input.fa} {output.cofold}
        """


rule fixendpair:
    input:
        cofold = os.path.join(ROOT, "jshyb", "cofold_{sample}_{direction}.csv"),
    output:
        cofold = os.path.join(ROOT, "jshyb", "{sample}_{direction}_hybcons.csv"),
    threads: 1
    shell:
        """
        rnatk_fixcofold -d {wildcards.direction} {input.cofold} {output.cofold}
        """
