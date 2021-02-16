import os
PXD_ID = ["PXD018241", "PXD018594", "PXD021120"]
rule Jss_pepToBed:
    input:
        pep = os.path.join("data", "JSs", "MS_results", "{PXD_id}", "peptides.txt"),
        fa = os.path.join("data", "annotation", "WIV04.fasta"),
    output:
        bed = os.path.join("results", "JSs", "analysis", "{PXD_id}.pepToBed.bed"),
    shell:
        """
        python scripts/Jss.pepToBed.py -i {input.pep} -g {input.fa} -o {output.bed}
        """

rule all:
    input:
        expand(rules.Jss_pepToBed.output.bed, PXD_id = PXD_ID),
