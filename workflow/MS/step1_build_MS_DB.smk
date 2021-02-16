import os

rule Jss_make_MS_db:
    input:
        fa_ORF1ab = os.path.join("data", "JSs", "WIV04.nanopore.ORF1ab.iso.fa"),
        fa_S2N = os.path.join("data", "JSs", "WIV04.nanopore.S2N.iso.fa"),
    output:
        fa = os.path.join("results", "JSs", "Jss.make_MS_db.fasta"),
    shell:
        """
        python scripts/Jss.make_MS_db.py -p data/JSs/ -o {output.fa}
        """

rule Jss_make_MS_db_bed:
    input:
        fa = rules.Jss_make_MS_db.output.fa,
    output:
        bed = os.path.join("results", "JSs", "Jss.make_MS_db.bed"),
    shell:
        """
        python scripts/Jss.make_MS_db_bed.py -i {input.fa} -p data/JSs/ -o {output.bed}
        """


rule Jss_genome_translation_6f:
    input:
        fa_genome = os.path.join("data", "annotation", "WIV04.fasta"),
    output:
        fa = os.path.join("results", "JSs", "Jss.genome_translation_6f.fasta"),
        bed = os.path.join("results", "JSs", "Jss.genome_translation_6f.bed"),
    shell:
        """
        python scripts/Jss.genome_translation_6f.py -i {input.fa_genome} -f {output.fa} -b {output.bed}
        """

rule DB_merge:
    input:
        fa_1 = rules.Jss_make_MS_db.output.fa,
        fa_2 = rules.Jss_genome_translation_6f.output.fa,
        fa_3 = os.path.join("data", "JSs", "uniprot-human-20200928-192656.fasta"),
        fa_4 = os.path.join("data", "JSs", "uniprot-greenMonkey-20201022-19525.fasta"),
    output:
        fa_human = os.path.join("results", "JSs", "Jss.human_merge.fasta"),
        fa_monkey = os.path.join("results", "JSs", "Jss.monkey_merge.fasta"),
    shell:
        """
        cat {input.fa_1} > {output.fa_human}
        cat {input.fa_2} >> {output.fa_human}
        cat {input.fa_3} >> {output.fa_human}
        cat {input.fa_1} > {output.fa_monkey}
        cat {input.fa_2} >> {output.fa_monkey}
        cat {input.fa_4} >> {output.fa_monkey}
        """

rule all:
    input:
        rules.Jss_make_MS_db.output.fa,
        rules.Jss_make_MS_db_bed.output.bed,
        rules.Jss_genome_translation_6f.output.fa,
        rules.DB_merge.output.fa_human,
        expand(rules.Jss_pepToBed.output.bed, PXD_id = PXD_ID),
