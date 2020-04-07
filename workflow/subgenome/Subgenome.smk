import os
DATA_DIR = os.path.join("data", "subgenome")
RESULT_DIR = os.path.join("analysis", "subgenome")

## Figure 1, S2 & S3
rule filter_jncsite:
    input:
        bw = expand(os.path.join(DATA_DIR, "bw", "{source}.Forward.bw"), source=["NGS", "Nanopore"]),
        js = expand(os.path.join(DATA_DIR, "js", "CoV_{rep}_jd.csv"), rep=["Rep1", "Rep2", "merged"]),
    output:
        consistent_gap = os.path.join(RESULT_DIR, "NGS_Nanopore", "all_consistent_gap.tsv"),
    shell:
        "Rscript scripts/filter_jncsite.R"

## Figure 2 & S4
rule build_subgenome:
    input:
        rules.filter_jncsite.output,
    output:
        subgenome = os.path.join(RESULT_DIR, "subgenome", "subgenome.group.tsv"),
    shell:
        "Rscript scripts/build_subgenome.R"

## Figure 3, 4 & S5
rule identify_rule:
    input:
        rules.filter_jncsite.output,
    output:
        pair_ht = os.path.join(RESULT_DIR, "seq_pair", "seq_pair.leader_consistent.by_gap.ht.pdf"),
        motif = os.path.join(RESULT_DIR, "seq_pair", "Motif.cnt.by_gap.pdf"),
    shell:
        "Rscript scripts/identify_rule.R"

## Figure 5
rule nanopore_long_JS:
    input:
        nanopore_bed = expand(os.path.join(DATA_DIR, "nanopore", "data", "WIV04.{rep}.gap.bed"), rep=["rep1", "rep2"]),
    output:
        ORF1ab_bed = os.path.join(RESULT_DIR, "long_JS", "ORF1ab.all.bed")
    shell:
        "Rscript scripts/nanopore_long_JS.R"

rule subgenome_other:
    input:
        rules.identify_rule.output,
    output:
        polyA = os.path.join(RESULT_DIR, "Nanopore", "polyA.pdf"),
        pair_curve = os.path.join(RESULT_DIR, "seq_pair", "JS.gap_type.curve.consistent.pdf")
    shell:
        "Rscript scripts/subgenome_other.R"
