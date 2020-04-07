import os
MERGE_FA = os.path.join("data", "ChlSab_WIV04_merge.fa")
DATA_BASE = "data"
RESULTS_BASE = "results"

rule minimap2_map:
    input:
        merge_fa = MERGE_FA,
        raw_fq = os.path.join(DATA_BASE, "{sample}.fq"),
    output:
        bam = temp(os.path.join(RESULTS_BASE, "minimap2_map", "{sample}.bam")),
        sort_bam = os.path.join(RESULTS_BASE, "minimap2_map", "{sample}.merge.sort.bam"),
        virus_bam = os.path.join(RESULTS_BASE, "minimap2_map", "{sample}.virus.sort.bam"),
        merge_flagstat = os.path.join(RESULTS_BASE, "minimap2_map", "{sample}.merge_flagstat.txt"),
        virus_flagstat = os.path.join(RESULTS_BASE, "minimap2_map", "{sample}.virus_flagstat.txt"),
    threads: 20,
    shell:
        """
        minimap2 -ax splice -un -k14 --no-end-flt --secondary=no \
                -t {threads} {input.merge_fa} {input.raw_fq} | samtools view -@ {threads} -bh > {output.bam}
        samtools sort -@ {threads} {output.bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        samtools view -bh {output.sort_bam} MN996528 > {output.virus_bam}
        samtools index {output.virus_bam}
        samtools flagstat {output.sort_bam} > {output.merge_flagstat}
        samtools flagstat {output.virus_bam} > {output.virus_flagstat}
        """
