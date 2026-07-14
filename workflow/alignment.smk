rule align_bowtie2:
    input:
        trimmed1="results/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed2="results/trimmed/{sample}_R2_val_2.fq.gz",
        index=bowtie2_index_files
    output:
        temp("results/bam/raw/{sample}.bam")
    threads: 8
    params:
        index=REF["bt2_index"],
        min_mapq=MIN_MAPQ,
        max_insert_size=MAX_INSERT_SIZE
    log:
        "results/logs/bowtie2_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/bam/raw results/logs
        bowtie2 -x {params.index:q} -1 {input.trimmed1:q} -2 {input.trimmed2:q} \
            -X {params.max_insert_size} \
            --no-mixed --no-discordant -p {threads} 2> {log:q} \
          | samtools view -b -q {params.min_mapq} -F 1804 -f 2 \
                -o {output:q} - 2>> {log:q}
        """


rule sort_markdup:
    input:
        bam="results/bam/raw/{sample}.bam"
    output:
        bam="results/bam/{sample}.sorted.bam",
        bai="results/bam/{sample}.sorted.bam.bai",
        flagstat="results/qc/samtools/{sample}.flagstat.txt",
        markdup_stats="results/qc/samtools/{sample}.markdup.txt"
    threads: 4
    log:
        "results/logs/sort_markdup_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/bam results/logs results/qc/samtools
        samtools sort -n -@ {threads} -o results/bam/{wildcards.sample}.name_sorted.bam {input.bam:q} > {log:q} 2>&1
        samtools fixmate -m results/bam/{wildcards.sample}.name_sorted.bam results/bam/{wildcards.sample}.fixmate.bam >> {log:q} 2>&1
        samtools sort -@ {threads} -o results/bam/{wildcards.sample}.coord_sorted.bam results/bam/{wildcards.sample}.fixmate.bam >> {log:q} 2>&1
        samtools markdup -r -s results/bam/{wildcards.sample}.coord_sorted.bam {output.bam:q} 2> {output.markdup_stats:q}
        cat {output.markdup_stats:q} >> {log:q}
        samtools index -@ {threads} {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        samtools flagstat -@ {threads} {output.bam:q} > {output.flagstat:q} 2>> {log:q}
        rm -f results/bam/{wildcards.sample}.name_sorted.bam results/bam/{wildcards.sample}.fixmate.bam results/bam/{wildcards.sample}.coord_sorted.bam
        """
