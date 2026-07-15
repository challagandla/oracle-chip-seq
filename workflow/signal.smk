rule bamcoverage:
    input:
        bam="results/bam/{sample}.sorted.bam",
        bai="results/bam/{sample}.sorted.bam.bai",
        blacklist=lambda wc: REF["black_list"]
    output:
        bw="results/bigwig/{sample}.rpgc.bw"
    threads: DEEPTOOLS_THREADS
    params:
        effective_size=REF["effective_genome_size"],
        bin_size=TRACK_BIN_SIZE
    log:
        "results/logs/bamcoverage_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/bigwig results/logs
        bamCoverage -b {input.bam:q} -o {output.bw:q} \
            --normalizeUsing RPGC --effectiveGenomeSize {params.effective_size} \
            --binSize {params.bin_size} --extendReads --samFlagInclude 64 \
            --blackListFileName {input.blacklist:q} \
            --numberOfProcessors {threads} > {log:q} 2>&1
        """


rule bamcompare:
    input:
        chip="results/bam/{sample}.sorted.bam",
        chip_bai="results/bam/{sample}.sorted.bam.bai",
        control=lambda wc: f"results/bam/{CONTROL_MAP[wc.sample]}.sorted.bam",
        control_bai=lambda wc: f"results/bam/{CONTROL_MAP[wc.sample]}.sorted.bam.bai",
        blacklist=lambda wc: REF["black_list"]
    output:
        bw="results/bigwig/{sample}.log2ratio.bw"
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in CHIP_SAMPLES)
    threads: DEEPTOOLS_THREADS
    params:
        normalization=bamcompare_normalization_args,
        pseudocount=LOG2_PSEUDOCOUNT,
        bin_size=TRACK_BIN_SIZE
    log:
        "results/logs/bamcompare_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/bigwig results/logs
        bamCompare -b1 {input.chip:q} -b2 {input.control:q} -o {output.bw:q} \
            --operation log2 --pseudocount {params.pseudocount} {params.normalization} \
            --binSize {params.bin_size} --extendReads --samFlagInclude 64 \
            --blackListFileName {input.blacklist:q} \
            --numberOfProcessors {threads} > {log:q} 2>&1
        """
