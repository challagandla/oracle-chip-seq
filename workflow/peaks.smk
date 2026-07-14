rule call_peaks_broad:
    input:
        bam="results/bam/{sample}.sorted.bam",
        control=lambda wc: f"results/bam/{CONTROL_MAP[wc.sample]}.sorted.bam"
    output:
        broadpeak="results/peaks/raw/{sample}_peaks.broadPeak"
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in BROAD_SAMPLES) if BROAD_SAMPLES else "$^"
    params:
        gsize=REF["gsize"],
        broad_cutoff=BROAD_CUTOFF,
        outdir=lambda wc, output: os.path.dirname(str(output.broadpeak))
    log:
        "results/logs/macs3_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        macs3 callpeak -t {input.bam:q} -c {input.control:q} --format BAMPE \
            --name {wildcards.sample} --broad --broad-cutoff {params.broad_cutoff} \
            --keep-dup all --gsize {params.gsize:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """


rule call_peaks_narrow:
    input:
        bam="results/bam/{sample}.sorted.bam",
        control=lambda wc: f"results/bam/{CONTROL_MAP[wc.sample]}.sorted.bam"
    output:
        narrowpeak="results/peaks/raw/{sample}_peaks.narrowPeak"
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in NARROW_SAMPLES) if NARROW_SAMPLES else "$^"
    params:
        gsize=REF["gsize"],
        qvalue=NARROW_QVALUE,
        outdir=lambda wc, output: os.path.dirname(str(output.narrowpeak))
    log:
        "results/logs/macs3_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        macs3 callpeak -t {input.bam:q} -c {input.control:q} --format BAMPE \
            --name {wildcards.sample} --gsize {params.gsize:q} --keep-dup all \
            -q {params.qvalue} --outdir {params.outdir:q} > {log:q} 2>&1
        """


rule filter_blacklist:
    input:
        peaks="results/peaks/raw/{sample}_peaks.{ext}",
        blacklist=lambda wc: REF["black_list"],
        chrom_sizes=lambda wc: REF["chrom_sizes"]
    output:
        "results/peaks/{sample}_peaks.{ext}"
    log:
        "results/logs/filter_blacklist_{sample}_{ext}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/peaks results/logs
        bedtools intersect -v -a {input.peaks:q} -b {input.blacklist:q} 2> {log:q} \
          | bedtools sort -i - -g {input.chrom_sizes:q} > {output:q} 2>> {log:q}
        """


rule condition_consensus_peaks:
    input:
        peaks=factor_condition_peaks,
        chrom_sizes=lambda wc: REF["chrom_sizes"]
    output:
        "results/peaks/consensus/{factor}/{condition}.bed"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output[0])),
        min_replicates=CONSENSUS_MIN_REPLICATES
    log:
        "results/logs/consensus_{factor}_{condition}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        bedtools multiinter -i {input.peaks:q} 2> {log:q} \
          | awk -v min={params.min_replicates} 'BEGIN {{OFS="\t"}} $4 >= min {{print $1, $2, $3}}' \
          | bedtools sort -i - -g {input.chrom_sizes:q} 2>> {log:q} \
          | bedtools merge > {output:q} 2>> {log:q}
        if [[ ! -s {output:q} ]]; then
            echo "No replicate-supported peaks remained for {wildcards.factor}/{wildcards.condition}; inspect replicate peak calls" >> {log:q}
            exit 1
        fi
        """


rule factor_consensus_peaks:
    input:
        peaks=factor_condition_consensus,
        chrom_sizes=lambda wc: REF["chrom_sizes"]
    output:
        "results/peaks/consensus/{factor}.bed"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output[0]))
    log:
        "results/logs/consensus_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        cat {input.peaks:q} \
          | bedtools sort -i - -g {input.chrom_sizes:q} 2> {log:q} \
          | bedtools merge > {output:q} 2>> {log:q}
        if [[ ! -s {output:q} ]]; then
            echo "No factor-level consensus peaks remained for {wildcards.factor}" >> {log:q}
            exit 1
        fi
        """


rule merge_peaks:
    input:
        peaks=FACTOR_CONSENSUS,
        chrom_sizes=lambda wc: REF["chrom_sizes"]
    output:
        "results/peaks/consensus_peaks.bed"
    log:
        "results/logs/merge_factor_consensus.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/peaks results/logs
        cat {input.peaks:q} \
          | bedtools sort -i - -g {input.chrom_sizes:q} 2> {log:q} \
          | bedtools merge > {output:q} 2>> {log:q}
        """
