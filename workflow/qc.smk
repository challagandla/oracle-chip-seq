rule fastqc_raw:
    input:
        R1=lambda wc: FASTQ[wc.sample][0],
        R2=lambda wc: FASTQ[wc.sample][1]
    output:
        "results/fastqc/raw/{sample}_R1_fastqc.html",
        "results/fastqc/raw/{sample}_R1_fastqc.zip",
        "results/fastqc/raw/{sample}_R2_fastqc.html",
        "results/fastqc/raw/{sample}_R2_fastqc.zip"
    threads: 2
    log:
        "results/logs/fastqc_raw_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/fastqc/raw results/logs
        tmpdir=$(mktemp -d results/fastqc/raw/.{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        ln -s "$(readlink -f {input.R1:q})" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -s "$(readlink -f {input.R2:q})" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        fastqc -o results/fastqc/raw -t {threads} \
            "$tmpdir/{wildcards.sample}_R1.fastq.gz" \
            "$tmpdir/{wildcards.sample}_R2.fastq.gz" > {log:q} 2>&1
        """


rule fastq_screen_raw:
    input:
        R1=lambda wc: FASTQ[wc.sample][0],
        R2=lambda wc: FASTQ[wc.sample][1],
        conf=FASTQ_SCREEN_CONF
    output:
        R1txt="results/contamination/fastq_screen/{sample}_R1_screen.txt",
        R1html="results/contamination/fastq_screen/{sample}_R1_screen.html",
        R2txt="results/contamination/fastq_screen/{sample}_R2_screen.txt",
        R2html="results/contamination/fastq_screen/{sample}_R2_screen.html"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.R1txt)),
        subset=FASTQ_SCREEN_SUBSET
    threads: 4
    log:
        "results/logs/fastq_screen_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        tmpdir=$(mktemp -d {params.outdir:q}/tmp.{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        ln -sf "$(readlink -f {input.R1:q})" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -sf "$(readlink -f {input.R2:q})" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        fastq_screen --conf {input.conf:q} --aligner bowtie2 --threads {threads} \
            --subset {params.subset} --outdir {params.outdir:q} \
            "$tmpdir/{wildcards.sample}_R1.fastq.gz" \
            "$tmpdir/{wildcards.sample}_R2.fastq.gz" > {log:q} 2>&1
        """


rule trim_galore:
    input:
        R1=lambda wc: FASTQ[wc.sample][0],
        R2=lambda wc: FASTQ[wc.sample][1]
    output:
        trimmed1="results/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed2="results/trimmed/{sample}_R2_val_2.fq.gz"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.trimmed1))
    threads: 4
    log:
        "results/logs/trim_galore_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        tmpdir=$(mktemp -d {params.outdir:q}/.{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        ln -s "$(readlink -f {input.R1:q})" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -s "$(readlink -f {input.R2:q})" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        trim_galore --paired --cores {threads} --output_dir {params.outdir:q} \
            "$tmpdir/{wildcards.sample}_R1.fastq.gz" \
            "$tmpdir/{wildcards.sample}_R2.fastq.gz" > {log:q} 2>&1
        """


rule fastqc_trimmed:
    input:
        R1="results/trimmed/{sample}_R1_val_1.fq.gz",
        R2="results/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        "results/fastqc/trimmed/{sample}_R1_val_1_fastqc.html",
        "results/fastqc/trimmed/{sample}_R1_val_1_fastqc.zip",
        "results/fastqc/trimmed/{sample}_R2_val_2_fastqc.html",
        "results/fastqc/trimmed/{sample}_R2_val_2_fastqc.zip"
    threads: 2
    log:
        "results/logs/fastqc_trimmed_{sample}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/fastqc/trimmed results/logs
        fastqc -o results/fastqc/trimmed -t {threads} {input.R1} {input.R2} > {log} 2>&1
        """
