

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']
OUTDIR = config['parameters']['outdir'] + "/"


rule all:
    input:
        expand(OUTDIR + "{sample}_trimmed.fastq", sample = SAMPLES),
        expand(OUTDIR + "{sample}_mapped.bam", sample = SAMPLES),
        expand(OUTDIR + "{sample}_mapped.bam.bai", sample = SAMPLES),
        expand(OUTDIR + "{sample}_consencus.fasta", sample = SAMPLES),
        expand(OUTDIR + "{sample}_fin.fasta", sample = SAMPLES),

rule trimming:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        OUTDIR + "{sample}_trimmed.fastq"
    conda:
        "envs/cutadapt.yaml"
    threads: 4
    shell:
        """
        cutadapt -u 30 -u -30 -o {output} {input} -m 75 -j {threads}
        """

rule mapping:
    input:
        trim=OUTDIR + "{sample}_trimmed.fastq",
        ref="reference.fasta"
    output:
        OUTDIR + "{sample}_mapped.bam"
    conda:
        "envs/refmap.yaml"
    threads: 10
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {input.ref} {input.trim} | samtools view -bF 4 - | samtools sort -@ {threads} - > {output}
        """

rule bamIndex:
    input:
        OUTDIR + "{sample}_mapped.bam"
    output:
        OUTDIR + "{sample}_mapped.bam.bai"
    conda:
        "envs/refmap.yaml"
    threads: 1
    shell:
        """
        samtools index -@ {threads} {input} 
        """

rule bam2concensus:
    input:
        bamfile = OUTDIR + "{sample}_mapped.bam",
        sortfile = OUTDIR + "{sample}_mapped.bam.bai"
    output:
        OUTDIR + "{sample}_consencus.fasta"
    conda:
        "envs/bam2concencus.yaml"
    shell:
        """
        bin/bam2consensus.py -i {input.bamfile} -o {output} -d 30 -g 1
        """

rule align_consensus:
    input:
        consen=OUTDIR + "{sample}_consencus.fasta",
        ref="reference.fasta"
    output:
        OUTDIR + "{sample}_fin.fasta"
    conda:
        "envs/alignconcencus.yaml"
    shell:
        """
        bin/align_to_ref.py -i {input.consen} -o {output} -r {input.ref} -n {wildcards.sample}
        """
