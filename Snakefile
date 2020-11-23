

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']
OUTDIR = config['parameters']['outdir'] + "/"


rule all:
    input:
        expand(OUTDIR + "{sample}_final.fasta", sample = SAMPLES),

localrules: combine, 

rule trimming:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        temp(OUTDIR + "{sample}_trimmed.fastq")
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
        temp(OUTDIR + "{sample}_mapped.bam")
    conda:
        "envs/refmap.yaml"
    threads: 5
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {input.ref} {input.trim} | samtools view -bF 4 - | samtools sort -@ {threads} - > {output}
        """

rule bamIndex:
    input:
        OUTDIR + "{sample}_mapped.bam"
    output:
        temp(OUTDIR + "{sample}_mapped.bam.bai")
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
        OUTDIR + "{sample}_consensus.fasta"
    conda:
        "envs/bam2concensus.yaml"
    params: 30
    shell:
        """
        bin/bam2consensus.py -i {input.bamfile} -o {output} -d {params} -g 1
        """

rule align_consensus:
    input:
        consen=OUTDIR + "{sample}_consensus.fasta",
        ref="reference.fasta"
    output:
        OUTDIR + "{sample}_final.fasta"
    conda:
        "envs/alignconcensus.yaml"
    shell:
        """
        bin/align_to_ref.py -i {input.consen} -o {output} -r {input.ref} -n {wildcards.sample}
        """

