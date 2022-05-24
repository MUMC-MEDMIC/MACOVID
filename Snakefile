"""
Pipeline to analyse COVID variants
MACOVID version 2.2.0
MUMC+
"""

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']
OUTDIR = config['parameters']['outdir'] + "/"
MERGE =  config['parameters']['merge_files']
COVERAGE = config['parameters']['coverage']
SCHEMEDIR = config['parameters']['scheme'] + "/"
SCHEMEPREFIX = config['parameters']['schemePrefix']
MINLENGTH = config['parameters']['minLength']
MAXLENGTH = config['parameters']['maxLength']
MAJORITY = config['parameters']['majority']


rule all:
    input:
        expand(OUTDIR + "{sample}-barplot.png", sample = SAMPLES),
        expand(OUTDIR + "{sample}-boxplot.png", sample = SAMPLES),
        expand(OUTDIR + "{sample}.consensus.fasta", sample = SAMPLES),
        expand(OUTDIR + "{sample}.primer.vcf",sample = SAMPLES)

rule readfilter:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        OUTDIR + "{sample}_trimmed.fastq"
    conda:
        "envs/cutadapt.yaml"
    params:
        minLength = MINLENGTH,
        maxLength = MAXLENGTH
    threads: 4
    shell:
        """
        cutadapt -o {output} {input} -m {params.minLength} -M {params.maxLength} -j {threads}
        """


rule mapping:
    input:
        trim = OUTDIR + "{sample}_trimmed.fastq",
    output:
        OUTDIR + "{sample}_mapped.bam"
    conda:
        "envs/refmap.yaml"
    threads: 4
    params:
        ref = SCHEMEDIR + SCHEMEPREFIX + ".reference.fasta"
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {params.ref} {input.trim} | samtools view -bF 4 - | samtools sort -@ {threads} - > {output}
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


rule trimAlignment:
    input:
        bamfile = OUTDIR + "{sample}_mapped.bam",
        bamindex = OUTDIR + "{sample}_mapped.bam.bai"
    output:
        report = OUTDIR + "{sample}.alignreport.txt",
        dropped = OUTDIR + "{sample}.alignreport.er",
        trimmedBamfile = OUTDIR + "{sample}.trimmed.rg.sorted.bam",
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam"
    conda:
        "envs/artic.yaml"
    params:
        prefix = "{sample}", 
        scheme = SCHEMEDIR + SCHEMEPREFIX + ".scheme.bed"
    shell:
        """
        align_trim --start --normalise 200 {params.scheme} --report {output.report} < {input.bamfile} 2> {output.dropped} | samtools sort -T {params.prefix} - -o {output.trimmedBamfile};
        align_trim --normalise 200 {params.scheme} --remove-incorrect-pairs --report {output.report} < {input.bamfile} 2> {output.dropped} | samtools sort -T {params.prefix} - -o {output.primertrimmedBamfile}
        """


rule indexTrimmedBam:
    input:
        trimmedBamfile = OUTDIR + "{sample}.trimmed.rg.sorted.bam",
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam"
    output:
        trimmedBamfileIndex = OUTDIR + "{sample}.trimmed.rg.sorted.bam.bai",
        primertrimmedBamfileIndex = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam.bai"
    conda:
        "envs/refmap.yaml"
    threads: 1
    shell:
        """
        samtools index -@ {threads} {input.trimmedBamfile};
        samtools index -@ {threads} {input.primertrimmedBamfile} 
        """


rule splitPrimerpool:
    input: 
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam",
        primertrimmedBamfileIndex = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam.bai"
    output:
        pool1Bam = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_1.sorted.bam",
        pool2Bam = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_2.sorted.bam",
        pool1Index = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_1.sorted.bam.bai",
        pool2Index = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_2.sorted.bam.bai"
    params:
        pool1 = SCHEMEPREFIX + "_1",
        pool2 = SCHEMEPREFIX + "_2"
    conda:
        "envs/refmap.yaml"
    threads: 1
    shell:
        """
        samtools view -b -r {params.pool1} {input.primertrimmedBamfile} > {output.pool1Bam};
        samtools index {output.pool1Bam};
        samtools view -b -r {params.pool2} {input.primertrimmedBamfile} > {output.pool2Bam};
        samtools index {output.pool2Bam}
        """


rule medakaConsensus:
    input:
        poolBam = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_{num}.sorted.bam",
        poolIndex = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_{num}.sorted.bam.bai"
    output:
        poolHdf = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_{num}.hdf"
    conda:
        "envs/medaka.yaml"
    threads: 1
    shell:
        """
        medaka consensus --model r941_min_hac_g507 {input.poolBam} {output.poolHdf};
        """


rule medakaVariant:
    input:
        poolHdf = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_{num}.hdf"
    output:
        poolVcf = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_{num}.vcf"
    conda:
        "envs/medaka.yaml"
    threads: 1
    params:
        ref = SCHEMEDIR + SCHEMEPREFIX + ".reference.fasta"
    shell:
        """
        medaka variant {params.ref} {input.poolHdf} {output.poolVcf};
        """
        

rule vcfMerge:
    input:
        pool1Vcf = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_1.vcf",
        pool2Vcf = OUTDIR + "{sample}.primertrimmed." + SCHEMEPREFIX + "_2.vcf"
    output:
        mergedVcf = OUTDIR + "{sample}.merged.vcf"
    conda:
        "envs/artic.yaml"
    threads: 1
    params: 
        prefix = "{sample}",
        outdir = OUTDIR,
        scheme = os.getcwd() + "/" + SCHEMEDIR + SCHEMEPREFIX + ".scheme.bed",
        schemeprefix = SCHEMEPREFIX
    shell:
        """
        cd {params.outdir}
        artic_vcf_merge {params.prefix} {params.scheme} {params.schemeprefix}_1:{input.pool1Vcf} {params.schemeprefix}_2:{input.pool2Vcf} 
        cd -
        """


rule vcfBGzip:
    input: OUTDIR + "{sample}.merged.vcf"
    output: OUTDIR + "{sample}.merged.vcf.gz"
    conda: "envs/artic.yaml"
    shell:
        """
        bgzip -f {input}
        """


rule vcfMergeTabix:
    input: OUTDIR + "{sample}.merged.vcf.gz"
    output: OUTDIR + "{sample}.merged.vcf.gz.tbi"
    conda: "envs/artic.yaml"
    shell:
        """
        tabix -p vcf {input}
        """


rule longshot:
    input:
        mergedVcfGz = OUTDIR + "{sample}.merged.vcf.gz",
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam",
        primertrimmedBamfileIndex = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam.bai",
        mergedTabix = OUTDIR + "{sample}.merged.vcf.gz.tbi"
    output: 
        recall = OUTDIR + "{sample}.recall.vcf",
        longshot = OUTDIR + "{sample}.longshot.vcf"
    conda:
        "envs/longshot.yaml"
    threads: 1
    params: 
        prefix = OUTDIR + "{sample}",
        ref = SCHEMEDIR + SCHEMEPREFIX + ".reference.fasta"
    shell:
        """
        longshot -P 0 -F -c 30 -C 850 --no_haps --bam {input.primertrimmedBamfile} --ref {params.ref} --out {output.recall} --potential_variants {params.prefix}.merged.vcf.gz;
        longshot -P 0 -F -c 30 -C 850 --no_haps --bam {input.primertrimmedBamfile} --ref {params.ref} --out {output.longshot} 
        """

rule vcfPreparation:
    input:
        medakaVcf = OUTDIR + "{sample}.merged.vcf.gz",
        recallVcf = OUTDIR + "{sample}.recall.vcf",
        longshotVcf = OUTDIR + "{sample}.longshot.vcf"
    output: 
        vcfPass = OUTDIR + "{sample}.pass.vcf",
        vcfFail = OUTDIR + "{sample}.fail.vcf"
    conda:
        "envs/longGapVcf.yaml"
    params: 
        majority = MAJORITY
    threads: 1
    shell:
        """
        Rscript scripts/long_gap_vcf.R {input.medakaVcf} {input.recallVcf} {input.longshotVcf} {params.majority} {output.vcfPass} {output.vcfFail}
        """

rule primerMutations:
    input: 
         vcfPass = OUTDIR + "{sample}.pass.vcf",
         vcfFail = OUTDIR + "{sample}.fail.vcf"
    output: 
         primerMutationVcf = OUTDIR + "{sample}.primer.vcf"
    conda:
        "envs/bedtools.yaml"
    params: 
        prefix = "{sample}", 
        scheme = SCHEMEDIR + SCHEMEPREFIX + ".scheme.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {params.scheme} -b {input.vcfPass} -wa -wb -header > {output.primerMutationVcf}
        """


rule depthMask:
    input:
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam",
        primertrimmedBamfileIndex = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam.bai"
    output: OUTDIR + "{sample}.coverage_mask.txt"
    conda:
        "envs/artic.yaml"
    params: 
        coverage = COVERAGE,
        ref = SCHEMEDIR + SCHEMEPREFIX + ".reference.fasta"
    threads: 1
    shell:
        """
        artic_make_depth_mask --depth {params.coverage} --store-rg-depths {params.ref} {input.primertrimmedBamfile} {output}
        """


rule plotAmpliconDepth:
    input: 
        maskdepth = OUTDIR + "{sample}.coverage_mask.txt",
        mergedTabix = OUTDIR + "{sample}.merged.vcf.gz.tbi"
    output:
        barplot = OUTDIR + "{sample}-barplot.png",
        boxplot = OUTDIR + "{sample}-boxplot.png"
    conda:
        "envs/artic.yaml"
    threads: 1
    params: 
        prefix = "{sample}",
        scheme = os.getcwd() + "/" + SCHEMEDIR + SCHEMEPREFIX + ".scheme.bed",
        outDir = OUTDIR,
    shell:
        """
        cd {params.outDir}
        artic_plot_amplicon_depth --primerScheme {params.scheme} --sampleID {params.prefix} --outFilePrefix {params.prefix} {params.prefix}*.depths
        cd -
        rm {params.outDir}{params.prefix}*.depths
        """


rule vcfPassBGzip:
    input: vcfPass = OUTDIR + "{sample}.pass.vcf",
           primerVcf = OUTDIR + "{sample}.primer.vcf"
    output: OUTDIR + "{sample}.pass.vcf.gz"
    conda: "envs/artic.yaml"
    shell:
        """
        bgzip -f {input.vcfPass}
        """


rule preconsensus:
    input: 
        mask = OUTDIR + "{sample}.coverage_mask.txt",
        vcfPassGz = OUTDIR + "{sample}.pass.vcf.gz",
        vcfFail = OUTDIR + "{sample}.fail.vcf"
    output: 
        preconsensus = OUTDIR + "{sample}.preconsensus.fasta",
        vcfPassGzIndex = OUTDIR + "{sample}.pass.vcf.gz.tbi"
    conda:
        "envs/artic.yaml"
    threads: 1
    params:
        ref = SCHEMEDIR + SCHEMEPREFIX + ".reference.fasta"
    shell:
        """
        tabix -p vcf {input.vcfPassGz};
        artic_mask {params.ref} {input.mask} {input.vcfFail} {output.preconsensus}
        """


rule consensus:
    input: 
        preconsensus = OUTDIR + "{sample}.preconsensus.fasta",  
        vcfPassGz = OUTDIR + "{sample}.pass.vcf.gz",
        mask = OUTDIR + "{sample}.coverage_mask.txt",
        vcfFail = OUTDIR + "{sample}.fail.vcf",
        vcfPassGzIndex = OUTDIR + "{sample}.pass.vcf.gz.tbi"

    output: OUTDIR + "{sample}.consensus.fasta"
    conda:
        "envs/artic.yaml"
    threads: 1
    params: "{sample}"
    shell:
        """
        bcftools consensus -f {input.preconsensus} {input.vcfPassGz} -m {input.mask} -o {output};
        artic_fasta_header {output} {params}
        """