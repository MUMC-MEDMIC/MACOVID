

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']
OUTDIR = config['parameters']['outdir'] + "/"
COVERAGE = config['parameters']['coverage']

rule all:
    input:
        expand(OUTDIR + "{sample}.consensus.fasta", sample = SAMPLES),
        expand(OUTDIR + "{sample}.helpfile2.txt", sample = SAMPLES)

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
        cutadapt -o {output} {input} -m 75 -j {threads}
        """

rule mapping:
    input:
        trim = OUTDIR + "{sample}_trimmed.fastq",
        ref = "primer_schemes/nCoV-2019.reference.fasta"
    output:
        OUTDIR + "{sample}_mapped.bam"
    conda:
        "envs/refmap.yaml"
    threads: 4
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

rule trimAlignment:
    input:
        scheme = "primer_schemes/nCoV-2019.scheme.bed",
        bamfile = OUTDIR + "{sample}_mapped.bam",
        sortfile = OUTDIR + "{sample}_mapped.bam.bai"
    output:
        report = OUTDIR + "{sample}.alignreport.txt",
        dropped = OUTDIR + "{sample}.alignreport.er",
        trimmedBamfile = OUTDIR + "{sample}.trimmed.rg.sorted.bam",
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam"
    conda:
        "envs/artic.yaml"
    params: "{sample}"  
    shell:
        """
        align_trim --start --normalise 200 {input.scheme} --report {output.report} < {input.bamfile} 2> {output.dropped} | samtools sort -T {params} - -o {output.trimmedBamfile};
        align_trim --normalise 200 {input.scheme} --remove-incorrect-pairs --report {output.report} < {input.bamfile} 2> {output.dropped} | samtools sort -T {params} - -o {output.primertrimmedBamfile}
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
        pool1Bam = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.sorted.bam",
        pool2Bam = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.sorted.bam",
        pool1Index = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.sorted.bam.bai",
        pool2Index = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.sorted.bam.bai"
    conda:
        "envs/refmap.yaml"
    threads: 1
    shell:
        """
        samtools view -b -r "nCoV-2019_1" {input.primertrimmedBamfile} > {output.pool1Bam};
        samtools index {output.pool1Bam};
        samtools view -b -r "nCoV-2019_2" {input.primertrimmedBamfile} > {output.pool2Bam};
        samtools index {output.pool2Bam}
        """

rule medakaConsensus:
    input:
        pool1Bam = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.sorted.bam",
        pool2Bam = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.sorted.bam",
        pool1Index = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.sorted.bam.bai",
        pool2Index = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.sorted.bam.bai"
    output:
        pool1Hdf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.hdf",
        pool2Hdf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.hdf"
    conda:
        "envs/artic.yaml"
    threads: 1
    shell:
        """
        medaka consensus --chunk_len 400 --chunk_ovlp 200 {input.pool1Bam} {output.pool1Hdf};
        medaka consensus --chunk_len 400 --chunk_ovlp 200 {input.pool2Bam} {output.pool2Hdf}
        """

rule medakaVariant:
    input:
        ref = "primer_schemes/nCoV-2019.reference.fasta",  
        pool1Hdf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.hdf",
        pool2Hdf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.hdf"
    output:
        pool1Vcf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.vcf",
        pool2Vcf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.vcf"
    conda:
        "envs/artic.yaml"
    threads: 1
    shell:
        """
        medaka variant {input.ref} {input.pool1Hdf} {output.pool1Vcf};
        medaka variant {input.ref} {input.pool2Hdf} {output.pool2Vcf}
        """
        
rule vcfMerge:
    input:
        scheme = "primer_schemes/nCoV-2019.scheme.bed",  
        pool1Vcf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_1.vcf",
        pool2Vcf = OUTDIR + "{sample}.primertrimmed.nCoV-2019_2.vcf"
    output: OUTDIR + "{sample}.helpfile.txt"
    conda:
        "envs/artic.yaml"
    threads: 1
    params: 
        prefix = "{sample}",
        outDir = OUTDIR,
        curDir = os.getcwd() + "/"
    shell:
        """
        bash vcfMerge.sh {input.scheme} {input.pool1Vcf} {input.pool2Vcf} {params.curDir} {params.outDir} {params.prefix};
        echo "ok" > {output}
        """

rule longshot:
    input:
        ref = "primer_schemes/nCoV-2019.reference.fasta",  
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam",
        primertrimmedBamfileIndex = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam.bai",
        helpfile = OUTDIR + "{sample}.helpfile.txt"
    output: OUTDIR + "{sample}.longshot.vcf"
    conda:
        "envs/artic.yaml"
    threads: 1
    params: 
        prefix = OUTDIR + "{sample}"
    shell:
        """
        longshot -P 0 -F -A --no_haps --bam {input.primertrimmedBamfile} --ref {input.ref} --out {output} --potential_variants {params.prefix}.merged.vcf.gz
        """

rule vcfFilter:
    input: OUTDIR + "{sample}.longshot.vcf"
    output: 
        vcfPass = OUTDIR + "{sample}.pass.vcf",
        vcfFail = OUTDIR + "{sample}.fail.vcf"
    conda:
        "envs/artic.yaml"
    threads: 1
    shell:
        """
        artic_vcf_filter --longshot {input} {output.vcfPass} {output.vcfFail}
        """

rule depthMask:
    input:
        ref = "primer_schemes/nCoV-2019.reference.fasta",  
        primertrimmedBamfile = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam",
        primertrimmedBamfileIndex = OUTDIR + "{sample}.primertrimmed.rg.sorted.bam.bai"
    output: OUTDIR + "{sample}.coverage_mask.txt"
    conda:
        "envs/artic.yaml"
    params: COVERAGE
    threads: 1
    shell:
        """
        artic_make_depth_mask --depth {params} --store-rg-depths {input.ref} {input.primertrimmedBamfile} {output}
        """

rule plotAmpliconDepth:
    input: 
        scheme = "primer_schemes/nCoV-2019.scheme.bed",
        helpfile = OUTDIR + "{sample}.helpfile.txt"
    output: OUTDIR + "{sample}.helpfile2.txt"
    conda:
        "envs/artic.yaml"
    threads: 1
    params: 
         prefix = "{sample}",
         outDir = OUTDIR,
         curDir = os.getcwd() + "/"
    shell:
        """
        bash plotAmplicon.sh {input.scheme} {params.prefix} {params.curDir} {params.outDir};
        echo "ok" > {output}
        """

rule preconsensus:
    input: 
        ref = "primer_schemes/nCoV-2019.reference.fasta",  
        mask = OUTDIR + "{sample}.coverage_mask.txt",
        vcfPass = OUTDIR + "{sample}.pass.vcf",
        vcfFail = OUTDIR + "{sample}.fail.vcf"

    output: 
        vcfPassGz = OUTDIR + "{sample}.pass.vcf.gz",
        preconsensus = OUTDIR + "{sample}.preconsensus.fasta"
    conda:
        "envs/artic.yaml"
    threads: 1
    shell:
        """
        bgzip -f {input.vcfPass};
        sleep 45s;
        tabix -p vcf {output.vcfPassGz};
        artic_mask {input.ref} {input.mask} {input.vcfFail} {output.preconsensus}
        """

rule consensus:
    input: 
        preconsensus = OUTDIR + "{sample}.preconsensus.fasta",  
        vcfPassGz = OUTDIR + "{sample}.pass.vcf.gz",
        mask = OUTDIR + "{sample}.coverage_mask.txt",
        vcfFail = OUTDIR + "{sample}.fail.vcf"

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
