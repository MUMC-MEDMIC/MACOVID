

configfile: "config/config.yaml"
SAMPLES = config['SAMPLES']
OUTDIR = config['parameters']['outdir'] + "/"


rule all:
    input:
        expand("{sample}_consensus.fasta", sample = SAMPLES)
