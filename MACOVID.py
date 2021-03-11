#!/usr/bin/env python3
import shutil
from argparse import ArgumentParser
import yaml
import os
import pandas as pd
import re
import glob
import sys


locationrepo = os.path.dirname(os.path.abspath(__file__)) 

def get_absolute_path(path):
    return os.path.abspath(path)

def file_name_generator(filepath):
    return os.path.splitext(os.path.basename(filepath))[0]

###################
# Define the input files
###################

def define_input(inputdir):
    ## Input dir treated as list. Remove square brackets
    if bool(re.search("\[*\]", str(inputdir))):
         inputdir = str(inputdir)[2:-2]

    ## Remove last slash if existed in the input folder
    if str(inputdir)[-1] == "/":
        inputdir = inputdir[:-1]

    ## Path of folder with fastq files
    inputPath = get_absolute_path(inputdir)

    ## Fastq files
    fastqFiles = glob.glob(f"{inputPath}/*fastq*")

    ## Remove unclassified files
    fastqFiles[:] = [x for x in fastqFiles if "unclassified" not in x]

    ## Check fastq files
    if len(fastqFiles) == 0:
        print ("Check directory no fastq file found")
        sys.exit()

    return fastqFiles

####################
# changing the names of barcoded samples
####################

def change_names(sampledir, manifest, reverse):
    """
    take samples and change name according to the manifest file.
    if -rev flag is used, the samples are converted back 
    """
    print(f'reverse is {reverse}')
    if reverse:
        pos = 0
    else:
        pos = 1 
    ## Read csv files:
    names = pd.read_csv(manifest, index_col = pos, sep = ",|;|\t", engine = 'python').dropna().to_dict()

#    print (names)
    ## Key is the sample_ID in this nested dict
    key = [x for x in names.keys()][0]
    barcodes = names[key]

#    print (barcodes)
    ## Identify fastq files
    fastqFiles = define_input(sampledir)
#    print (fastqFiles)

    runFiles = []
    ## Loop through fastq file found
    for q in fastqFiles:
#        print (q)
        ## Loop through barcode dictionary
        for k,v in barcodes.items():
#            print (q,k,v)
            basename = q.split("/")[-1].split(".fastq")[0]
#            print ("fastq",basename,k, v)
            if basename == k:
#                print (basename, v)
                location = os.path.dirname(q)
#                print (location)
                newFile = location + "/" + f"{v}.fastq"
#                print (newFile)
                shutil.move(q, newFile)
                runFiles.append(newFile)

    return runFiles
###################
# Snakemake pipeline for generating consensus fastas
###################

def snakemake_in(samplesin, outdir, coverage):

    ## Sample dictionary
    samplesdic = {}
    ## parameter for the ourDir
    samplesdic['parameters'] = {}
    samplesdic['parameters']["outdir"] = get_absolute_path(outdir)
    samplesdic['parameters']["coverage"] = coverage
    samplesdic["SAMPLES"] = {}
    
    # generate the samples dictionary as input for snakemake 
    for i in samplesin:
        samplename = file_name_generator(i)
        samplesdic["SAMPLES"][samplename] = get_absolute_path(i)
    data = yaml.dump(samplesdic, default_flow_style=False)
    
    # make and write config file location
    os.system(f"mkdir -p {locationrepo}/config")
    with open(f"{locationrepo}/config/config.yaml", 'w') as f:
        f.write(data)

####################
# Command line Parsers initialization
####################

def main(command_line = None):
    #add main parser object
    parser = ArgumentParser(description = "COVID seq toolkit...")

    #add sub parser object
    subparsers = parser.add_subparsers(dest = "mode")
    #add snakemake pipeline to completely run fasta to clustered output
    mapreads = subparsers.add_parser("mapreads", help = "run full pipeline from fastq to consensus sequence")
    mapreads.add_argument("-i", required = True, dest = "input_directory", nargs = "+", help = 'give fastq files, and use basename to shoot into the pipeline')
    mapreads.add_argument("--cores", dest = 'cores', required = True, type = int, help = 'Number of CPU cores to use')
    mapreads.add_argument("-o", required = True, dest = "outdir")
    mapreads.add_argument("-m", required = True, dest = "manifest")
    mapreads.add_argument("-l", required = False, dest = "local", action = "store_false")
    mapreads.add_argument("--cov", required = False, dest = "coverage", default = 30, type = int, help = "min coverage per base (default: 30)")

    namechange = subparsers.add_parser("namechanger", help = "change barcode names")
    namechange.add_argument("-i", required = True, nargs = "+", dest = "input_directory")
    namechange.add_argument("-m", required = True, dest = "manifest")
    namechange.add_argument("-rev", required = False, dest = "rev", action = "store_true") 

    rerun = subparsers.add_parser("rerun", help = "Rerun samples")
    rerun.add_argument("-i", required = True, nargs = "+", dest = "input_directory")
    rerun.add_argument("-o", required = True, dest = "outdir")
    rerun.add_argument("--cores", dest = 'cores', required = True, type = int, help = 'Number of CPU cores to use')
    rerun.add_argument("-l", required = False, dest = "local", action = "store_false")
    rerun.add_argument("--cov", required = False, dest = "coverage", default = 30, type = int, help = "min coverage per base (default: 30)")

####################
# parsing part
####################

    args = parser.parse_args(command_line)
    if args.mode == "mapreads":
        ## Set running location where MACOVID is located
        os.chdir(f"{locationrepo}")

        ## Change name according to manifest
        samplesin = change_names(
                sampledir = args.input_directory,
                manifest = args.manifest,
                reverse ="FALSE",
                coverage = args.coverage
                )
        snakemake_in(
                samplesin = samplesin,
                outdir = args.outdir
                )
        if not args.local:
                os.system(f"snakemake --cores {args.cores} --use-conda --latency-wait 30 -k")
        else:
                os.system(f"snakemake --cluster 'sbatch --output=/dev/null' --jobs 100 --latency-wait 90 --cores {args.cores} --use-conda -k")
    elif args.mode == "namechanger":

       change_names(
                sampledir = args.input_directory,
                manifest = args.manifest,
                reverse = args.rev 
                )
    elif args.mode == "rerun":

        samplesin = define_input(
                inputdir = args.input_directory,
                )
        snakemake_in(
                samplesin = samplesin,
                outdir = args.outdir,
                coverage = args.coverage
                )
        if not args.local:
                os.system(f"snakemake --cores {args.cores} --use-conda --latency-wait 30 -k")
        else:
                os.system(f"snakemake --cluster 'sbatch --output=/dev/null' --jobs 100 --latency-wait 90 --cores {args.cores} --use-conda -k")
    else:
        parser.print_usage()

if __name__ == "__main__":
    main()
