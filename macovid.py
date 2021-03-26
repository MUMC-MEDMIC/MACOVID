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

def define_inDir(inputdir):
    ## Input dir treated as list. Remove square brackets
    if bool(re.search("\[*\]", str(inputdir))):
         inputdir = str(inputdir)[2:-2]

    ## Remove last slash if existed in the input folder
    if str(inputdir)[-1] == "/":
        inputdir = inputdir[:-1]

    return (inputdir)
    
def define_inFastq(inputDir):
    ## Path of folder with fastq files
    inputPath = get_absolute_path(inputDir)

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

    ## Read csv files into nested doct:
    names = pd.read_csv(manifest, index_col = pos, sep = ",|;|\t", engine = 'python').dropna().to_dict()

    ## Extract keys (sample_ID) from manifest file
    key = [x for x in names.keys()][0]

    ## Extract the values (barcodes and new names) from dict
    barcodes = names[key]

    ## fastq folder
    folderLoc = define_inDir(sampledir)

    ## Identify fastq files
    fastqFiles = define_inFastq(folderLoc)

    runFiles = []

    ## Incase of float converts back to string
    formatNumber = lambda n: int(n) if isinstance(n,float) and n.is_integer() else n

    ## Loop through found fastq files
    for fstq in fastqFiles:
        ## Check through barcode dictionary
        for bark,barv in barcodes.items():
            
            basename = fstq.split("/")[-1].split(".fastq")[0]
            ## Basename found on the barcode dictionary
            if basename == bark:
                ## Solution if the new name is numeric or float
                if barv.isnumeric():
                    barv = formatNumber(barv)

                ## Rename of found fastq file
                location = os.path.dirname(fstq)
                newFile = location + "/" + f"{barv}.fastq"
                shutil.move(fstq, newFile)
                ## Writes run files to new list
                runFiles.append(newFile)

    return runFiles

###################
# Snakemake pipeline for generating consensus fastas
###################

def snakemake_in(samplesin, folderin,  outdir, coverage, scheme, schemePrefix):

    ## Sample dictionary
    samplesdic = {}
    ## parameter for the ourDir
    samplesdic['parameters'] = {}
    samplesdic['parameters']['merge_files'] = folderin
    samplesdic['parameters']["outdir"] = get_absolute_path(outdir)
    samplesdic['parameters']["coverage"] = coverage
    samplesdic['parameters']["scheme"] = scheme
    samplesdic['parameters']["schemePrefix"] = schemePrefix
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
    mapreads.add_argument("--trim_start", required = False, dest = "trimStart", default = 54, type = int, help = "trim start of consensus that is not sequenced (default: 54)")	
    mapreads.add_argument("--trim_end", required = False, dest = "trimEnd", default = 79, type = int, help = "trim end of consensus that is not sequenced (default: 79)")
    mapreads.add_argument("--scheme", required = False, dest = "schemeDir", default = "primer_schemes/EMC", type = str, help = "path to scheme directory (default: primer_schemes/EMC)")
    mapreads.add_argument("--scheme_prefix", required = False, dest = "schemePrefix", default = "nCoV-2019", type = str, help = "prefix of primer scheme (default: nCoV-2019)")


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
    rerun.add_argument("--trim_start", required = False, dest = "trimStart", default = 54, type = int, help = "trim start of consensus that is not sequenced (default: 54)")	
    rerun.add_argument("--trim_end", required = False, dest = "trimEnd", default = 79, type = int, help = "trim end of consensus that is not sequenced (default: 79)")
    rerun.add_argument("--scheme", required = False, dest = "schemeDir", default = "primer_schemes/EMC", type = str, help = "path to scheme directory (default: primer_schemes/EMC)")
    rerun.add_argument("--scheme_prefix", required = False, dest = "schemePrefix", default = "nCoV-2019", type = str, help = "prefix of primer scheme (default: nCoV-2019)")


####################
# parsing part
####################

    args = parser.parse_args(command_line)
    if args.mode == "mapreads":
        ## Set running location where MACOVID is located
        os.chdir(f"{locationrepo}")
        ## Identify fastq input files, change name according to manifest
        samplesIn = change_names(
                sampledir = args.input_directory,
                manifest = args.manifest,
                reverse ="FALSE",
                )

        ## fastq folder
        folderLoc = define_inDir(inputdir = args.input_directory)

        snakemake_in(
                samplesin = samplesIn,
                folderin = folderLoc,
                outdir = args.outdir,
                coverage = args.coverage,
                scheme = args.schemeDir,
                schemePrefix = args.schemePrefix
                )
        if not args.local:
                print ("Running MACOVID locally")
                os.system(f"snakemake --cores {args.cores} --use-conda --latency-wait 30 -k -p ")
                os.system(f"cat {args.outdir}/*.consensus.fasta | cutadapt -u {args.trimStart} -u -{args.trimEnd} - > {args.outdir}/merged_trimmed.fasta")
        else:
                print ("Running MACOVID on the cluster")
                os.system(f"snakemake --cluster 'sbatch --output=/dev/null' --jobs 100 --latency-wait 90 --cores {args.cores} --use-conda -k -p ")
                os.system(f"cat {args.outdir}/*.consensus.fasta | cutadapt -u {args.trimStart} -u -{args.trimEnd} - > {args.outdir}/merged_trimmed.fasta")

    elif args.mode == "namechanger":

       change_names(
                sampledir = args.input_directory,
                manifest = args.manifest,
                reverse = args.rev 
                )

    ## Rerun fastq files from specific folder
    elif args.mode == "rerun":

        ## Identify fastq folder
        folderLoc = define_inDir(inputdir = args.input_directory)

        ## Detect fastq files inside folder
        samplesIn = define_inFastq(inputDir = folderLoc)

        ## Rerun commands
        snakemake_in(
                samplesin = samplesIn,
                folderin = folderLoc,
                outdir = args.outdir,
                coverage = args.coverage,
                scheme = args.schemeDir,
                schemePrefix = args.schemePrefix
                )
 
        if not args.local:
                print ("Re-running MACOVID locally")
                os.system(f"snakemake --cores {args.cores} --use-conda --latency-wait 30 -k -p ")
                os.system(f"cat {args.outdir}/*.consensus.fasta | cutadapt -u {args.trimStart} -u -{args.trimEnd} - > {args.outdir}/merged_trimmed.fasta")
        else:
                print ("Re-running MACOVID on the cluster")
                os.system(f"snakemake --cluster 'sbatch --output=/dev/null' --jobs 100 --latency-wait 90 --cores {args.cores} --use-conda -k -p ")
                os.system(f"cat {args.outdir}/*.consensus.fasta | cutadapt -u {args.trimStart} -u -{args.trimEnd} - > {args.outdir}/merged_trimmed.fasta")
    else:
        parser.print_usage()

if __name__ == "__main__":
    main()
