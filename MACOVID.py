#!/usr/bin/env python3
import shutil
from argparse import ArgumentParser
import yaml
import os
import pandas as pd


locationrepo = os.path.dirname(os.path.abspath(__file__)) 

def get_absolute_path(path):
    return os.path.abspath(path)

def file_name_generator(filepath):
    return os.path.splitext(os.path.basename(filepath))[0]

###################
# Snakemake pipeline for generating consensus fastas
###################

def snakemake_in(samples, outdir):
    samplesdic = {}
    samplesdic['parameters'] = {}
    samplesdic['parameters']["outdir"] = get_absolute_path(outdir)
    samplesdic["SAMPLES"] = {}
    
    # generate the samples dictionary as input for snakemake 
    for i in samples:
        samplename = file_name_generator(i)
        samplesdic["SAMPLES"][samplename] = get_absolute_path(i)
    data = yaml.dump(samplesdic, default_flow_style=False)
    
    # make and write config file location
    os.system(f"mkdir -p {locationrepo}/config")
    with open(f"{locationrepo}/config/config.yaml", 'w') as f:
        f.write(data)

####################
# changing the names of barcoded samples
####################

def basenamechanger(namedict, path):
    """
    take full path a sample and change the base name of that sample accordig to the added dictionary
    """
    basename = path.split("/")[-1].split(".fastq")[0]
    location = "/".join(path.split("/")[:-1])
    

    return location + "/" + namedict[basename]

def change_names(samples, manifest, reverse ):
    """
    take samples and change name according to the manifest file.
    if -rev flag is used, the samples are converted back 
    """
    print(f'reverse is {reverse}')
    if reverse:
        pos = 1 
    else:
        pos = 0 
    names = pd.read_csv(manifest, index_col = pos, sep = ";").dropna().to_dict()
    # take the one and only key in this nested dict
    key = [x for x in names.keys()][0]
    names = names[key]
    for i in samples:
        try:
            newname = basenamechanger(names, i)
            shutil.move(i, newname + ".fastq")
        except:
            print(f'Could not change name of {i}, not present in manifest')
    



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
    mapreads.add_argument("-i", required = True, dest = "input_files", nargs = "+", help = 'give fastq files, and use basename to shoot into the pipeline')
    mapreads.add_argument("--cores", dest = 'cores', required = True, type = int, help = 'Number of CPU cores to use')
    mapreads.add_argument("-o", required = True, dest = "outdir")

    namechange = subparsers.add_parser("namechanger", help = "change barcode names")
    namechange.add_argument("-i", required = True, nargs = "+", dest = "input_samples")
    namechange.add_argument("-csv", required = True, dest = "manifest")
    namechange.add_argument("-rev", required = False, dest = "rev", action = "store_true") 

####################
# parsing part
####################

    args = parser.parse_args(command_line)
    if args.mode == "mapreads":
        snakemake_in(
                samples = args.input_files,
                outdir = args.outdir,
                )
        os.chdir(f"{locationrepo}")
        os.system(f"snakemake --cores {args.cores} --use-conda")

    elif args.mode == "namechanger":
        change_names(
                samples = args.input_samples,
                manifest = args.manifest,
                reverse = args.rev 
                )
    else:
        parser.print_usage()

if __name__ == "__main__":
    main()
