# MACOVID
Maastricht MUMC+ Covid pipeline, adapted from [Artic-ncovid2019](https://github.com/artic-network/artic-ncov2019) pipeline.


## Setting up

#### Prerequisites

Required the following to be installed

```
- Anaconda/Miniconda
- git
```

#### Step 1: Obtain a copy of this workflow

[Clone](https://github.com/MUMC-MEDMIC/MACOVID.git) the Macovid into your local system, where you want to do the analysis.

```
git clone https://github.com/MUMC-MEDMIC/MACOVID.git
```

#### Step 2: Configure workflow

Configure the workflow by creating and activating the environment.

```
conda env create -f macovid.yaml
conda activate macovid
```

## Running MACOVID

### Manifest input

The manifest file could either be tab, comma or semicolon seperated. This is how the table should be:

| barcode_ID | sample_ID |
| ---------- |:---------:|
| barcode01  | SRR01     |
| barcode02  | SRR234    |
| barcode04  |           |
| barcode09  | SRR678_1  |
| barcode32  | SRR0234   |


### Analysis using manifest file

    -i Input directory  
    -m Manifest file  
    -o Output directory  
    --cores Number of threads   
    --cov min coverage per base (default: 30)  
    --cluster Run on a slurm cluster  
    --trim_start trim start of consensus that is not sequenced (default: 54)  
    --trim_end trim end of consensus that is not sequenced (default: 79)  
    --scheme path to scheme directory (default: primer_schemes/EMC)
    --scheme_prefix prefix of primer scheme (default: nCoV-2019)

MACOVID scans for fastq files in the input directory. Please avoid space and symbols in the input folder. Found fastq files are automatically renamed based on information from the manifest file. The output directory must be specified for the results. All the final fasta files are concatenated into a single file using the name of the input folder.

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -m MACOVID_MANIFEST.csv -o OUTPUT_DIRECTORY --cores X 
```

### Files renaming

    -i Input directory  
    -rev Reverse name changed

To rename from barcode to the sample id based on the manifest file:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -m MACOVID_manifest.csv 
```

To reverse sample ID back to barcode based on the manifest file:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -m MACOVID_manifest.csv -rev
```

### Rerun samples
 
    -i Input directory  
    -o Output directory  
    --cores Number of threads 
    --cov Set coverage (Default 30)  
    --cluster Run on a slurm cluster  
    --trim_start trim start of consensus that is not sequenced (default: 54)  
    --trim_end trim end of consensus that is not sequenced (default: 79)
    --scheme path to scheme directory (default: primer_schemes/EMC)
    --scheme_prefix prefix of primer scheme (default: nCoV-2019)

To rerun samples from specific folder. Input files could be in gz format. The output directory and number of cores must be specified. 


```
python MACOVID.py rerun -i FASTQ_DIRECTORY -o OUTPUT_DIRECTORY --cores X -l -cov 30
```

Note: The program makes use of Snakemake so if the output directory contains the final files of the rerun, those files will not be reanalysed. Solution: remove or move the old files or input a new output directory. 


## Add other primer schemes

To add other schemes one could add an additional folder in the primer_schemes/ directory.
The required files are *.primer.bed, *.reference.fasta, *.reference.fasta.fai, and *.scheme.bed.
These files can be obtained using primalscheme https://github.com/aresti/primalscheme
The path to the primer scheme and prefix should be defined in the pipeline:

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -m MACOVID_MANIFEST.csv -o OUTPUT_DIRECTORY --cores X --scheme path primer_schemes/YOUR_SCHEME --scheme_prefix YOUR_SCHEME_PREFIX
```	

The merged consensus file can be trimmed to remove not sequenced parts according to your scheme:

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -m MACOVID_MANIFEST.csv -o OUTPUT_DIRECTORY --cores X --trim_start XX --trim_end XX
```	