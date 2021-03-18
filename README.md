# MACOVID
Maastricht MUMC+ Covid pipeline, adapted from [ENA_SARS_Cov2_nanopore](https://github.com/dnieuw/ENA_SARS_Cov2_nanopore) and [Artic-ncovid2019](https://github.com/artic-network/artic-ncov2019).


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
conda env create -f mainEnvs.yaml
conda activate envsMacovid
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
 --c Number of cores use  
 --cov Set coverage (Default 30)  
 -l Run locally  


MACOVID scans for fastq files in the input directory. Please avoid space and symbols in the input folder. Found fastq files are automatically renamed based on information from the manifest file. The output directory must be specified for the results. All the final fasta files are concatenated into a single file using the name of the input folder. For locally run use -l command (optional).

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -m MACOVID_MANIFEST.csv -o OUTPUT_DIRECTORY --cores X -l --cov
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
 --c Number of cores use  
 --cov Set coverage (Default 30)  
 -l Run locally  


To rerun samples from specific folder. Input files could be in gz format. The output directory and number of cores must be specify. Local run use -l command (optional).

Note: The program makes use of Snakemake so if the output directory contains the final files of the rerun, those files will not be analysis. Solution: remove old files or pick a new output directory. 

```
python MACOVID.py rerun -i FASTQ_DIRECTORY -o OUTPUT_DIRECTORY --cores X -l -cov 30
```
