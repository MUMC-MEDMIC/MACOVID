# MACOVID
Covid processing pipeline (Paragraph of explaination)


## Getting Started

To get the pipeline up and running you need to following.


### Prerequisites

```
Anaconda/Miniconda:
Git
```

### Installing

```
git clone https://github.com/MUMC-MEDMIC/MACOVID.git
```

## Create and activate the main environment

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


### Core run


**-i** Input directory
**-m** Manifest file
**-o** Output directory
**--c** Number of cores use
**-cov** Set coverage (Default 30)
**-l** Run locally


Requires input directory, manifest file (csv), output directory and the number of cores for running.
MACOVID scans for fastq files in the input directory. Found fastq files are automatically renamed from barcode to sample IDs based on information in the manifest. Fasta and concensus files are generated in the output directory. Local run use -l command (optional).

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -m MACOVID_MANIFEST.csv -o OUTPUT_DIRECTORY --cores X -l
```

### Standalone renaming

To rename from barcode to the sample id based on the manifest file:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -m MACOVID_manifest.csv 
```

To reverse sample ID back to barcode based on the manifest file:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -m MACOVID_manifest.csv -rev
```

To rerun samples from specific folder. It is possible to run files in gz format. The output directory and number of cores must be specify. Local run use -l command (optional).

Note: The program makes use of Snakemake so if the output directory contains the final files of the rerun, those files will not be analysis. Solution: remove old files or pick a new output directory. 

```
python MACOVID.py rerun -i FASTQ_DIRECTORY -o OUTPUT_DIRECTORY --cores X -l
```

Adapted from: https://github.com/dnieuw/ENA_SARS_Cov2_nanopore
