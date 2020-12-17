# MACOVID
Maastricht Covid processing pipeline

## Requirements:

### Python packages:

shutil
argparse
yaml
os
pandas
re
glob
sys

and Anaconda/Miniconda

## How to run

### Name changer

The content of the manifest.csv file:

| barcode_ID | sample_ID |
| ---------- |:---------:|
| barcode1   | SRR01     |
| barcode2   | SRR234    |
| barcode4   |           |
| barcode9   | SRR678_1  |
| barcode32  | SRR0234   |


### Core run

Requires input directory, manifest file (csv), output directory and the number of cores for running.
MACOVID scans for fastq files in the input directory. Found fastq files are automatically renamed from barcode to sample IDs based on information in the manifest. Fasta and concensus files are generated in the output directory.

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -m MACOVID_MANIFEST.csv -o OUTPUT_DIRECTORY --cores X
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

To rerun samples from specific folder. It is possible to run files in gz format. The output directory and number of cores must be specify.

Note: The program makes use of Snakemake so if the output directory contains the final files of the rerun, those files will not be analysis. Solution: remove old files or pick a new output directory. 

```
python MACOVID.py rerun -i FASTQ_DIRECTORY -o OUTPUT_DIRECTORY --cores X
```

