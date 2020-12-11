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

To run the analysis:

```
python MACOVID.py mapreads -i FASTQ_DIRECTORY -csv MACOVID_MANIFEST -o OUTPUTDIR --cores X
```

MACOVID will scan for fastq files in the input directory. The script will only run if the directory contain the correct fastq files.
Found barcode files are renamed to sample IDs based on information from the manifest. Analysis is carried out right after and final files are stored in the output directory.


### Standalone renaming

Rename barcode to sample ID:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -csv MACOVID_manifest.csv 
```

To reverse sample ID back to barcode:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -csv MACOVID_manifest.csv -rev
```



