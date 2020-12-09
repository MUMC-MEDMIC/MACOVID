# MACOVID
Maastricht Covid processing pipeline

## How to run

### Name changer

The content of the manifest.csv file:

| barcode_ID | sample_ID |
| ---------- |:---------:|
| barcode1   | SRR01     |
| barcode2   | SRR234    |
| barcode99  | SRR678_1  |

To rename sample ID to barcode:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -csv MACOVID_manifest.csv 
```

To reverse barcode back to sample ID:

```
python MACOVID.py namechanger -i FASTQ_DIRECTORY -csv MACOVID_manifest.csv -rev
```

To run the analysis:
```
python MACOVID.py mapreads -i ALLMYcovidDATA.fastq -o OUTPUTDIR --cores X
```
