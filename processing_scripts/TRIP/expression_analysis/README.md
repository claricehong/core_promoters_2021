The script in this folder is used to get barcodes from DNA and RNA barcode sequencing. The fastq files are already demultiplexed by the sequencing core so I just have to get the barcodes out and count them. 

```
usage: bc2_counts_from_fastq_collapsed.py [-h] fastq primer TRIP_BCs output

Gets insulator/TRIP barcode counts from fastq files

positional arguments:
  fastq       input fastq file (either gzip or uncompressed)
  primer      primer sequence used for PCR  
  TRIP_BCs    insulator names and barcodes  
  output      output file basenameoptional arguments:
  -h, --help  show this help message and exit
```

example usage:

```
python3 bc2_counts_from_fastq_collapsed.py <fastq_file> CTCGCTTCGAGTCTAGA core_promoter_TRIP_bcs.txt <output_basename>
```