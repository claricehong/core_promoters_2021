The script in this folder is used to get barcodes from DNA and RNA barcode sequencing. The fastq files are already demultiplexed by the sequencing core so I just have to get the barcodes out and count them. Sample files are provided in this folder. Full sequencing files can be found in GEO.

```
usage: bc2_counts_from_fastq_collapsed.py [-h] fastq primer TRIP_BCs output

Gets promoter/TRIP barcode counts from fastq files

positional arguments:
  fastq       input fastq file (either gzip or uncompressed)
  primer      primer sequence used for PCR  
  TRIP_BCs    promoter names and barcodes  
  output      output file basename
  optional arguments:
  -h, --help  show this help message and exit
```

Example:

```
python3 bc2_counts_from_fastq_collapsed.py expression_DNA_sample.fastq.gz CTCGCTTCGAGTCTAGA core_promoter_TRIP_bcs.txt <output_basename>
```