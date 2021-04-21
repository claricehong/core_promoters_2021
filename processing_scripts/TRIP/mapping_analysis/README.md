These scripts are used for mapping TRIP integrations. 

Before starting, make sure there is a TRIP barcodes file that contains the REVERSE COMPLEMENT barcode because the barcodes are in the R2 read. The files needed to run these scripts are in this folder, along with sample fastq files. Full fastq files are provided in GEO.

1. Demultiplex the fastq file using Illumina_demultiplex_mapping_R1_only.py

```
usage: Illumina_demultiplex_mapping_R1_only.py [-h] [-o OUTPUT] mpBC fastq_R1 fastq_R2

Demultiplexes mapping fastq files

positional arguments:
  mpBC                  file containing multiplexing barcodes
  fastq_R1              input fastq file R1  
  fastq_R2              input fastq file R2
  optional arguments:  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT                        output file basename

  ```

  Example:

  ```
<<<<<<< HEAD
  python3 Illumina_demultiplex_mapping_R1_only.py MPBC_file.txt mapping_R1_sample.fastq.gz mapping_R2_sample.fastq.gz
=======
  python3 Illumina_demultiplex_mapping_R1_only.py MPBC_file.txt 
>>>>>>> parent of a026744... Updated README with examples and added sample files to run the commands.
  ```

2. Run bwa_align.sbatch

The demultiplexing file automatically generates a lookup file for input into bwa_align.sbatch for parallel processing. For other servers/non-parallel processing, run 

```
bwa mem <path to genome> <read1> > <output>.sam
```
The sbatch file also converts the sam file to a bam file for more efficient storage. 

3. Run annotate_integrations_R1_only.py using annotate_integrations.sbatch. 

Required python packages: pysam 

```
usage: annotate_integrations_R1_only.py [-h] bam insulator_bcs output

Read bam file and output barcode with annotated locations

positional arguments:
  bam            bamfile
  promoter_bcs   promoter BC file
  output         output file name

optional arguments:
  -h, --help     show this help message and exit
```

This generates a bed file that I usually intersect with the chromHMM file (K562_chromHMM_05032020.bed) using bedtools to get the chromatin anntations. However this can also be done in R if preferred. 

