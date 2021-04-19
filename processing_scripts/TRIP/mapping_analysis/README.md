These scripts are used for mapping TRIP integrations. 

Before starting, make sure there is a TRIP barcodes file that contains the REVERSE COMPLEMENT barcode because the barcodes are in the R2 read. 

1. Demultiplex the fastq file using Illumina_demultiplex_mapping_R1_only.py

```
usage: Illumina_demultiplex_mapping_R1_only.py [-h] [-o OUTPUT]
                                               mpBC fastq_R1 fastq_R2

Demultiplexes mapping fastq files

positional arguments:
  mpBC                  file containing multiplexing barcodes
  fastq_R1              input fastq file R1  fastq_R2              input fastq file R2optional arguments:  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT                        output file basename

  ```

2. Run bwa_align.sbatch
3. Run annotate_integrations.py using annotate_integrations.sbatch. The sbatch file will intersect the locations with chromHMM, but that function can be disabled by simply commenting it out. 