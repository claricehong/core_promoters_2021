These scripts are used for mapping TRIP integrations. 

Before starting, make sure there is a TRIP barcodes file that contains the REVERSE COMPLEMENT barcode because the barcodes are in the R2 read. 

1. Run Illumina_demultiplex_mapping.py using Illumina_demultiplex_mapping.sbatch. 
2. Run bwa_align.sbatch
3. Run annotate_integrations.py using annotate_integrations.sbatch. The sbatch file will intersect the locations with chromHMM, but that function can be disabled by simply commenting it out. 