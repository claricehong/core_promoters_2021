#!/usr/bin/env bash

#SBATCH -n 1
#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@htcf.wustl.edu
#sleep 1 # pause for 1 second between job submissions

##Change number of arrays to number of files you need to run
#SBATCH --array=1-12%12

eval "$(conda shell.bash hook)"
conda activate free_barcodes

R1=$( sed -n ${SLURM_ARRAY_TASK_ID}p R1_lookup.txt )
R2=$( sed -n ${SLURM_ARRAY_TASK_ID}p R2_lookup.txt )
basename=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )

python3 /home/hongkityeeclarice/cohen_lab/promoter_screen/freebarcodes/fastq_processing/generate_freebarcodes_input.py ${R1} ${R2} LP_BCs -o ${basename}

freebarcodes decode /home/hongkityeeclarice/cohen_lab/promoter_screen/freebarcodes/fastq_processing/barcodes12-1.txt ${basename}.fastq 

python3 /home/hongkityeeclarice/cohen_lab/promoter_screen/freebarcodes/fastq_processing/process_freebarcodes_output.py ${basename}_decoded.txt 
