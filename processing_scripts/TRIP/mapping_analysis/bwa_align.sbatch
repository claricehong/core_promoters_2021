#!/usr/bin/env bash

#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@htcf.wustl.edu
#SBATCH --mem=8196
#SBATCH --cpus-per-task=4

#SBATCH --array=1-8%8

module load bwa
module load samtools

read read1 read2 output < <( sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt )

bwa mem /home/hongkityeeclarice/cohen_lab/hg38_ensembl/hg38 ${read1} ${read2} > ${output}.sam

samtools view -bS ${output}.sam > ${output}.bam

samtools sort ${output}.bam | samtools index - ${output}.bam.bai

rm ${output}.sam
