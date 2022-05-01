#!/bin/bash
#SBATCH -p schumer,hns,normal,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 3:00:00
#SBATCH --mem 48000
~/miniconda3/envs/TTev/bin/kallisto quant -i xbirchmanni_10x_inferred_transcriptome.idx -o S9_TSR1-22C-malxbirchF1-brain_9-C4_kallisto ./rnaseq_reads/S9_TSR1-22C-malxbirchF1-brain_9-C4_R1_val_1.fq.gz ./rnaseq_reads/S9_TSR1-22C-malxbirchF1-brain_9-C4_R2_val_2.fq.gz
