#!/bin/bash

#SBATCH --job-name=demo
#SBATCH --mem=16000M
#SBATCH -o parfum_%A.out # Standard output
#SBATCH -e parfum_%A.err # Standard error
#SBATCH --comment=Demo

module load parfums

srun /scratch/gdlab/manish/parfums/software/1.1.0/parfum_wrapper.pl --config 03_config_file.txt --step 0
