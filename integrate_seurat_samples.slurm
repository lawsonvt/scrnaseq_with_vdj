#!/bin/bash

#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -o seurat_integrate.out
#SBATCH -p standard
#SBATCH -A lukenslab
#SBATCH --mem=150GB

module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.5.0

Rscript src/integrate_seurat_samples.r
