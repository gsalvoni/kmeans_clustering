#!/bin/bash
# Submission script for NIC4
#SBATCH --job-name=hw2hpc
#SBATCH --time=00:10:00 # hh:mm:ss
#
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=2625 # megabytes
#SBATCH --partition=defq
#
#SBATCH --mail-user=g.salvoni@student.ulg.ac.be
#SBATCH --mail-type=ALL

module load openmpi/1.6.4/gcc-4.8.1
mpirun ./kmeans_book.exe 4
