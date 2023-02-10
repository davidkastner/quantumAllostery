#!/bin/bash
#SBATCH --job-name=mc6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

# Load modules
module load mpi/openmpi-4.1.1 icc/2019.5 anaconda/2022a

# Execute analyze script
OMP_NUM_THREADS=20 python -u /home/gridsan/dkastner/src/quantumAllostery/qa/analyze.py > analyze.out
