#!/bin/bash
#SBATCH --job-name=mc6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

# Load modules
module load mpi/openmpi-4.1.1 icc/2019.5 anaconda/2022a
source activate /home/gridsan/dkastner/.conda/envs/AmberTools22

# Execute analyze script
cpptraj -i cpptraj_cacovar.in