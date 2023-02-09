#!/bin/bash
#SBATCH --job-name=mc6s
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=20

module load icc/2019.5 cuda/11.0
#---TC setup---
export TeraChem=/data1/groups/HJKgroup/src/terachem/build
export terachem=$TeraChem/bin/terachem
export PATH=$TeraChem/bin:$PATH
export LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH

# Each job takes about 10 minutes, 12 GPUs allocated. Could run about 120 jobs per node per day, so 1440 jobs per user
# your command to run terachem
for i in $(seq 0 100 39999); do
cd $i
terachem ../inputfiles/${i}.input > ${i}.out
echo I just completed $i job
sleep 10
cd ../
done
