#!/bin/bash
cat *.out | awk 'f{print $4;f=0} /Segmentation/{f=1}' > ids_to_resubmit
filename="ids_to_resubmit"
# Read in as a variable the number of coordinates (always first line in a file, then add 2 for xyz format)
while read -r line; do
    index="$line"
    echo "Your failed job is $index, Let's resubmit that"
    for i in $(seq $index 1 $index); do
    echo "#!/bin/bash
#SBATCH --job-name=straggler
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=20

module load icc/2019.5
module load cuda/11.0
module load terachem/1.9-2021.10-dev

# Each job takes about 10 minutes, 12 GPUs allocated. Could run about 120 jobs per node per day, so 1440 jobs per user
# your command to run terachem
cd $i
terachem ../inputfiles/${i}.input > ${i}.out
echo "I just completed $i job"
cd ../
" > straggler_${i}.sh
done
sbatch straggler_*.sh
done < "$filename"
echo "I'm cleaning up after myself"
rm ids_to_resubmit
rm straggler_*.sh
