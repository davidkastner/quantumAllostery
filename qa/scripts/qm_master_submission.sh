#!/bin/bash

########################
# QM MASTER SUBMISSION #
########################

# TO USE THIS SCRIPT, RUN IT IN THE SAME DIRECTORY AS YOUR .XYZ TRAJECTORY
# 1. ACTIVATE --> chmod +x qm_master_submission.sh
# 2. RUN --> ./qm_master_submission trajectory.xyz frame_offset

# EXAMPLE RUN
# ./qm_master_submission all_coors.xyz 100

echo "A"
filename="$1" # Put your xyz file here, such that "./cut_coords.sh file.xyz stride"
echo "B"
stride="$2"
# written by CRR 12/19
# number of atoms is always the first line of the file, so use awk to grab that, and add 2 for the total length of one coordinate block
mkdir coordinates
echo "D"
rm -rf cut_xyz_commands.txt
echo "E"
# Read in as a variable the number of coordinates (always first line in a file, then add 2 for xyz format)
read xyz_size <<< $(awk 'NR==1 {print $1+2}' $filename)
echo "F"
echo "your xyz file is $xyz_size lines, +2 from the number of coordinates"
# get number of lines in the file in total with wc -l 
read file_size <<< $(wc -l $filename | awk '{print $1}')
echo "Your file is $file_size lines"
# Use BC to get the number of frames by dividing file line count by frame line count
read number_of_frames <<< $(bc <<< "$file_size/$xyz_size")
echo "You have $number_of_frames frames"
# We need an index for our loop to cut out the files. Note csplit can be a great alternative to this, but needs a specific string or delimiter
read index <<< $(bc <<< "$number_of_frames-1")
echo $index
for i in $(seq 0 $stride $index); do
# Getting indices for file cutting. File lines start at 1 not 0, hence add 1.
let j=$i+1
let k=$i*$xyz_size+1
let l=$j*$xyz_size
echo "sed -n '$k,${l}p' $filename > ./coordinates/$i.xyz" >> cut_xyz_commands.txt
# This currently makes xyz files in the working directory. Can be redirected by changing the line above
chmod +x cut_xyz_commands.txt
# Uncomment the last line to run and generate the files
#./cut_xyz_commands.txt
done
echo "I wrote the file to make the coordinates"
# Run after loop to avoid issues. I just moved all the files after they were generated
./cut_xyz_commands.txt

# INPUTFILE SECTION
mkdir inputfiles
for i in $(seq 0 $stride $index); do
mkdir $i
# Needed to read file size like in first script for index here. You can use echo like I do here, or you could just put a string in and used sed to replace it
# This is currently set up for TC on supercloud
echo "coordinates ../coordinates/${i}.xyz
basis lacvps_ecp
method wpbeh
charge 0
spinmult 1
guess ../opt-wfn/c0
maxit 1000
scf diis+a
scrdir ./scr
pcm cosmo
epsilon 80
pcm_radii read
pcm_radii_file /home/gridsan/dkastner/src/terachem/pcm_radii
end" > $i.input
done
mv *.input ./inputfiles/

echo "#!/bin/bash
#SBATCH --job-name=mc6s_8
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=20

module load icc/2019.5 cuda/11.0
#---TC setup---
export TeraChem=/data1/groups/HJKgroup/src/terachem/build
export terachem=\$TeraChem/bin/terachem
export PATH=\$TeraChem/bin:\$PATH
export LD_LIBRARY_PATH=\$TeraChem/lib:\$LD_LIBRARY_PATH

# Each job takes about 10 minutes, 12 GPUs allocated. Could run about 120 jobs per node per day, so 1440 jobs per user
# your command to run terachem
for i in \$(seq 0 ${stride} ${index}); do
cd \$i
terachem ../inputfiles/\${i}.input > \${i}.out
echo "I just completed \$i job"
sleep 10
cd ../
done" > job_submission.sh



