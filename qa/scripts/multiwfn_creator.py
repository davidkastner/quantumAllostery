import os

# Resource allocations
cpus = 4
nodes = 1
task = 1

# Lists to iterate
replicates = [1,2,3,4,5,6,7,8]
job_dirs = [str(dir) for dir in range(0, 39900, 100)]

# Divide the lists into pices
step = 100 # How far apart each qm job was
divide = 20 # Submit as junks of five jobs
final_value = job_dirs[-1]
sliced = job_dirs[::divide] 
sliced.append(final_value)
# Turn the list into tuples of form (first, last)
first_last_tuples = [(sliced[i], sliced[i+1]) for i in range(len(sliced)-1)]

# Generate the sbatch files and execute them
for replicate in replicates:
    for first_last in first_last_tuples:
        first = first_last[0]
        last = first_last[1]

        sbatch_file_name = f"j{replicate}_{last}.sh"
        with open(sbatch_file_name, "w") as job:
            lines = ["#!/bin/bash\n",
                    f"#SBATCH --job-name=j{replicate}_{last}\n",
                    f"#SBATCH --nodes={nodes}\n",
                    f"#SBATCH --cpus-per-task={task}\n",
                    f"#SBATCH --ntasks={cpus}\n",
                    f"python /home/gridsan/dkastner/src/quantumAllostery/cli.py -q {replicate} {first} {last} {step}\n"
                    ]
            job.writelines(lines)

        os.system(f"sbatch {sbatch_file_name}")
