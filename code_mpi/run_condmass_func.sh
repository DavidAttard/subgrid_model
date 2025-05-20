#!/bin/bash
#SBATCH --account=lgreion
#SBATCH --job-name=testing
#SBATCH --output=log/%x.%j.o
#SBATCH --error=log/%x.%j.e
#SBATCH --cpus-per-task=1          # each task gets 1 CPU core
#SBATCH --nodes=1       # allocate 2 nodes
#SBATCH --ntasks-per-node=4     # up to 48 tasks per node
#SBATCH --time=03:00:00
#SBATCH --mem=0
#SBATCH --constraint=largedata
#SBATCH --array=00-00
module load Intel

# Print memory usage at the start
echo "Memory usage before job starts:"
free -h

# This will run every checkpoint as a seperate task
srun /p/project/lgreion/david/subgrid_MultiGrid/LMACH/code_mpi/read_slurm_info.sh

# Print memory usage after the job finishes
echo "Memory usage after job finishes:"
free -h