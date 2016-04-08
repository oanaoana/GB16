#!/bin/bash

#SBATCH --account=k1123
#SBATCH --job-name=nektester_3072
#SBATCH --error=/scratch/project/k1123/GB16/tester/tester/nektester_3072.%j.err
#SBATCH --output=/scratch/project/k1123/GB16/tester/tester/nektester_3072.%j.out
#SBATCH --workdir=/scratch/project/k1123/GB16/tester/
#SBATCH --nodes=3072
#SBATCH --time=00:30:00

export MPICH_USE_DMAPP_COLL=1
srun --ntasks=3072 --ntasks-per-node=1 numactl -H
srun --ntasks=98304 --ntasks-per-node=32 --cpus-per-task=2 ./nektester
