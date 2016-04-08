#!/bin/bash

#SBATCH --account=k1123
#SBATCH --job-name=nektester_6144
#SBATCH --error=/scratch/project/k1123/GB16/tester/tester/nektester_6144.%j.err
#SBATCH --output=/scratch/project/k1123/GB16/tester/tester/nektester_6144.%j.out
#SBATCH --workdir=/scratch/project/k1123/GB16/tester/
#SBATCH --nodes=196608
#SBATCH --time=00:30:00

export MPICH_USE_DMAPP_COLL=1
srun --ntasks=6144 --ntasks-per-node=1 numactl -H
srun --ntasks=196608 --ntasks-per-node=32 --cpus-per-task=2 ./nektester
