#!/bin/bash

#SBATCH -J dmontiel_job
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -A prisms_project1
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mail-user=dmontiel@umich.edu
#SBATCH --mail-type=all

# Loading env3 
module purge

# Set environment variables
export DEAL_II_DIR=$HOME/dealii_9p5p2

# Load required modules
module load intel openmpi cmake tbb

# Set compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpifort

# Print message to confirm activation
echo "Switched to Environment 3 (Intel + Deal.II 9.5.2)"

module list
pwd
date

#Creating target directory on scratch and copying files into it
locdir=${PWD##*/}
targpath=/scratch/prisms_project_root/prisms_project1/dmontiel/env3
fulltargpath=$targpath/$locdir

echo $fulltargpath

mkdir $fulltargpath
cp parameters.prm $fulltargpath
cp main $fulltargpath

cd $fulltargpath
srun ./main -i parameters.prm
