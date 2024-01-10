#!/bin/bash 

#SBATCH -A m2360
#SBATCH -J myjob
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes 1
#SBATCH --tasks-per-node=64
#SBATCH -t 12:00:00
#SBATCH --mail-user=dmontiel@umich.edu
#SBATCH --mail-type=all

module list
pwd
date

#Creating target directory on scratch and copying files into it
locdir=${PWD##*/}
targpath=$SCRATCH
fulltargpath=$targpath/$locdir

echo $fulltargpath

mkdir $fulltargpath
cp parameters.prm $fulltargpath
cp main $fulltargpath
cp microstructure_US.vtk $fulltargpath
cd $fulltargpath

srun -n 64 --cpu_bind=cores ./main
