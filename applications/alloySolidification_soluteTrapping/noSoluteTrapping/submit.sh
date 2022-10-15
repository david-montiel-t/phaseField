#!/bin/bash -l
# created: Jan 24, 2019 10:37 AM
# author: tpinomaa
#SBATCH -J          PFdir
#SBATCH -o       log
#SBATCH -e error_log
#SBATCH --partition test
#SBATCH -n 16
#SBATCH --time 0-00:15:00
#SBATCH --mem-per-cpu=2000
#SBATCH --account=project_2003810
##SBATCH --mail-type=END
##SBATCH --mail-user=tatu.pinomaa@vtt.fi


# partitions: https://docs.csc.fi/computing/running/batch-job-partitions/
# test: time limit 15 min, maximum 2 nodes, 80 tasks
# small: time limit 3 d, maximum 1 nodes, 40 tasks
# large: time limit 3 d, maximum 100 nodes, 4000 tasks

# For more information
#   man sbatch
#   more examples in Taito guide in
#   http://research.csc.fi/taito-user-guide

# example run commands

#-------------------------------------------------
# module load cmake
#root=/scratch/project_2003810/tpinomaa/ExaCA_gcc
#input_file=Inp_SingleScanTrack.txt
executable=../main

srun $executable $input_file 
#srun amitex_fftp  -nz $MATEVTK  -a $ALGOXML -c $LOADXML -m $MATEXML  -s ${results_folder}/res
# mpiexec -n 1 ../build/install/bin/ExaCA-Kokkos ../examples/Inp_DirSolidification.txt

