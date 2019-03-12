#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-nocomm-72
#SBATCH --output="apf-strong-scale-nocomm-72.out"
#SBATCH --partition="compute"
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=x1xie@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 72 ../apf -n 1800 -i 2000 -x 2 -y 36 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 3 -y 24 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 4 -y 18 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 6 -y 12 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 8 -y 9 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 9 -y 8 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 12 -y 6 -k
ibrun -np 72 ../apf -n 1800 -i 2000 -x 18 -y 4 -k
