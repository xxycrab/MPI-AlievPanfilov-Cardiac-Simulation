#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-nocomm-480
#SBATCH --output="apf-strong-scale-nocomm-480.out"
#SBATCH --partition="compute"
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=x1xie@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 480 ../apf -n 8000 -i 2000 -x 2 -y 240 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 3 -y 160 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 4 -y 120 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 6 -y 80 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 8 -y 60 -k
# ibrun -np 480 ../apf -n 8000 -i 2000 -x 10 -y 48 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 12 -y 40 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 15 -y 32 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 16 -y 30 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 20 -y 24 -k
ibrun -np 480 ../apf -n 8000 -i 2000 -x 24 -y 20 -k
