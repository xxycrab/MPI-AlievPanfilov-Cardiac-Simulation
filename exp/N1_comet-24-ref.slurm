#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-ref-24
#SBATCH --output="apf-strong-scale-ref-24.out"
#SBATCH --partition="compute"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=x1xie@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 24 $PUB/HW/pa3/apf-ref -n 1800 -i 2000 -x 2 -y 12
ibrun -np 24 $PUB/HW/pa3/apf-ref -n 1800 -i 2000 -x 3 -y 8
ibrun -np 24 $PUB/HW/pa3/apf-ref -n 1800 -i 2000 -x 4 -y 6
ibrun -np 24 $PUB/HW/pa3/apf-ref -n 1800 -i 2000 -x 6 -y 4
ibrun -np 24 $PUB/HW/pa3/apf-ref -n 1800 -i 2000 -x 8 -y 3
