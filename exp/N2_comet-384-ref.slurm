#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-ref-384
#SBATCH --output="apf-strong-scale-ref-384.out"
#SBATCH --partition="compute"
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=x1xie@ucsd.edu
#SBATCH -t 00:04:30
# Commands go here

ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 2 -y 192
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 3 -y 128
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 4 -y 96
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 6 -y 64
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 8 -y 48
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 12 -y 32
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 16 -y 24
ibrun -np 384 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 24 -y 16
