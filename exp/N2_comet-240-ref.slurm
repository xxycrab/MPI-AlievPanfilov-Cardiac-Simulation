#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-ref-240
#SBATCH --output="apf-strong-scale-ref-240.out"
#SBATCH --partition="compute"
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=x1xie@ucsd.edu
#SBATCH -t 00:04:30
# Commands go here

ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 2 -y 120
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 3 -y 80
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 4 -y 60
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 6 -y 40
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 8 -y 30
# ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 10 -y 24
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 12 -y 20
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 15 -y 16
ibrun -np 240 $PUB/HW/pa3/apf-ref -n 8000 -i 2000 -x 16 -y 15
