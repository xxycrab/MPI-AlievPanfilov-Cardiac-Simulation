#!/bin/bash
#### this a bash script which we write out and then submit to the batch queue
# This script is intepreted by the Bourne Shell, sh
#
#SBATCH --account=csd562
#SBATCH --job-name=apf-strong-scale-48
#SBATCH --output="apf-strong-scale-48.out"
#SBATCH --partition="compute"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
####   #SBATCH --mail-type=BEGIN,END,FAIL
####   #SBATCH --mail-user=your_email@ucsd.edu
#SBATCH -t 00:03:00
# Commands go here

ibrun -np 48 ./apf -n 1800 -i 2000 -x 2 -y 24
ibrun -np 48 ./apf -n 1800 -i 2000 -x 4 -y 12
ibrun -np 48 ./apf -n 1800 -i 2000 -x 6 -y 8
ibrun -np 48 ./apf -n 1800 -i 2000 -x 8 -y 6
ibrun -np 48 ./apf -n 1800 -i 2000 -x 12 -y 4
ibrun -np 48 ./apf -n 1800 -i 2000 -x 24 -y 2


