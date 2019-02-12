#!/bin/bash

#SBATCH -A csd562
#SBATCH --job-name="pa3"
#SBATCH --output="pa3.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH -t 01:30:00

#This job runs with 2  nodes, 24 cores per node for a total of 48 cores.

echo
echo " *** Current working directory"
pwd
echo
echo " *** Compiler"
# Output which  compiler are we using and the environment
mpicc -v 
echo
echo " *** Environment"
printenv

echo

echo ">>> Job Starts"
date

mpirun -np 1 ./apf -n 400 -i 2000 -x 1 -y 1
mpirun -np 2 ./apf -n 400 -i 2000 -x 1 -y 2
mpirun -np 4 ./apf -n 400 -i 2000 -x 1 -y 4
mpirun -np 8 ./apf -n 400 -i 2000 -x 1 -y 8
mpirun -np 2 ./apf -n 400 -i 2000 -x 1 -y 2 -k
mpirun -np 4 ./apf -n 400 -i 2000 -x 1 -y 4 -k
mpirun -np 8 ./apf -n 400 -i 2000 -x 1 -y 8 -k

date
echo ">>> Job Ends"
