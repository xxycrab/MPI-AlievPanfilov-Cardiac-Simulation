#!/bin/bash

#SBATCH -A csd562
#SBATCH --job-name="pa3-serial"
#SBATCH --output="pa3-serial.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=10
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
./apf -n 150 -i 400
date
echo ">>> Job Ends"
