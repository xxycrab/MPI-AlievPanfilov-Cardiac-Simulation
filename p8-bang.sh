#!/bin/bash
# This script is intepreted by the Bourne Shell, sh
#
# Documentation for SGE is found in:
# http://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html
#
# Tell SGE which shell to run the job script in rather than depending
# on SGE to try and figure it out.
#$ -S /bin/bash
#
# Export all my environment variables to the job
#$ -V
# Tun the job in the same directory from which you submitted it
#$ -cwd
#
# Give a name to the job
#$ -N apf-strong-scale-8
#
# *** Specify the number of cores
#$ -pe orte 8

# Specify a time limit for the job
#$ -l h_rt=00:04:00
#
# Join stdout and stderr so they are reported in job output file
#$ -j y
#
# There are 4 queues
#
# 1. Debug queue (debug.q): only 1 node may be used at a time
#                           for up to 10 minutes
#
# 2. mpi queue (mpi.q): job may use up to 8 nodes for up to 10 minutes
#
# 3. debug-mpi.q:job may use up to 2 nodes for up to 20 minutes
#
# 4.  Interactive queue (qlogiN) batch jobs, up to 8 cores for 30 minutes
#     Maximum of 1 job per user running at a time

# Run on the debug queue
#$ -q debug.q
#
# Specifies the circumstances under which mail is to be sent to the job owner
# defined by -M option. For example, options "bea" cause mail to be sent at the 
# begining, end, and at abort time (if it happens) of the job.
# Option "n" means no mail will be sent.
#$ -m aeb
#
# *** Change to the address you want the notification sent to
#$ -M yourLogin@ucsd.edu
#

# Change to the directory where the job was submitted from
cd $SGE_O_WORKDIR

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

# Commands go here

mpirun -np 1 ./apf -n 400 -i 2000 -x 1 -y 1
mpirun -np 2 ./apf -n 400 -i 2000 -x 1 -y 2
mpirun -np 4 ./apf -n 400 -i 2000 -x 1 -y 4
mpirun -np 8 ./apf -n 400 -i 2000 -x 1 -y 8
mpirun -np 2 ./apf -n 400 -i 2000 -x 1 -y 2 -k
mpirun -np 4 ./apf -n 400 -i 2000 -x 1 -y 4 -k
mpirun -np 8 ./apf -n 400 -i 2000 -x 1 -y 8 -k
