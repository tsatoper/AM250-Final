#!/bin/bash

#SBATCH -p 128x24
#SBATCH -J test
#SBATCH -e testc.err
#SBATCH -o testc.out
#SBATCH -N 1
#SBATCH -c 1234
#SBATCH -t 00:01:00 

mpirun -np 1234 column.exe
