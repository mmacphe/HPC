#!/bin/bash
#PBS -l nodes=1:ppn10      #Number of nodes and processors per node
#PBS -l walltime=48:00:00  #Maximum wall time
#PBS -N Tyrannidae_subsp   #Job name; default is name of script
#PBS -o <file name>        #File name for standard output; default system will name
#PBS -e <file name>        #File name for standard error; default system will name
#PBS -q single             #Queue name for serial job; most common is checkpt or workq
#PBS -A hpc_masonlab03     #Allocation name
#PBS -m e                  #Send mail when job ends
#PBS -M maggie.macpherson@gmail.com #Send mail to this address

module load r/4.0.3/intel-19.0.5

cd $PBS_O_WORKDIR

R --file=/home/maggie/HPC/2b_Migration_Trait_Evolution.R
