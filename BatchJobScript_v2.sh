#!/bin/bash
#PBS -l nodes=1:ppn=20    
#PBS -l walltime=48:00:00  
#PBS -N Tyrannidae_subsp   
#PBS -o <file name>        
#PBS -e <file name>        
#PBS -q checkpt            
#PBS -A hpc_masonlab03     
#PBS -m e                  
#PBS -M maggie.macpherson@gmail.com 

module load r/4.0.3/intel-19.0.5

cd $PBS_O_WORKDIR

R --file=/home/maggie/HPC/2b_Migration_Trait_Evolution.R

