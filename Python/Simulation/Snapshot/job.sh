#!/bin/bash
#PBS -N JobName
#PBS -m abe
#PBS -M 2000cbh@sjtu.edu.cn
#PBS -l nodes=1:ppn=72
#PBS -l mem=360gb
#PBS -l walltime=48:00:00
#PBS -q normal
#PBS -j oe
#PBS -o /home/bhchen/my_output.log
source /usr/share/Modules/init/bash
echo running with ${PBS_NP} processes on ${PBS_NUM_NODES} nodes
source /home/bhchen/.cbh_set/cbh_set.sh
# run your own program!!!
 
 
echo Done!
 
