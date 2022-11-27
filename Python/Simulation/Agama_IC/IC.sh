#!/bin/bash
#PBS -N agama_ic
#PBS -m abe
#PBS -M 2000cbh_@sjtu.edu.cn
#PBS -l nodes=1:ppn=1
#PBS -l mem=80gb
#PBS -l walltime=2:00:00
#PBS -q debug
#PBS -j oe
#PBS -o $HOME/Pipeline/Python/Simulation/Agama_IC/output.log

echo running with ${PBS_NP} processes on ${PBS_NUM_NODES} nodes
source $HOME/.cbh_set/cbh_set.sh
# run your own program!!!
cd $Agama
# Create IC:
python ic.py -id=$VMGIC/test_ic -od=$VMGIC/test_ic -if=diskhalo -of=model_ -st=Nbody -mt='two component' -ir=$INI

# Txt to hdf5:
python txt2hdf5.py -id=$VMGIC/test_ic -od=$VMGIC/test_ic -if=model_ -of=diskhalo -mt='two component' -ir=$INI
echo Done!
