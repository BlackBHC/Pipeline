#!/bin/bash
#PBS -N see_ic
#PBS -m abe
#PBS -M 2000cbh_@sjtu.edu.cn
#PBS -l nodes=1:ppn=1
#PBS -l mem=80gb
#PBS -l walltime=2:00:00
#PBS -q debug
#PBS -j oe
#PBS -o $HOME/Pipeline/Python/Simulation/Snapshot/output.log

echo running with ${PBS_NP} processes on ${PBS_NUM_NODES} nodes
source /usr/share/Modules/init/bash
source $HOME/.cbh_set/shell.sh
cd $Agama
# Create IC:
python ic.py -id=$VMGIC/test_ic -od=$VMGIC/test_ic -if=Sandeeps -of=model_ -st=Nbody -mt='two component' -ir=$INI

# Txt to hdf5:
python txt2hdf5.py -id=$VMGIC/test_ic -od=$VMGIC/test_ic -if=model_ -of=diskhalo -mt='two component' -ir=$INI

cd $Snap
# For hdf5 to txt file:
python hdf52txt.py -id=$VMGIC/test_ic -od=$VMGIC/Fig -if=diskhalo -of=model -ir=$INI

# For radial profile
python radial_profile.py -id=$VMGIC/Fig -od=$VMGIC/Fig -if=model -of=radial_profile -ir=$INI
echo Done!
