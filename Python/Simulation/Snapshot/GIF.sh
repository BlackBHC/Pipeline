#!/bin/bash
#PBS -N mini_movie_1
#PBS -m abe
#PBS -M 2000cbh@sjtu.edu.cn
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=3:00:00
#PBS -q debug
#PBS -j oe
#PBS -o /home/bhchen/FeGradient/Fig/gif.txt
 
echo running with ${PBS_NP} processes on ${PBS_NUM_NODES} nodes
source /home/bhchen/configure/shell.sh
# run your own program!!!

cd $Snap

fmt=jpg

i=1000
for ((i=1000;i<=1400;i++))
do
	var=${i:1}
	indir=$VMG/Nbody/output 
	outdir=$VMGFig
	python cartesian_map.py -id=$indir -od=$outdir -if=snapshot_$var -of=cartesian_map_1 -ir=$INI
	mv $outdir/cartesian_map_1_column_density.${fmt} $outdir/cartesian_map_1_column_density$[$i-1000].${fmt}
	echo ${i:1}
done

# For gif images from some figs
#python movie.py -id=$outdir -od=$outdir -if=cartesian_map_1_column_density -of=cartesian_map_1 -ir=$INI

echo Done!
 
