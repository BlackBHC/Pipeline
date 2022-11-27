import numpy as np
from matplotlib import pyplot as plt
import h5py
from scipy.stats import binned_statistic as bin1d
from scipy.stats import binned_statistic_2d as bin2d
import imageio
import sys
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
import random
import argparse


# Parse the parameters.
def parser_args():
	usage       = "Specify the directory and the filename of the output hdf5 file."
	description = "Parameters transfer from the command line into the program."
	parser      = argparse.ArgumentParser(description = description, conflict_handler = 'resolve')
	parser.add_argument("-indir", "-in", type=str, default="/home/bhchen/BPX_Populations/IC/", help='Input directory.')
	parser.add_argument("-outdir", "-out", type=str, default="./", help='Output directory.')
	parser.add_argument("-infile", "-if", type=str, default="snapshot", help=\
			'Name of the input hdf5 file without suffix')
	parser.add_argument("-outfile", "-of", type=str, default="model", help='Name of the output txt file')
	parser.add_argument("-recenter", "-r", type=np.uint32, default=1, help='Whether recenter the particles: 1=Yes, 0=Not.')
	return parser.parse_args()
args     = parser_args()
indir    = args.indir
outdir   = args.outdir
infile   = args.infile
outfile  = args.outfile
recenter = args.recenter

if recenter==1:
	print("\033[1;33m**** WILL RECENTER THE SYSTEM FOR EACH COMPONENT!\033[0m")
	print("\033[1;33m**** NOTE THE CENTER OF EACH COMPOENT MAY NOT COINCIDE WITH EACH OTHER!\033[0m")

components = ('disk', 'bulge', 'halo')
print("\033[1;33m**** {}.\033[0m".format("Read in the "+infile+'.hdf5'+" file"))
f = h5py.File(indir+infile+'.hdf5', 'r')
# the time of the snapshot
time = f['Header'].attrs['Time']
# number of PartType
num_comps = 0
comp_id   = []
for i in f.keys():
	if i[:8]=='PartType':
		if i[8]!='0':
			num_comps += 1
			comp_id   += [i[8:]]
	else: pass

comps = []
for i in range(num_comps): comps += [f['PartType'+comp_id[i]]]


# Recenter the components.
for n in range(num_comps):
	coord = comps[n]['Coordinates'][...]
	if recenter==1:
		coord -= np.mean(coord, axis=0)


for i in range(3):
	print("\033[1;33m**** Start transform the {} particles.\033[0m".format(components[i]))
	comp  = f["PartType{}".format(i+1)]
	vels  = comp['Velocities'][...]
	coord = comp['Coordinates'][...]
	mass  = comp['Masses'][...]

	data = np.column_stack((coord, vels))
	data = np.column_stack((data,  mass))
	np.savetxt(outdir+outfile+"_"+components[i], data)
	
print("\033[1;33m**** FINISHED!\033[0m")
f.close()

