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
	parser.add_argument("-infilename", "-if", type=str, default="snapshot", help='Name of the output hdf5 file without suffix')
	parser.add_argument("-outfilename", "-of", type=str, default="model", help='Name of the output hdf5 file without suffix')
	return parser.parse_args()
args      = parser_args()
indir     = args.indir
outdir    = args.outdir
infile    = args.infilename
outfile   = args.outfilename


components = ('disk', 'bulge', 'halo')
raw = []
for comp in components:
    raw += [np.loadtxt(directory+"model_"+comp+"_final")]
print("Read in the txt file of different components...")


f = h5py.File(directory+filename+'.hdf5', 'w')
f.create_group("Header")
print("Create the hdf5 file '{}' in {}...".format(filename+'.hdf5', directory))

f["Header"].attrs['NumPart_ThisFile'] = np.array([0, len(raw[0]), len(raw[1]), len(raw[2])], dtype=np.uint32)
f["Header"].attrs['NumPart_Total'] = np.array([0, len(raw[0]), len(raw[1]), len(raw[2])], dtype=np.uint64)
f["Header"].attrs['MassTable'] = np.array([0, 0, 0, 0])
f["Header"].attrs['Time'] = 0.0
f["Header"].attrs['Redshift'] = 0.0
f["Header"].attrs['BoxSize'] = 0.0
f["Header"].attrs["NumFilesPerSnapshot"] = 0
print("Create the Header...")

cmps = []
cmps += [f.create_group("PartType1")]
cmps += [f.create_group("PartType2")]
cmps += [f.create_group("PartType3")]
cmps += [f.create_group("PartType0")]

i = 0 
count = 0
while(i < 3):
    cmps[i].create_dataset(name="Coordinates", data = raw[i][:,:3])
    cmps[i].create_dataset(name="Velocities", data = raw[i][:,3:6])
    cmps[i].create_dataset(name="Masses", data = raw[i][:,6])
    cmps[i].create_dataset(name="ParticleIDs", dtype=np.uint32, data = np.arange(count, count + len(raw[i])))
    count += len(raw[i])
    i += 1
print("Create the different PartType...")


f.close()
print("Finish the transformation!")

