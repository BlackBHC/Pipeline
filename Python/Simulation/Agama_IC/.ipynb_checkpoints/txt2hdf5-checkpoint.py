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
	usage       = "Specify the directory and the filename of the input text file and output hdf5 file."
	description = "Parameters transfer from the command line into the program."
	parser      = argparse.ArgumentParser(description = description, conflict_handler = 'resolve')
	parser.add_argument("-indir", '-id', type=str, default="/home/bhchen/BPX_Populations/IC/test/", help='Input directory where \
			has the txt file for different components.')
	parser.add_argument("-outdir", '-od', type=str, default="/home/bhchen/BPX_Populations/IC/test/", help='Output directory where \
			to output the IC hdf5 file.')
	parser.add_argument("-infile", '-if', type=str, default="model_", help='Prefix of filename for different txt file\
			of different components')
	parser.add_argument("-outfile", '-of', type=str, default="model", help='Prefix of the output hdf5 file without sufix\
			.hdf5')
	parser.add_argument("-iniroot", '-ir', type=str, default="/home/bhchen/Codes/Simulation/INI/", help='Root directory of the\
			paramter INI files')
	parser.add_argument("-simulationtype", '-st', type=str, default="Nbody", help='Simulation type: Nbody or SPH et al.')
	parser.add_argument("-modeltype", '-mt', type=str, default="two component", help='Model type')
	return parser.parse_args()

# Command line arguments parsing
args            = parser_args()
in_dir          = args.indir
out_dir         = args.outdir
in_file         = args.infile
out_file        = args.outfile
ini_root        = args.iniroot
type_simulation = args.simulationtype
type_model      = args.modeltype


# Read in Model.ini configuration file
config             = cp.ConfigParser(allow_no_value=True)
config.read(ini_root+'Model.ini')
model_type         = config[type_simulation][type_model]
components         = model_type.split()
num_iteration      = config['Iteration'].getint('iterations')
num_particle_comps = []


components = ('disk', 'halo')
raw = []
for comp in components:
    raw += [np.loadtxt(indir+"model_"+comp+"_final")]
print("Read in the txt file of different components...")


f = h5py.File(outdir+outfile+'.hdf5', 'w')
f.create_group("Header")
print("Create the hdf5 file '{}' in {}...".format(outfile+'.hdf5', outdir))

f["Header"].attrs['NumPart_ThisFile'] = np.array([0, len(raw[0]), len(raw[1])], dtype=np.uint32)
f["Header"].attrs['NumPart_Total'] = np.array([0, len(raw[0]), len(raw[1])], dtype=np.uint64)
#f["Header"].attrs['MassTable'] = np.array([.1, np.mean(np.sum(raw[0][:, -1]))/len(raw[0])\
#		, np.mean(np.sum(raw[1][:, -1]))/len(raw[1]), np.mean(np.sum(raw[2][:, -1]))/len(raw[2])])
f["Header"].attrs['MassTable'] = np.array([0, 0, 0])
f["Header"].attrs['Time'] = 0.0
f["Header"].attrs['Redshift'] = 0.0
f["Header"].attrs['BoxSize'] = 0.0
f["Header"].attrs["NumFilesPerSnapshot"] = 0
print("Create the Header...")

cmps = []
cmps += [f.create_group("PartType1")]
cmps += [f.create_group("PartType2")]
cmps += [f.create_group("PartType0")]

i = 0 
count = 0
while(i < 2):
    cmps[i].create_dataset(name="Coordinates", data = raw[i][:,:3])
    cmps[i].create_dataset(name="Velocities", data = raw[i][:,3:6])
    cmps[i].create_dataset(name="Masses", data = raw[i][:,6])
    cmps[i].create_dataset(name="ParticleIDs", dtype=np.uint32, data = np.arange(count, count + len(raw[i])))
    count += len(raw[i])
    i += 1
print("Create the different PartType...")


f.close()
print("Finish the transformation!")

