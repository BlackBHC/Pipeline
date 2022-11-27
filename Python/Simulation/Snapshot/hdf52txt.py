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
from matplotlib import font_manager, rcParams
config = {
    "font.family":"DejaVu Serif",
	"font.size": 22,
	"mathtext.fontset": "stix"
}
rcParams.update(config)
import configparser as cp



# Parse the parameters.
def parser_args():
	usage       = "Specify directory and filename of related parameter files."
	description = "Parameters transfer from the command line into the program."
	parser      = argparse.ArgumentParser(description = description, conflict_handler = "resolve")
	parser.add_argument("-indir", "-id", type=str, default="/home/bhchen/BPX_Populations/Nbody/output/", help="Input directory where \
			has the Config.ini file.")
	parser.add_argument("-outdir", "-od", type=str, default="/home/bhchen/Fig/", help="Output directory.")
	parser.add_argument("-infile", "-if", type=str, default="snapshot_400", help="File name of without sufix.")
	parser.add_argument("-outfile", "-of", type=str, default="cartesian_map", help="Prefix of the output file.")
	parser.add_argument("-iniroot", "-ir", type=str, default="/home/bhchen/Codes/Simulation/INI/", help="Root directory of the\
			paramter INI files")
	parser.add_argument("-simulationtype", "-st", type=str, default="Nbody", help="Simulation type: Nbody or SPH et al.")
	parser.add_argument("-modeltype", "-mt", type=str, default="two component", help="Model type")
	return parser.parse_args()



# Command line arguments parsing
args            = parser_args()
in_dir          = args.indir + '/'
out_dir         = args.outdir + '/'
in_file         = args.infile
out_file        = args.outfile
ini_root        = args.iniroot + '/'
type_simulation = args.simulationtype
type_model      = args.modeltype



# Read in configuration files
config             = cp.ConfigParser(allow_no_value=True)
config.read(ini_root+"Model.ini")
config.read(ini_root+"Figure.ini")
config.read(ini_root+"Statistic.ini")
print(config.keys())
model_type         = config[type_simulation][type_model]
components         = model_type.split()
# Whether has gas in this model?
gas_in_model         = False
for comp in components:
	if comp.lower() == "gas":
		gas_in_model = True
		break
# Whether needs to recenter the system?
recenter = config["Recenter"].getboolean("recenter")
if recenter:
	recenter_size       = config["Recenter"].getint("recenter size")
	recenter_diff_comps = config["Recenter"].getboolean("different components")
	# If need recenter, whether needs to recenter different components w.r.t. their CoM?
	# If not, recenter the system w.r.t. the CoM of all particles.
	if recenter_diff_comps:
		print("\033[1;33m**** WILL RECENTER THE SYSTEM FOR EACH COMPONENT!\033[0m")
	else:
		print("\033[1;33m**** WILL RECENTER THE SYSTEM AS A WHOLE!\033[0m")



# Read in the snapshot
print("\033[1;33m**** READING IN THE SNAPSHOT {}.hdf5\033[0m".format(in_file))
snapshot = h5py.File(in_dir+in_file+".hdf5", "r")
# the time of the snapshot
time     = snapshot["Header"].attrs["Time"]

comp_ids = []
if gas_in_model==False:
	for i in range(len(components)):
		comp_ids += [i+1]
else:
	for i in range(len(components)):
		comp_ids += [i]
comp_ids = tuple(comp_ids)
# PartType{id} list for all components



print("\033[1;33m**** READING IN DIFFERENT COMPONENTS IN THE SNAPSHOT\033[0m")
comps = []
for comp_id in comp_ids:
	comps += [snapshot["PartType{}".format(comp_id)]]
# store the data of different components



if recenter and not(recenter_diff_comps):
	denominator = 0
	numerator   = np.zeros(3)
	for i in range(len(comps)):
		masses       = comps[i]["Masses"][...]
		coords       = comps[i]["Coordinates"][...]
		index        = np.where(np.linalg.norm(coords, axis=1, ord=2)<recenter_size)[0]
		masses       = masses[index]
		coords       = coords[index]
		denominator += np.sum(masses)
		numerator   += np.sum(np.column_stack((masses, masses, masses))*coords, axis=0)
	CoM = numerator/denominator
	print("Center of mass for all particles:\n", CoM)
	# The center of mass for all particles
elif recenter:
	CoMs = []
	for i in range(len(comps)):
		masses      = comps[i]["Masses"][...]
		coords      = comps[i]["Coordinates"][...]
		index        = np.where(np.linalg.norm(coords, axis=1, ord=2)<recenter_size)[0]
		masses      = masses[index]
		coords      = coords[index]
		denominator = np.sum(masses)
		numerator   = np.sum(np.column_stack((masses, masses, masses))*coords, axis=0)
		com         = numerator/denominator
		CoMs       += [com]
	print("CoMs for each component in the model:\n", CoMs)
	# The center of mass for each component



for comp_id in range(len(comp_ids)):
	print("\033[1;33m**** Tansforming the {} component ...\033[0m".format(components[comp_id]))
	if recenter and not(recenter_diff_comps):
		n   = len(comps[comp_id]["Coordinates"][...])
		com = np.column_stack((np.ones(n)*CoM[0], np.ones(n)*CoM[1], np.ones(n)*CoM[2]))
		cartesian_coordinates = comps[comp_id]["Coordinates"][...] - com
	elif recenter:
		n   = len(comps[comp_id]["Coordinates"][...])
		com = np.column_stack((np.ones(n)*CoMs[comp_id][0], np.ones(n)*CoMs[comp_id][1], np.ones(n)*CoMs[comp_id][2]))
		cartesian_coordinates = comps[comp_id]["Coordinates"][...] - com
	else:
		cartesian_coordinates = comps[comp_id]["Coordinates"][...]	
	cartesian_velocities = comps[comp_id]["Velocities"][...] 
	masses               = comps[comp_id]["Masses"][...]
	roundn				 = lambda x, n=3: int(x*10**n)/10**n
	header				 = "t(Gyr): {}".format(roundn(time))
	data				 = np.column_stack((cartesian_coordinates, cartesian_velocities, masses))
	dir_file			 = out_dir+out_file+components[comp_id]
	np.savetxt(fname=dir_file, X=data, header=header)



print("\033[1;33m**** JOB HAS BEEN FINISHED!\033[0m")
print("Input hdf5 file:{}.hdf5 in directory:{}".format(in_file, in_dir))
print("Output text file:{}* in directory:{}".format(out_file, out_dir))

