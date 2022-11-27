import agama
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
# Set the unit system
agama.setUnits(mass=1e10, length=1, velocity=1)


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
in_dir          = args.indir
out_dir         = args.outdir
in_file         = args.infile
out_file        = args.outfile
ini_root        = args.iniroot
type_simulation = args.simulationtype
type_model      = args.modeltype



# Read in configuration files
config             = cp.ConfigParser(allow_no_value=True)
config.read(ini_root+"Model.ini")
config.read(ini_root+"Figure.ini")
config.read(ini_root+"Statistic.ini")
model_type         = config[type_simulation][type_model]
components         = model_type.split()
# Whether has gas in this model?
gas_in_model       = False
for comp in components:
	if comp.lower() == "gas":
		gas_in_model = True
		break
# Whether needs to recenter the system?
recenter = config["Recenter"].getboolean("recenter")
if recenter:
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
		index        = np.where(abs(coords[:,0])<1000)[0]
		index        = index[np.where(abs(coords[index,1])<1000)[0]]
		index        = index[np.where(abs(coords[index,2])<1000)[0]]
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
		index       = np.where(abs(coords[:,0])<1000)[0]
		index       = index[np.where(abs(coords[index,1])<1000)[0]]
		index       = index[np.where(abs(coords[index,2])<1000)[0]]
		masses      = masses[index]
		coords      = coords[index]
		denominator = np.sum(masses)
		numerator   = np.sum(np.column_stack((masses, masses, masses))*coords, axis=0)
		com         = numerator/denominator
		CoMs       += [com]
	print("CoMs for each component in the model:\n", CoMs)
	# The center of mass for each component



# Parameters of the figures.
scale       = config["Figure Size"].getfloat("basic length") # the witdh and height of the x-y fig
margin      = scale * config["Figure Size"].getfloat("margin")       # the empty margin to avoid axis and axis' labels being truncated
height      = scale * config["Imshow"].getfloat("plot height")	     # the vertical height of fig, in range (0, 1)
W           = scale + 2*margin
H           = height + 2*margin
# final figsize of the canvus
fig_format  = config["Fmt"]["fig format"]
# fig format
curve_types = config["Curve"]["curve type sequence"]



# Calculate the basic quantities.
p_snapshot   = agama.Potential(pots[0], pots[1])
r            = np.linspace(rmin, rmax, rbins)
points       = np.column_stack((r, r*0, r*0))
force, deriv = p_snapshot.forceDeriv(points)
kappa        = np.sqrt(-deriv[:,0] - 3*force[:,0]/r)
nu           = np.sqrt(-deriv[:,2])
Sigma        = p_snapshot.projectedDensity(points[:, :2])




print("\033[1;33m**** Plotting RC...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
# the total RC
ax.plot(r, np.sqrt(-r*p_snapshot.force(points)[:,0]), color="k", label="Total")
# the RC for different components
colors = ("r", "y", "grey")
i = 0
labels = ("Disk", "DM Halo")
for pot in p_snapshot:  ax.plot(r, np.sqrt(-r*pot.force(points)[:,0]), color=colors[i], label=labels[i]); i+=1
ax.set_xlabel("$R$ [kpc]")
plt.ylabel("$V$ [km/s]")
plt.legend(loc="upper right", fontsize=12)
if fmt=="pdf":
	plt.savefig(outdir+"RCs"+".pdf") 
elif fmt=="png":
	plt.savefig(outdir+"RCs"+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** Plotting bulk density profile...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
# the total density profile
ax.plot(r, p_snapshot.density(points), color="k", label="Total")
# the density profile for different components
colors = ("r", "grey")
i = 0
labels = ("Disk", "DM Halo")
for pot in p_snapshot:  ax.plot(r, pot.density(points), color=colors[i], label=labels[i]); i+=1
ax.set_xlabel("$R$ [kpc]")
plt.ylabel(r"$\rho$ [???]")
plt.legend(loc="upper right", fontsize=12)
plt.yscale("log")
if fmt=="pdf":
	plt.savefig(outdir+"Bulk_densitys"+".pdf") 
elif fmt=="png":
	plt.savefig(outdir+"Bulk_densitys"+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** Plotting surface density profile...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
# the total surface density profile
ax.plot(r, p_snapshot.projectedDensity(points[:, :2]), color="k", label="Total")
# the surface density profile for different components
colors = ("r", "grey")
i = 0
labels = ("Disk", "DM Halo")
for pot in p_snapshot:  ax.plot(r, pot.projectedDensity(points[:, :2]), color=colors[i], label=labels[i]); i+=1
ax.set_xlabel("$R$ [kpc]")
plt.ylabel(r"$\Sigma$ [???]")
plt.legend(loc="upper right", fontsize=12)
plt.yscale("log")
if fmt=="pdf":
	plt.savefig(outdir+"Surface_densitys"+".pdf") 
elif fmt=="png":
	plt.savefig(outdir+"Surface_densitys"+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** Plotting Toomre parameter profile...\033[0m")
fig, ax    = plt.subplots(figsize=(W, H))
# calculate the radial velocity dispersion
binsize    = (rmax - rmin) / (rbins - 1)
rbin_edges = np.linspace(rmin-binsize/2., rmax+binsize/2, rbins+1)
# the rbins corresponding to the previous r values
data       = np.loadtxt(indir+"model_disk.txt")
coord      = data[:, :3]
vels       = data[:, 3:6]
rs         = np.linalg.norm(coord, ord=2, axis=1)
vr         = np.sum(vels*coord/np.column_stack((rs, rs, rs)), axis=1)
sigma_r, _,_ = bin1d(x = rs, values=vr, statistic=np.std, bins=rbin_edges)
# the radial velocity dispersion
ToomreQ    = sigma_r * kappa / 3.36 / pots[0].projectedDensity(points[:, :2]) / agama.G
# the Toomre parameter
ax.plot(r, ToomreQ, color="k")
ax.plot(r, np.ones(len(r)), color="r")
# the surface density profile for different components
ax.set_xlabel("$R$ [kpc]")
plt.ylabel("Q")
if fmt=="pdf":
	plt.savefig(outdir+"ToomreQ"+".pdf") 
elif fmt=="png":
	plt.savefig(outdir+"ToomreQ"+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")
