import agama
import numpy as np
from matplotlib import pyplot as plt
import h5py
from scipy.stats import binned_statistic as bin1d
from scipy.stats import binned_statistic_2d as bin2d
import scipy as sp
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
	parser.add_argument("-indir", "-id", type=str, default="/home/bhchen/Temp/", help="Input directory where \
			has the Config.ini file.")
	parser.add_argument("-outdir", "-od", type=str, default="/home/bhchen/Fig/", help="Output directory.")
	parser.add_argument("-infile", "-if", type=str, default="model", help="File name of without sufix.")
	parser.add_argument("-outfile", "-of", type=str, default="radial_profile", help="Prefix of the output file.")
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
model_type         = config[type_simulation][type_model]
components         = model_type.split()
# Whether has gas in this model?
gas_in_model       = False
for comp in components:
	if comp.lower() == "gas":
		gas_in_model = True
		break



# Parameters of the figures.
scale       = config["Figure Size"].getfloat("basic length") # the witdh and height of the x-y fig
margin      = scale * config["Figure Size"].getfloat("margin")       # the empty margin to avoid axis and axis' labels being truncated
height      = scale * config["Figure Size"].getfloat("plot height")	     # the vertical height of fig, in range (0, 1)
W           = scale + 2*margin
H           = height + 2*margin
# final figsize of the canvus
fig_format  = config["Fmt"]["fig format"]
# fig format



# Parameters of the curves
curve_types  = config["Curve"]["curve type sequence"].split()
curve_colors = config["Curve"]["curve color sequence"].split()



# Fontsize
fontsize = config["Text Legend"]["fontsize"]



# Calculate the potential with multipole and cylspine expansion.
pots = []
for i,component in enumerate(components):
	dir_file=in_dir+in_file+component
	print("\033[1;33m**** Reading in the {} component...\033[0m".format(component))
	if component[:4].lower() == "disk":
		pots += [agama.Potential(type="CylSpline", file=dir_file, lmax=10)]
	elif gas_in_model and i==0:
		pots += [agama.Potential(type="CylSpline", file=dir_file, lmax=10)]
	else:
		pots += [agama.Potential(type="Multipole", file=dir_file, lmax=10)]



# Creat the potential instance of the snapshot
args = ''
for i in range(len(components)):
	args += 'pots[{}],'.format(i) 
exec("p_snapshot = agama.Potential({})".format(args))



# Calculate the basic quantities.
rmin		 = config["Radial Profile"].getfloat("r min")
rmax		 = config["Radial Profile"].getfloat("r max")
rbins		 = config["Radial Profile"].getint("bins") 
r            = np.linspace(rmin, rmax, rbins)
points       = np.column_stack((r, r*0, r*0))
force, deriv = p_snapshot.forceDeriv(points)
kappa        = np.sqrt(-deriv[:,0] - 3*force[:,0]/r)
nu           = np.sqrt(-deriv[:,2])
Sigma        = p_snapshot.projectedDensity(points[:, :2])



if config["Radial Profile"].getboolean("rotation curve"):
	print("\033[1;33m**** Plotting RC...\033[0m")
	fig, ax = plt.subplots(figsize=(W, H))
	# the total RC
	ax.plot(r, np.sqrt(-r*p_snapshot.force(points)[:,0]), color=curve_colors[-1], linestyle=curve_types[0], label="Total")
	# the RC for different components
	for i,pot in enumerate(p_snapshot):
		ax.plot(r, np.sqrt(-r*pot.force(points)[:,0]), color=curve_colors[i], linestyle=curve_types[0], label=components[i])
	ax.set_xlabel("$R$ [kpc]", fontsize=fontsize)
	ax.set_ylabel("$V_C$ [km/s]", fontsize=fontsize)
	plt.legend(loc="upper right", fontsize="medium")
	dir_file = out_dir+out_file+"_RC"
	plt.savefig(dir_file+".pdf")


if config["Radial Profile"].getboolean("bulk density"):
	print("\033[1;33m**** Plotting bulk density profile...\033[0m")
	fig, ax = plt.subplots(figsize=(W, H))
	# the total density profile
	ax.plot(r, p_snapshot.density(points), color=curve_colors[-1], linestyle=curve_types[0], label="Total")
	# the density profile for different components
	for i,pot in enumerate(p_snapshot):
		ax.plot(r, pot.density(points), color=curve_colors[i], linestyle=curve_types[0], label=components[i])
	ax.set_xlabel("$R$ [kpc]", fontsize=fontsize)
	ax.set_ylabel(r"$\rho\ \rm{[M_\odot/pc^3]}$", fontsize=fontsize)
	plt.legend(loc="upper right", fontsize=fontsize)
	plt.yscale("log")
	dir_file = out_dir+out_file+"_rho"
	plt.savefig(dir_file+".pdf")
	

if config["Radial Profile"].getboolean("surface density"):
	print("\033[1;33m**** Plotting surface density profile...\033[0m")
	fig, ax = plt.subplots(figsize=(W, H))
	# the total surface density profile
	ax.plot(r, p_snapshot.projectedDensity(points[:,:2]), color=curve_colors[-1], linestyle=curve_types[0], label="Total")
	# the surface density profile for different components
	labels = ("Disk", "DM Halo")
	for i,pot in enumerate(p_snapshot): 
		ax.plot(r, pot.projectedDensity(points[:,:2]), color=curve_colors[i], linestyle=curve_types[0], label=components[i])
	ax.set_xlabel("$R$ [kpc]", fontsize=fontsize)
	ax.set_ylabel(r"$\Sigma\ \rm{[k\,M_\odot/pc^2]}$", fontsize=fontsize)
	plt.legend(loc="upper right", fontsize=fontsize)
	plt.yscale("log")
	dir_file = out_dir+out_file+"_Sigma"
	plt.savefig(dir_file+".pdf")



if config["Radial Profile"].getboolean("toomre parameter"):
	print("\033[1;33m**** Plotting Toomre parameter profile...\033[0m")
	fig, ax    = plt.subplots(figsize=(W, H))
	# calculate the radial velocity dispersion
	binsize    = (rmax - rmin) / (rbins - 1)
	rbin_edges = np.linspace(rmin-binsize/2., rmax+binsize/2, rbins+1)
	# the rbins corresponding to the previous r values
	dir_file   = in_dir+in_file+"Disk"
	data       = np.loadtxt(dir_file)
	coord      = data[:, :3]
	vels       = data[:, 3:5]
	rs         = np.linalg.norm(coord, ord=2, axis=1)
	vr         = np.sum(vels*coord[:, :2]/np.column_stack((rs, rs)), axis=1)
	sigma_r, _,_ = bin1d(x = rs, values=vr, statistic="std", bins=rbin_edges)
	# the radial velocity dispersion
	ToomreQ    = sigma_r * kappa / 3.36 / pots[0].projectedDensity(points[:, :2]) / agama.G
	# the Toomre parameter
	if config["Radial Profile"].getboolean("smooth toomre"):
		ToomreQ = sp.signal.savgol_filter(x=ToomreQ, window_length=int(rbins/2), polyorder=8)
	ax.plot(r, ToomreQ, color="k")
	ax.plot(r, np.ones(len(r)), color="r")
	# the surface density profile for different components
	ax.set_xlabel("$R$ [kpc]", fontsize=fontsize)
	ax.set_ylabel("Q", fontsize=fontsize)
	dir_file = out_dir+out_file+"_Q"
	plt.savefig(dir_file+".pdf")

