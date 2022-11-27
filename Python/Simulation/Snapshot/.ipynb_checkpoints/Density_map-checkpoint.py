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
    "font.family":'DejaVu Serif',
	"font.size": 22,
	"mathtext.fontset": 'stix'
}
rcParams.update(config)


# Parse the parameters.
def parser_args():
	usage       = "Specify directory and filename of related parameter files."
	description = "Parameters transfer from the command line into the program."
	parser      = argparse.ArgumentParser(description = description, conflict_handler = 'resolve')
	parser.add_argument("-indir", '-id', type=str, default="/home/bhchen/BPX_Populations/IC/", help='Input directory where \
			has the Config.ini file.')
	parser.add_argument("-outdir", '-od', type=str, default="/home/bhchen/BPX_Populations/IC/test/", help='Output directory where \
			to output the IC txt file.')
	parser.add_argument("-infile", '-if', type=str, default="Config", help='File name of the \
			configuration file without sufix .ini')
	parser.add_argument("-outfile", '-of', type=str, default="model_", help='Prefix of the output txt file, for example \
			if has disk and halo components with the default value, the output files are model_disk.txt and\
			model_halo.txt.')
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
for key in config[model_type]: num_particle_comps += [config[model_type].getint(key)]


if recenter==1:
	print("\033[1;33m**** WILL RECENTER THE SYSTEM FOR EACH COMPONENT!\033[0m")
	print("\033[1;33m**** NOTE THE CENTER OF EACH COMPOENT MAY NOT COINCIDE WITH EACH OTHER!\033[0m")


# Read in the snapshot
print("\033[1;33m**** READ IN SNAPSHOT\033[0m")
f = h5py.File(directory+infile+'.hdf5', 'r')
# the time of the snapshot
time = f['Header'].attrs['Time']
# number of PartType
num_comps = 0
comp_id   = []
for i in f.keys():
	if i[:8]=='PartType':
		if gas==1:
			num_comps += 1
			comp_id   += [i[8:]]
		else:
			if i[8]!='0':
				num_comps += 1
				comp_id   += [i[8:]]
	else: pass

comps = []
for i in range(num_comps): comps += [f['PartType'+comp_id[i]]]


print("\033[1;33m**** PLOTTING SPATIAL DISTRIBUTION FIGS \033[0m")
# Parameters of the subfigures.
edge   = 2.0   # the edge to avoid axis, axis' labels being truncation
scale  = 8.0   # the witdh and height of the face-on fig
ratio  = .4    # the axis ratio of the edge-on fig, in range (0, 1)
v_gap  = .5    # the vertical gap between the face-on and edge-on figs
h_gap  = 2.5   # the horizontal gap between images of different components
w_cbar = .5    # the width of the color bar (in this script, all fig share the identical color bar)
nrows  = 2     # the number of rows for subfigs
ncols  = 2     # the number of cols for subfigs
W      = ncols*(scale+h_gap) - h_gap + 2*edge
H      = scale*(1+ratio) + v_gap + 2*edge
# final figsize


# Create the canvus and set positions.
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(W, H))
for j in range(ncols):
	axes[0, j].set_position([(j*(scale+h_gap)+edge)/W, (v_gap + ratio*scale + edge)/H, scale/W, scale/H])
	axes[1, j].set_position([(j*(scale+h_gap)+edge)/W, edge/H, scale/W, scale*ratio/H])


# plt.imshow the images
binnum = [100, 100] # the binnum for bin2d of different comps
size   = [10,  100] # the box to show the different comps
for n in range(num_comps):
	coord = comps[n]['Coordinates'][...]
	vels  = comps[n]['Velocities'][...]
	if recenter==1 and n==0:
		coord -= np.mean(coord, axis=0)
	num_p  = len(coord)
	image_xy, _,_,_ = bin2d(x=coord[:, 1], y=coord[:, 0], values=coord[:, 0], statistic='count',\
			range = [[-size[n], size[n]], [-size[n], size[n]]], bins=binnum[n])
	image_xz, _,_,_ = bin2d(x=coord[:, 2], y=coord[:, 0], values=coord[:, 0], statistic='count',\
			range = [[ratio*-size[n], ratio*size[n]], [-size[n], size[n]]], bins=[binnum[n]*ratio, binnum[n]])
	
	# calculate the relative density
	index_xy           = np.where(image_xy<1)
	index_xz           = np.where(image_xz<1)
	image_xy[index_xy] = 1
	image_xy           = np.log10(image_xy)
	image_xy          /= np.max(image_xy)
	#image_xy[index_xy] = None # save the no particle bins
	image_xz[index_xz] = 1
	image_xz           = np.log10(image_xz)
	image_xz          /= np.max(image_xz)
	#image_xz[index_xz] = None # save the no particle bins

	im = axes[0, n].imshow(image_xy, cmap='jet', origin='lower', interpolation='gaussian')
	im = axes[1, n].imshow(image_xz, cmap='jet', origin='lower', interpolation='gaussian')

	# set the axis labels
	axes[0, n].set_xticks(np.linspace(0, binnum[n]-1, 9))
	axes[0, n].set_xticklabels(np.around(np.linspace(-size[n], size[n], 9), 2))
	axes[1, n].set_xticks(np.linspace(0, binnum[n]-1, 9))
	axes[1, n].set_xticklabels(np.around(np.linspace(-size[n], size[n], 9), 2))
	axes[1, n].set_xlabel("$X$ [kpc]")
	axes[0, n].set_yticks(np.linspace(0, binnum[n]-1, 9))
	axes[0, n].set_yticklabels(np.around(np.linspace(-size[n], size[n], 9), 2))
	axes[1, n].set_yticks(np.linspace(0, ratio*binnum[n]-1, 5))
	axes[1, n].set_yticklabels(np.around(np.linspace(ratio*-size[n], ratio*size[n], 5), 2))
	axes[0, n].set_ylabel("$Y$ [kpc]")
	axes[1, n].set_ylabel("$Z$ [kpc]")
	
	# the info of the components
	axes[0, n].text(.05*binnum[n], .85*binnum[n], "Parcicle Num: {num}\nT: {t}Gyr\nMass: {m}e10$M_\odot$"\
			.format(num=num_p, t=np.around(time, 3), m=int(comps[n]['Masses'][...].sum()*1000)/1000),\
			color='r')

cax = fig.add_axes([(edge+ncols*(scale+h_gap)-h_gap*.6)/W, edge/H, w_cbar/W, (scale*(1+ratio) + v_gap)/H])
plt.colorbar(mappable = im, cax = cax)

print("\033[1;33m**** SAVE FIG\033[0m")
fig.suptitle("Density Map")
if fmt=='pdf':
	plt.savefig(out_directroy+outfile+".pdf") 
elif fmt=='png':
	plt.savefig(out_directroy+outfile+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** CLOSE SNAPSHOT FILE\033[0m")
f.close()
print("\033[1;33m**** DONE!\033[0m")
