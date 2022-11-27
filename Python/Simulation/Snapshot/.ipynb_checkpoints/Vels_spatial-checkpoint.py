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
	usage       = "Specify the directory and the filename of the output hdf5 file."
	description = "Parameters transfer from the command line into the program."
	parser      = argparse.ArgumentParser(description = description, conflict_handler = 'resolve')
	parser.add_argument("-indir", type=str, default="/home/bhchen/BPX_Populations/IC/", help='Output directory.')
	parser.add_argument("-outdir", type=str, default="./", help='Output directory.')
	parser.add_argument("-filename", type=str, default="output", help='Name of the snapshot file without suffix .hdf5')
	parser.add_argument("-recenter", type=np.uint32, default=1, help='Whether recenter the particles: 1=Yes, 0=Not.')
	parser.add_argument("-gas", type=np.uint32, default=0, help='Whether show gas particles: 1=Yes, 0=Not.')
	return parser.parse_args()
args          = parser_args()
directory     = args.indir
out_directroy = args.outdir
filename      = args.filename
recenter      = args.recenter
gas           = args.gas


# Read in the snapshot
print("\033[1;33m**** READ IN SNAPSHOT\033[0m")
f = h5py.File(directory+filename+'.hdf5', 'r')
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


print("\033[1;33m**** PLOTTING VELOCITY SPATIAL DISTRIBUTION FIGS \033[0m")

vtype = ('x', 'y', 'z', '')
for m in range(len(vtype)):
	# Parameters of the subfigures.
	edge   = 2.0   # the edge to avoid axis, axis' labels being truncation
	scale  = 8.0   # the witdh and height of the face-on fig
	ratio  = .4    # the axis ratio of the edge-on fig, in range (0, 1)
	v_gap  = .5    # the vertical gap between the face-on and edge-on figs
	h_gap  = 3.0   # the horizontal gap between images of different components
	w_cbar = .4    # the width of the color bar (in this script, all fig share the identical color bar)
	nrows  = 2     # the number of rows for subfigs
	ncols  = 3     # the number of cols for subfigs
	W      = ncols*(scale+h_gap) - h_gap + 2*edge
	H      = scale*(1+ratio) + v_gap + 2*edge
	# final figsize
	
	
	# Create the canvus and set positions.
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols*2, figsize=(W, H))
	for j in range(ncols):
		axes[0, 2*j].set_position([(j*(scale+h_gap)+edge)/W, (v_gap + ratio*scale + edge)/H, scale/W, scale/H])
		axes[1, 2*j].set_position([(j*(scale+h_gap)+edge)/W, edge/H, scale/W, scale*ratio/H])
		axes[0, 2*j+1].set_position([(j*(scale+h_gap)+edge+.025*h_gap+scale)/W, (v_gap + ratio*scale + edge)/H, w_cbar/W, scale/H])
		axes[1, 2*j+1].set_position([(j*(scale+h_gap)+edge+.025*h_gap+scale)/W, edge/H, w_cbar/W, scale*ratio/H])

	# plt.imshow the images
	binnum = [100, 100, 100] # the binnum for bin2d of different comps
	size   = [10,   2,  50] # the box to show the different comps
	for n in range(num_comps):
		coord  = comps[n]['Coordinates'][...]
		vels   = comps[n]['Velocities'][...]
		vs     = np.linalg.norm(vels, ord=2, axis=1)
		num_p  = len(coord)
		
		# calculate the velocity maps
		if vtype[m] != '':
			image_xy, _,_,_ = bin2d(x=coord[:, 1], y=coord[:, 0], values=vels[:, m], statistic='mean',\
					range = [[-size[n], size[n]], [-size[n], size[n]]], bins=binnum[n])
			image_xz, _,_,_ = bin2d(x=coord[:, 2], y=coord[:, 0], values=vels[:, m], statistic='mean',\
					range = [[ratio*-size[n], ratio*size[n]], [-size[n], size[n]]], bins=[binnum[n]*ratio, binnum[n]])
		else:
			image_xy, _,_,_ = bin2d(x=coord[:, 1], y=coord[:, 0], values=vs, statistic='mean',\
					range = [[-size[n], size[n]], [-size[n], size[n]]], bins=binnum[n])
			image_xz, _,_,_ = bin2d(x=coord[:, 2], y=coord[:, 0], values=vs, statistic='mean',\
					range = [[ratio*-size[n], ratio*size[n]], [-size[n], size[n]]], bins=[binnum[n]*ratio, binnum[n]])
		
		im1 = axes[0, 2*n].imshow(image_xy, cmap='jet', origin='lower')
		im2 = axes[1, 2*n].imshow(image_xz, cmap='jet', origin='lower')
		plt.colorbar(mappable=im1, cax=axes[0, 2*n+1])
		plt.colorbar(mappable=im2, cax=axes[1, 2*n+1])

		# set the axis labels
		axes[0, 2*n].set_xticks(np.linspace(0, binnum[n]-1, 9))
		axes[0, 2*n].set_xticklabels(np.around(np.linspace(-size[n], size[n], 9), 2))
		axes[1, 2*n].set_xticks(np.linspace(0, binnum[n]-1, 9))
		axes[1, 2*n].set_xticklabels(np.around(np.linspace(-size[n], size[n], 9), 2))
		axes[1, 2*n].set_xlabel("$X$ [kpc]")
		axes[0, 2*n].set_yticks(np.linspace(0, binnum[n]-1, 9))
		axes[0, 2*n].set_yticklabels(np.around(np.linspace(-size[n], size[n], 9), 2))
		axes[1, 2*n].set_yticks(np.linspace(0, ratio*binnum[n]-1, 5))
		axes[1, 2*n].set_yticklabels(np.around(np.linspace(ratio*-size[n], ratio*size[n], 5), 2))
		axes[0, 2*n].set_ylabel("$Y$ [kpc]")
		axes[1, 2*n].set_ylabel("$Z$ [kpc]")
        	
        # the info of the components
        axes[0, 2*n].text(.05*binnum[n], .85*binnum[n], "Parcicle Num: {num}\nT: {t}Gyr\nMass: {m}e10$M_\odot$"\
			.format(num=len(vels), t=np.round(time, 2), m=np.round(comps[n]['Masses'][...].sum(), 2)),\
			color='r')
	
	print("\033[1;33m**** SAVE FIG FOR V{}\033[0m".format(vtype[m]))
	plt.savefig(out_directroy+'v'+vtype[m]+'_map.pdf') 


print("\033[1;33m**** CLOSE SNAPSHOT FILE\033[0m")
f.close()
print("\033[1;33m**** DONE!\033[0m")
