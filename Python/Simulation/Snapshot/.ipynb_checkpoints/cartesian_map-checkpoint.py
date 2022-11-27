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
in_dir          = args.indir
out_dir         = args.outdir
in_file         = args.infile
out_file        = args.outfile
ini_root        = args.iniroot
type_simulation = args.simulationtype
type_model      = args.modeltype



# Read in configuration files
config       = cp.ConfigParser(allow_no_value=True)
config.read(ini_root+"Model.ini")
config.read(ini_root+"Figure.ini")
config.read(ini_root+"Statistic.ini")
model_type   = config[type_simulation][type_model]
components   = model_type.split()
# Whether has gas in this model?
gas_in_model = False
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
# Interpolation
interpolation       = config["Imshow"]["interpolation type"]
interpolation_stage = config["Imshow"]["interpolation stage"]



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



# Parameters of the subfigures.
scale      = config["Figure Size"].getfloat("basic length") # the witdh and height of the x-y fig
margin     = scale * config["Figure Size"].getfloat("margin")       # the empty margin to avoid axis and axis' labels being truncated
ratio      = config["Imshow"].getfloat("z ratio")	     # the axis ratio of the x-z fig, in range (0, 1)
x_spacing  = scale * config["Figure Size"].getfloat("x spacing")		 # the horizontal gap between images of different components
y_spacing  = scale * config["Figure Size"].getfloat("y spacing")		 # the vertical gap between the x-y and x-z figs
w_cbar     = scale * config["Figure Size"].getfloat("color bar")    # the width of the color bar (in this script, all fig share the identical color bar)
cbar_shift = scale * config["Color Bar"].getfloat("shift") # the shift horizontal spacing between the color bar and imshow figs
nrows      = 2												 # the number of rows for subfigs
ncols      = len(components)								 # the number of cols for subfigs
W          = ncols*(scale+x_spacing) - x_spacing + 2*margin
H          = scale*(1+ratio) + y_spacing + 2*margin
# final figsize of the canvus
fig_format = config["Fmt"]["fig format"]



# read in parameters of the imshow images
bins         = config["Imshow"].getint("bins") # basic bin number used for imshow, bin number at z direction=ratio*bins
size	     = config["Imshow"].getfloat("size") # the box size of data to show, (-size, size) for x/y, (-size*ratio, size*ratio) for z
bin_factors  = config[model_type]["bins weight factors"].split() # bin weighting factors for different components, weight*bins(*ratio) is the final bin number used in imshow
size_factors = config[model_type]["sizes weight factors"].split() # similar as bin_factors but for size weighting
binnums      = [] # the binnum for bin2d of different comps
sizes        = [] # the box to show the different comps
for i in range(len(components)):
	binnums += [int(bins * eval(bin_factors[i]))]
	sizes   += [size * eval(size_factors[i])]
minor_ticks         = config["Axis"].getboolean("minor ticks") # whether show minor ticks
minor_ticks_density = config["Axis"]["density of minor ticks"] # density of minor ticks
basic_ticknums      = config["Axis"].getint("basic ticknums")
z_ticknums          = config["Axis"].getint("x-z vertical ticknums")


# title, text and legend parameter
fontsize      = config["Text Legend"]["fontsize"]
fontstyle     = config["Text Legend"]["fontstyle"]
fontcolor     = config["Text Legend"]["fontcolor"]
title_switch  = config["Text Legend"]["title switch"]
title_spacing = scale*config["Text Legend"].getfloat("title spacing")
	


datasets = ("column density", "cartesian velocity", "cylindrical velocity", "spherical velocity", "cartesian dispersion",\
		"cylindrical dispersion", "spherical dispersion")
# allowed color-code data types, it"s also corresponding to different cmap 
# cmap for density: afmhot for gas, jet for disk, gray_r for spherical components
# cmap for velocity: jet_r, positive are blue, negative are red
# cmap for dispersion/temperature: jet/plasma
# cmap for SFR: Reds, has not been implemented at present
spherical_velocities_ord   = ("$v_r$", "$v_\phi$", "$v_\Theta$")
cylindrical_velocities_ord = ("$v_R$", "$v_\phi$", "$v_z$")
cartesian_velocities_ord   = ("$v_x$", "$v_y$", "$v_z$")
cmaps                      = config["Imshow"]["cmap type sequence"].split()



def statistic_component(comp_id, dataset, velocity_id):
	# the function to calculate the image matrices of different component with different target dataset
	# only for internal use for show_component() function
	# comp_id for subscript of comps
	# dataset in datasets
	
	# diagnostic value check
	if (dataset in datasets)==False:
		raise ValueError("dataset should be one of,\n{}\nBut give dataset={}.".format(datasets,dataset))
	if (velocity_id in (0, 1, 2)) == False:
		raise ValueError("Velocity id should be 0, 1 or 2, but give value={}.".format(velocity_id))
	
	if recenter and not(recenter_diff_comps):
		n = len(comps[comp_id]["Coordinates"][...])
		com = np.column_stack((np.ones(n)*CoM[0], np.ones(n)*CoM[1], np.ones(n)*CoM[2]))
		cartesian_coordinates = comps[comp_id]["Coordinates"][...] - com
	elif recenter:
		n = len(comps[comp_id]["Coordinates"][...])
		com = np.column_stack((np.ones(n)*CoMs[comp_id][0], np.ones(n)*CoMs[comp_id][1], np.ones(n)*CoMs[comp_id][2]))
		cartesian_coordinates = comps[comp_id]["Coordinates"][...] - com
	else:
		cartesian_coordinates = comps[comp_id]["Coordinates"][...]	
	cartesian_velocities  = comps[comp_id]["Velocities"][...] 
	if dataset=="column density":
		image_xy, _,_,_ = bin2d(x=cartesian_coordinates[:,1], y=cartesian_coordinates[:,0],\
				values=cartesian_coordinates[:,0], bins=binnums[comp_id], range=[[-sizes[comp_id], \
				sizes[comp_id]],[-sizes[comp_id], sizes[comp_id]]], statistic="count")
		image_xz, _,_,_ = bin2d(x=cartesian_coordinates[:,2], y=cartesian_coordinates[:,0],\
				values=cartesian_coordinates[:,0], bins=[ratio*binnums[comp_id], binnums[comp_id]],\
				range=[[-sizes[comp_id]*ratio, sizes[comp_id]*ratio], [-sizes[comp_id], sizes[comp_id]]], statistic="count")
		index_xy  = np.where(image_xy<1)
		image_xy[index_xy] = 1
		image_xy  = np.log10(image_xy)
		image_xy /= image_xy.max()
		index_xz  = np.where(image_xz<1)
		image_xz[index_xz] = 1
		image_xz  = np.log10(image_xz)
		image_xz /= image_xz.max()

		if config["Imshow"].getboolean("show empty"):
			# show the pixel without data as write
			image_xy[index_xy] = None
			image_xz[index_xz] = None
		return image_xy, image_xz

	# statistic method: dispersion or mean for velocity
	if dataset[-4:]=="sion":
		method = "std"
	else:
		method = "mean"

	# in which frame to be shown
	if dataset[:3]=="car":
		# cartesian velocity
		value = cartesian_velocities[:, velocity_id]
	elif dataset[:3]=="sph":
		# spherical velocity
		radius   = np.linalg.norm(cartesian_coordinates, axis=1, ord=2)
		theta    = np.arcsin(cartesian_coordinates[:,2]/radius)
		phi      = np.arccos(cartesian_coordinates[:,0]/
				np.linalg.norm(cartesian_coordinates[:,:2], ord=2, axis=1))
		unit_r   = cartesian_coordinates / np.column_stack((radius, radius, radius))
		rotation = np.array([[0, -1],\
				[1, 0]])
		unit_phi = np.matmul(rotation, unit_r[:,:2].T).T
		unit_phi = np.column_stack((unit_phi[:,0], unit_phi[:,1], np.zeros(len(unit_phi))))
		# the unit azimuthal vector
		if velocity_id==0:
			v_r     = np.sum(unit_r*cartesian_velocities, axis=1)
			value   = v_r
		elif velocity_id==1:
			v_phi   = np.sum(unit_phi*cartesian_velocities, axis=1)
			value   = v_phi
		else:
			v_r     = np.sum(unit_r*cartesian_velocities, axis=1)
			v_phi   = np.sum(unit_phi*cartesian_velocities, axis=1)
			v_theta = np.sum(cartesian_velocities - unit_r*np.column_stack((v_r, v_r, v_r))-\
					unit_phi*np.column_stack((v_phi, v_phi, v_phi)), axis=1)
			value   = v_theta
	else:
		# cylindrical velocity
		radius   = np.linalg.norm(cartesian_coordinates[:,:2], axis=1, ord=2)
		phi      = np.arccos(cartesian_coordinates[:,0]/radius)
		unit_r   = cartesian_coordinates[:,:2] / np.column_stack((radius, radius))
		rotation = np.array([[0, -1], [1, 0]])
		unit_phi = np.matmul(rotation, unit_r[:,:2].T).T
        # the unit azimuthal vector
		if velocity_id==0:
			v_r   = np.sum(unit_r*cartesian_velocities[:,:2], axis=1)
			value = v_r
		elif velocity_id==1:
			v_phi = np.sum(unit_phi*cartesian_velocities[:, :2], axis=1)
			value = v_phi
		else:
			value = cartesian_velocities[:,2]

	image_xy, _,_,_ = bin2d(x=cartesian_coordinates[:,1], y=cartesian_coordinates[:,0], values=value,\
			bins=binnums[comp_id], range=[[-sizes[comp_id], sizes[comp_id]],[-sizes[comp_id], sizes[comp_id]]],\
			statistic=method)
	image_xz, _,_,_ = bin2d(x=cartesian_coordinates[:,2], y=cartesian_coordinates[:,0], values=value,\
			bins=[ratio*binnums[comp_id], binnums[comp_id]], range=[[-sizes[comp_id]*ratio, sizes[comp_id]*ratio],\
			[-sizes[comp_id], sizes[comp_id]]], statistic=method)

	return image_xy, image_xz
	


def show_component(comp_id, dataset, axes, axes_cbars, velocity_id=0, interpolation=interpolation, \
		interpolation_stage=interpolation_stage):
	# show the image matrices created by statistic_component() function

	# the cmap usde for imshow
	if (dataset in datasets)==False:
		raise ValueError("dataset should be one of:\n{}\n but given dataset={}".format(datasets, dataset))
	elif dataset=="column density":
		if gas_in_model and comp_id==0:
			cmap=cmaps[4]
		elif components[comp_id][:4].lower()=="disk":
			cmap=cmaps[0]
		else:
			cmap=cmaps[3]
	elif dataset[-4:]=="sion":
		cmap = cmaps[6]
	else:
		cmap = cmaps[0]
	
	
	# calculate the statistical images
	image_xy, image_xz = statistic_component(comp_id=comp_id, dataset=dataset, velocity_id=velocity_id)
	
	# inshow and colorbar
	# the color bar range
	min_xy = np.min(image_xy[np.where(image_xy>-1e20)])
	max_xy = np.max(image_xy[np.where(image_xy>-1e20)])	
	min_xz = np.min(image_xz[np.where(image_xz>-1e20)])
	max_xz = np.max(image_xz[np.where(image_xz>-1e20)])	
	vmin   = min(min_xy, min_xz)
	vmax   = max(max_xy, max_xz)
	axes[0, comp_id].imshow(image_xy, origin="lower", cmap=cmap, interpolation=interpolation,\
			interpolation_stage=interpolation_stage, vmin=vmin, vmax=vmax)
	im = axes[1, comp_id].imshow(image_xz, origin="lower", cmap=cmap, interpolation=interpolation,\
			interpolation_stage=interpolation_stage, vmin=vmin, vmax=vmax)
	plt.colorbar(mappable=im, cax=axes_cbars[comp_id])
	# labels
	axes[0, comp_id].set_xlabel("$X$ [kpc]")
	axes[0, comp_id].set_ylabel("$Y$ [kpc]")
	axes[1, comp_id].set_xlabel("$X$ [kpc]")
	axes[1, comp_id].set_ylabel("$Z$ [kpc]")
	# ticks 
	axes[0, comp_id].set_xticks(np.linspace(0, binnums[comp_id]-1, basic_ticknums))
	axes[0, comp_id].set_xticklabels(np.linspace(-sizes[comp_id], sizes[comp_id],basic_ticknums))
	axes[0, comp_id].set_yticks(np.linspace(0, binnums[comp_id]-1, basic_ticknums))
	axes[0, comp_id].set_yticklabels(np.linspace(-sizes[comp_id], sizes[comp_id],basic_ticknums))
	axes[1, comp_id].set_xticks(np.linspace(0, binnums[comp_id]-1, basic_ticknums))
	axes[1, comp_id].set_xticklabels(np.linspace(-sizes[comp_id], sizes[comp_id],basic_ticknums))
	axes[1, comp_id].set_yticks(np.linspace(0, ratio*binnums[comp_id]-1, z_ticknums))
	axes[1, comp_id].set_yticklabels(np.linspace(-sizes[comp_id]*ratio, sizes[comp_id]*ratio, z_ticknums))
	# minor ticks
	if minor_ticks:
		axes[0, comp_id].minorticks_on()
		axes[0, comp_id].tick_params(axis="both", which="both", direction="in", labelsize=minor_ticks_density)
		axes[1, comp_id].minorticks_on()
		axes[1, comp_id].tick_params(axis="both", which="both", direction="in", labelsize=minor_ticks_density)
	# information text
	mass    = comps[comp_id]["Masses"][...]
	par_num = len(mass)
	mass    = mass.sum()
	decimal = lambda x, n=3: int(x*10**n)/10**n
	axes[0, comp_id].text(binnums[comp_id]*0.05,binnums[comp_id]*0.85,\
			"Time: {}Gyr\nMass: {}e10$\ M_\odot$\nParticle Numbders: {}".format(time, decimal(x=mass), par_num),\
			fontsize=fontsize, fontstyle=fontstyle, color=fontcolor)
	pass



# the function to create canvus and imshow the desired quantities
def fig_comps_method(dataset, velocity_id=0, interpolation=interpolation, interpolation_stage=interpolation_stage):
	# Create the canvus and set positions of different panels and colorbars.
	fig, axes  = plt.subplots(nrows=nrows, ncols=ncols, figsize=(W, H))
	axes_cbars = []
	for i in range(ncols):
		axes[0, i].set_position([(i*(scale+x_spacing)+margin)/W, (y_spacing + ratio*scale + margin)/H, scale/W, scale/H])
		axes[1, i].set_position([(i*(scale+x_spacing)+margin)/W, margin/H, scale/W, scale*ratio/H])
		axes_cbars += [fig.add_axes([(i*(scale+x_spacing)+margin+scale+cbar_shift)/W, margin/H, w_cbar/W, (H-2*margin)/H])]
	
	for i in range(len(comps)):
		show_component(comp_id=i, dataset=dataset, velocity_id=velocity_id, axes=axes, axes_cbars=axes_cbars)
	if dataset=="column density":
		print("\033[1;33m**** Plotting cartesian maps of {}.\033[0m".format(dataset))
		if title_switch: plt.suptitle("Relative Column Density", x=.5, y=(margin+(1+ratio)*scale+y_spacing+title_spacing)/H,\
				ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
		plt.savefig(out_dir+out_file+"_"+dataset.replace(" ", "_")+"."+fig_format)
	else:
		if dataset[-4:]=="city":
			if dataset[:3]=="car" and title_switch:
				print("\033[1;33m**** Plotting cartesian maps of {}.\033[0m".format(\
						cartesian_velocities_ord[velocity_id][1:-1]))
				plt.suptitle("{}".format(cartesian_velocities_ord[velocity_id]), x=.5, y=(margin+(1+ratio)*scale+y_spacing+\
						title_spacing)/H, ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
			elif dataset[:3]=="cyl" and title_switch:
				print("\033[1;33m**** Plotting cartesian maps of {}.\033[0m".format(\
						cylindrical_velocities_ord[velocity_id][1:-1]))
				plt.suptitle("{}".format(cylindrical_velocities_ord[velocity_id]), x=.5, y=(margin+(1+ratio)*scale+y_spacing+\
						title_spacing)/H, ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
			elif dataset[:3]=="sph" and title_switch:
				print("\033[1;33m**** Plotting cartesian maps of {}.\033[0m".format(\
						spherical_velocities_ord[velocity_id][1:-1]))
				plt.suptitle("{}".format(spherical_velocities_ord[velocity_id]), x=.5, y=(margin+(1+ratio)*scale+y_spacing+\
						title_spacing)/H, ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
		else:
			if dataset[:3]=="car" and title_switch:
				print("\033[1;33m**** Plotting cartesian maps of dispersion of {}.\033[0m".format(\
						cartesian_velocities_ord[velocity_id][1:-1]))
				plt.suptitle("$\sigma$({})".format(cartesian_velocities_ord[velocity_id]), x=.5, y=(margin+(1+ratio)*scale+\
						y_spacing+title_spacing)/H, ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
			elif dataset[:3]=="cyl" and title_switch:
				print("\033[1;33m**** Plotting cartesian maps of dispersion of {}.\033[0m".format(\
						cylindrical_velocities_ord[velocity_id][1:-1]))
				plt.suptitle("$\sigma$({})".format(cylindrical_velocities_ord[velocity_id]), x=.5, y=(margin+(1+ratio)*scale\
						+y_spacing+title_spacing)/H, ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
			elif dataset[:3]=="sph" and title_switch:
				print("\033[1;33m**** Plotting cartesian maps of dispersion of {}.\033[0m".format(\
						spherical_velocities_ord[velocity_id][1:-1]))
				plt.suptitle("$\sigma$({})".format(spherical_velocities_ord[velocity_id]), x=.5, y=(margin+(1+ratio)*scale\
						+y_spacing+title_spacing)/H, ha="center", va="center", fontsize=fontsize, fontstyle=fontstyle)
		plt.savefig(out_dir+out_file+"_"+dataset.replace(" ", "_")+"_{}.".format(velocity_id)+fig_format)



# imshow different components
for dataset in datasets:
	if config["Imshow"].getboolean(dataset)==False:
		continue
	if dataset=="column density":
		fig_comps_method(dataset)
	else:
		for i in range(3):fig_comps_method(dataset, velocity_id=i)


