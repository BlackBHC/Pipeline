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
    "font.family":'DejaVu Serif',
	"font.size": 22,
	"mathtext.fontset": 'stix'
}
rcParams.update(config)


# Set the unit system
agama.setUnits(mass=1e10, length=1, velocity=1)


# Parse the parameters.
def parser_args():
	usage       = "Specify the directory and the filename of the output hdf5 file."
	description = "Parameters transfer from the command line into the program."
	parser      = argparse.ArgumentParser(description = description, conflict_handler = 'resolve')
	parser.add_argument("-indir", "-in", type=str, default="/home/bhchen/Temp/", help='Input directory contains the txt snapshots.')
	parser.add_argument("-outdir", "-out", type=str, default="/home/bhchen/Fig/", help='Output directory.')
	parser.add_argument("-fmt", type=str, default="pdf", help='Format of the output figs')
	return parser.parse_args()
args      = parser_args()
indir     = args.indir
outdir    = args.outdir
fmt       = args.fmt


components = ('disk', 'bulge', 'halo')

pots = []
for i in components:
	filename = 'model_'+i
	print("Read in the {} component".format(i))
	if i!='disk':
		pots += [agama.Potential(type="Multipole", file=indir+filename)]
	else:
		pots += [agama.Potential(type="CylSpline", file=indir+filename)]



# Calculate the basic quantities.
p_snapshot = agama.Potential(pots[0], pots[1], pots[2])
r            = np.linspace(0.01, 20, 100)
points       = np.column_stack((r, r*0, r*0))
force, deriv = p_snapshot.forceDeriv(points)
kappa        = np.sqrt(-deriv[:,0] - 3*force[:,0]/r)
nu           = np.sqrt(-deriv[:,2])


# Parameters of the subfigures.
edge   = 1.0   # the edge to avoid axis, axis' labels being truncation
width  = 10.0   # the width of the figs
height = 6   # the height of the figs
W      = width + 2*edge
H      = height + 2*edge
# final figsize


print("Plotting RC...")
fig, ax = plt.subplots(figsize=(14, 10))
# the total RC
ax.plot(r, np.sqrt(-r*p_snapshot.force(points)[:,0]), color='k', label='Total')
# the RC for different components
colors = ('r', 'y', 'grey')
i = 0
labels = ('Disk', 'Bulge', 'DM Halo')
for pot in p_snapshot:  ax.plot(r, np.sqrt(-r*pot.force(points)[:,0]), color=colors[i], label=labels[i]); i+=1
ax.set_xlabel('$R$ [kpc]')
plt.ylabel('$V$ [km/s]')
plt.legend(loc='upper right', fontsize=12)
plt.savefig(outdir+'RC.pdf')
if fmt=='pdf':
	plt.savefig(out_directroy+outfile+".pdf") 
elif fmt=='png':
	plt.savefig(out_directroy+outfile+".png") 
else:
	raise ValueError("The allowed file are only pdf and png")

'''print("Toomre paramter")
fig, ax = plt.subplots()
# calculate the Toomre parameter
Sigma = p_snapshot.projectedDensity(points)
ToomreQ = sigma[:,0]**0.5 * kappa / 3.36 / Sigma / agama.G
ax.plot(r, np.sqrt(-r*p_snapshot.force(points)[:,0]), color='k', label='Total')
# the RC for different components
colors = ('r', 'y', 'grey')
i = 0
labels = ('Disk', 'Bulge', 'DM Halo')
for pot in p_snapshot:  ax.plot(r, np.sqrt(-r*pot.force(points)[:,0]), color=colors[i], label=labels[i]); i+=1
plt.legend(loc='upper right', fontsize=12)
plt.savefig(outdir+'ToomreQ.pdf')

'''
