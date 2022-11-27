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
	print("\033[1;33m**** Reading in the {} component...\033[0m".format(i))
	if i!='disk':
		pots += [agama.Potential(type="Multipole", file=indir+filename, lmax=10)]
	else:
		pots += [agama.Potential(type="CylSpline", file=indir+filename, lmax=10)]



# Calculate the basic quantities.
p_snapshot   = agama.Potential(pots[0], pots[1], pots[2])
_rmin, _rmax = 0.5, 10 #!!!
rbins        = 100
r            = np.linspace(_rmin, _rmax, rbins)
points       = np.column_stack((r, r*0, r*0))
force, deriv = p_snapshot.forceDeriv(points)
kappa        = np.sqrt(-deriv[:,0] - 3*force[:,0]/r)
nu           = np.sqrt(-deriv[:,2])
Sigma        = p_snapshot.projectedDensity(points[:, :2])


# Parameters of the subfigures.
edge   = 1.0   # the edge to avoid axis, axis' labels being truncation
width  = 10.0   # the width of the figs
height = 6   # the height of the figs
W      = width + 2*edge
H      = height + 2*edge
# final figsize


print("\033[1;33m**** Plotting RC...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
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
if fmt=='pdf':
	plt.savefig(outdir+'RCs'+".pdf") 
elif fmt=='png':
	plt.savefig(outdir+'RCs'+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** Plotting bulk density profile...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
# the total density profile
ax.plot(r, p_snapshot.density(points), color='k', label='Total')
# the density profile for different components
colors = ('r', 'y', 'grey')
i = 0
labels = ('Disk', 'Bulge', 'DM Halo')
for pot in p_snapshot:  ax.plot(r, pot.density(points), color=colors[i], label=labels[i]); i+=1
ax.set_xlabel('$R$ [kpc]')
plt.ylabel(r'$\rho$ [???]')
plt.legend(loc='upper right', fontsize=12)
plt.yscale('log')
if fmt=='pdf':
	plt.savefig(outdir+'Bulk_densitys'+".pdf") 
elif fmt=='png':
	plt.savefig(outdir+'Bulk_densitys'+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** Plotting surface density profile...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
# the total surface density profile
ax.plot(r, p_snapshot.projectedDensity(points[:, :2]), color='k', label='Total')
# the surface density profile for different components
colors = ('r', 'y', 'grey')
i = 0
labels = ('Disk', 'Bulge', 'DM Halo')
for pot in p_snapshot:  ax.plot(r, pot.projectedDensity(points[:, :2]), color=colors[i], label=labels[i]); i+=1
ax.set_xlabel('$R$ [kpc]')
plt.ylabel(r'$\Sigma$ [???]')
plt.legend(loc='upper right', fontsize=12)
plt.yscale('log')
if fmt=='pdf':
	plt.savefig(outdir+'Surface_densitys'+".pdf") 
elif fmt=='png':
	plt.savefig(outdir+'Surface_densitys'+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")


print("\033[1;33m**** Plotting Toomre parameter profile...\033[0m")
fig, ax = plt.subplots(figsize=(W, H))
# calculate the radial velocity dispersion
binsize = (_rmax - _rmin) / (rbins - 1)
rbin_edges = np.linspace(_rmin-binsize/2., _rmax+binsize/2, rbins+1)
# the rbins corresponding to the previous r values
data = np.loadtxt(indir+'model_disk')
coord = data[:, :3]
vels = data[:, 3:6]
rs = np.linalg.norm(coord, ord=2, axis=1)
vr = np.sum(vels*coord/np.column_stack((rs, rs, rs)), axis=1)
sigma_r, _,_ = bin1d(x = rs, values=vr, statistic=np.std, bins=rbin_edges)
# the radial velocity dispersion
ToomreQ = sigma_r * kappa / 3.36 / pots[0].projectedDensity(points[:, :2]) / agama.G
# the Toomre parameter
ax.plot(r, ToomreQ, color='k')
#ax.plot(r, np.ones(len(r)), color='grey')
# the surface density profile for different components
ax.set_xlabel('$R$ [kpc]')
plt.ylabel('Q')
if fmt=='pdf':
	plt.savefig(outdir+'ToomreQ'+".pdf") 
elif fmt=='png':
	plt.savefig(outdir+'ToomreQ'+".png") 
else:
	raise ValueError("The allowed file formats are only pdf and png")
