import agama
import numpy as np
import sys
from configparser import RawConfigParser
import argparse
import configparser as cp
agama.setUnits(mass=1e10, velocity=1, length=1)


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
in_dir          = args.indir + '/'
out_dir         = args.outdir + '/'
in_file         = args.infile
out_file        = args.outfile
ini_root        = args.iniroot + '/'
type_simulation = args.simulationtype
type_model      = args.modeltype


# Read in Model.ini configuration file
config             = cp.ConfigParser(allow_no_value=True)
config.read(ini_root+'Model.ini')
model_type         = config[type_simulation][type_model]
components         = model_type.split()
num_iteration      = config['Iteration'].getint('iterations')
num_particle_comps = []
for key in config[model_type]: 
	if key[:3]=="num":
		num_particle_comps += [config[model_type].getint(key)]
# Check whether different component are disk, in later code disk components will be dealt with axis symmetry related method
# and external DF, otherwise spherical symmetry related method with internal DF
whether_disk       = []
for key in components:
	prefix = (key[:4]).lower()
	if prefix == 'disk':
		whether_disk += [True]
	else:
		whether_disk += [False]


# Calculation the sofening length for stars and dm.
def printSoftengingLength(components):
	for i in range(len(components)):
		densComp  = model.components[i].getDensity()
		mass       = densComp.totalMass()
		num        = num_particle_comps[i]
		if components[i].lower()=='halo':
			ratio          = mass/num / (50/1e6)
			softing_length = ratio**(1/3) * 0.05
		else:
			ratio          = mass/num / (1.5/1.5e5)
			softing_length = ratio**(1/3) * 0.02

		print("Recommended value for the softening length of {} particles: {} kpc.".format(components[i],\
				softing_length))


# Print out the mass and density at some points.
def printoutInfo(model, iteration):
	print("Information of iteration {}:".format(iteration))
	for i in range(len(components)):
		densComp = model.components[i].getDensity()
		pt0 = (2.0, 0, 0)
		pt1 = (2.0, 0, 0.25)
		print("%s  total mass=%g, rho(R=2,z=0)=%g, rho(R=2,z=0.25)=%g" % \
				(components[i], densComp.totalMass(), densComp.density(pt0), densComp.density(pt1)))


# Read in the Configuration file of the galaxy model.
iniFileName = in_dir+in_file+'.ini'
ini = RawConfigParser()
ini.read(iniFileName)

iniPotenComps = []
# Potential parameters of different components 
iniSCMComps   = []
# Self consistent model parameters of different components
for comp in components:
	iniPotenComps += [dict(ini.items("Potential "+comp.lower()))]
	iniSCMComps   += [dict(ini.items("SelfConsistentModel "+comp.lower()))]
iniDFDisk          = dict(ini.items("DF disk"))

iniSCM             = dict(ini.items("SelfConsistentModel"))


# Creat the model instance.
model = agama.SelfConsistentModel(**iniSCM)

# create initial density profiles of all components
densityComps = []
for i in range(len(components)):
	densityComps += [agama.Density(**iniPotenComps[i])]


# add components to SCM - at first, all of them are static density profiles
for i in range(len(components)):
	if whether_disk[i]:
		model.components.append(agama.Component(density=densityComps[i],  disklike=True))
	else:
		model.components.append(agama.Component(density=densityComps[i],  disklike=False))


# Initialize the model, namely calculate the potential of the system.
model.iterate()
printoutInfo(model,'init')


dfComps = []
for i in range(len(components)):
	if whether_disk[i]:
		dfComps += [agama.DistributionFunction(potential=model.potential, **iniDFDisk)]
		# construct the DF of the disk component, using the initial (non-spherical) potential
	else:
		dfComps += [agama.DistributionFunction(type='QuasiSpherical', potential=model.potential, density=densityComps[i])]
		# initialize the DFs of spheroidal components using the Eddington inversion formula
		# for their respective density profiles in the initial potential


# Iteration.
massinfo =  ""
for i in range(len(components)):
	massinfo += f"\nM{components[i]}={dfComps[i].totalMass()}"
print("\033[1;33m**** STARTING ITERATIVE MODELLING ****\033[0m\nMasses (computed from DF): "+massinfo)
# replace the initially static SCM components with the DF-based ones
for i in range(len(components)):
	if whether_disk[i]:
		model.components[i] = agama.Component(df=dfComps[i],  disklike=True,  **iniSCMComps[i])
	else:
		model.components[i] = agama.Component(df=dfComps[i],  disklike=False, **iniSCMComps[i])

# do a few more iterations to obtain the self-consistent density profile for both disks
for iteration in range(1,num_iteration+1):
    print("\033[1;37mStarting iteration #%d\033[0m" % iteration)
    model.iterate()
    printoutInfo(model, 'iter%d'%iteration)
    

# export model to an N-body snapshot
print("\033[1;33mCreating an N-body representation of the model\033[0m")
format = 'text'

# now create genuinely self-consistent models of both components,
# by drawing positions and velocities from the DF in the given (self-consistent) potential
for i in range(len(components)):
	comp = components[i].lower()
	print("Sampling {} DF".format(components[i]))
	agama.writeSnapshot(out_dir+"model_{}".format(comp), \
	    agama.GalaxyModel(potential=model.potential, df=dfComps[i],  af=model.af).sample(num_particle_comps[i]), format)
# note: use a 10x larger particle mass for halo than for bulge/disk


printSoftengingLength(components)

print("Input INI file:{}.ini in directory:{}".format(in_file, in_dir))
print("Output txt file:{}* in directory:{}".format(out_file, out_dir))

