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
	parser.add_argument("-infile", "-if", type=str, default="cer", help="File name of without sufix.")
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



# Parameters from INI parameters' file
config       = cp.ConfigParser(allow_no_value=True)
config.read(ini_root+"Model.ini")
config.read(ini_root+"Figure.ini")
config.read(ini_root+"Statistic.ini")
start  = config["Gif"].getint("start subscript")
end    = config["Gif"].getint("end subscript")
fps    = config["Gif"].getint("fps")
fmt    = config["Fmt"]["fig format"]



n = start
imgs = []
while(n <= end):
	imgs += [in_dir + in_file + f"{n}."+fmt]
	n += 1



gif = []
for i in imgs:
	gif.append(imageio.imread(i, format='jpg'))
	pass

imageio.mimsave(out_dir+out_file+".gif", gif, fps=fps)

