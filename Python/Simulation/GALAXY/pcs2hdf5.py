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

file = "run200.pcs0"



def read_block(residual_data):
    residual_data.read(4)
    curr = np.fromfile(residual_data, dtype=np.float32, count=num[3]*num[4])
    residual_data.read(4)
    return curr
    

with open(file, 'rb') as data:
    data.read(4)
    num = np.fromfile(data, dtype=np.int32, count=5)
    time = np.fromfile(data, dtype=np.float32, count=1)[0]
    particle_mass = np.fromfile(data, dtype=np.float32, count=1)[0]
    whether_perturber = np.fromfile(data, dtype=np.bool_, count=1)[0]
    print("Particle numbers of components: ", num[:3])
    print("Info dimension: ", num[3])
    print("Lenghth of block: ", num[4])
    print("Nominal particle mass: ", particle_mass)
    print("Whether perturber: ", str(whether_perturber))
    data.seek(21, 0)

    block_test = read_block(data)

block_test = np.array(block_test).reshape(num[4], num[3])
np.savetxt("First_block.tmp", block_test)
