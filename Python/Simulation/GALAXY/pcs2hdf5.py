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

file = "run999.pcs0"


raw = np.fromfile(file, dtype=np.float32, count=-1)
raw = raw[21:]

i = 0
n = (1e5+1e6) / 5e3
index = 0
while(i<n):
    if i==0:
        data = np.array(raw[index+1:index+1+5000*6]).reshape(5000,6)
    else:
        data = np.row_stack((data, np.array(raw[index+1:index+1+5000*6]).reshape(5000,6)))
    i += 1

print(len(data))
np.savetxt("initial.tmp", data)
