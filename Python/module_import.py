import numpy as np
import time
from matplotlib import pyplot as plt
from scipy.stats import binned_statistic_2d as bin2d
from scipy.stats import binned_statistic as bin1d
import h5py
import scipy as sp
from matplotlib import font_manager, rcParams
import random
import sys
fig_config = {
	    "font.family":'DejaVu Serif',
		"font.size": 22,
		"mathtext.fontset": 'stix'
}
rcParams.update(fig_config)
from astropy import units as u
from astropy import coordinates 
from astropy import constants as const
import astropy as ap
import configparser as cp


