import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline, interp1d
import scipy.stats 
import sys
import os
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import lognorm

# Code to extract infor from the high res simulation 4096^3 coarsened to 256^3
# IMPORTANT: The part of reading the *halo_mass256.bin file might be editted depending on
# the way it has been saved. In the below form we are assuming that the information about each
# subvolume is taking up two lines in the file.

count = 0
sz = '%.3f' %7.570
txt_file = open('/home/davida/Documents/Sussex PhD/research/subgrid_sources/testing/7.570halo_mass256.bin')
fcoll = []
dens = []

for line in txt_file:
    # we can process file line by line here, for simplicity I am taking count of lines
    elements = line.split()
    if len(elements)==5:
        dens_value = float(elements[0])
        dens.append(dens_value)
    else:
        fcoll_value = float(elements[0])
        fcoll.append(fcoll_value)
    
txt_file.close()

# The density values are stored in the 1st entry of every two lines
dens = dens[0::2]
overdense = np.array(dens/np.average(dens))
fcoll = np.array(fcoll)


# When fitting the lognormals we need to remove any fcoll=0 values
delta = overdense-1
delta_nonzero = delta[fcoll!=0]
fcoll_nonzero = fcoll[fcoll!=0]

# Define the bins for overdensities, starting from -1 to 10 with a bin width of 0.1
bin_edges = np.arange(-1, 10 + 0.1, 0.1)
# For Larger overdensities we have wider bins for sufficient statistics
bin_edges = np.append(bin_edges, 15)
bin_edges = np.append(bin_edges, 30)
bin_edges = np.append(bin_edges, 60)
bin_edges = np.append(bin_edges, 80)
bin_edges = np.append(bin_edges, 100)

bin_edges = np.append(bin_edges, 100)


# Digitize the overdensity values to determine which bin they fall into
bin_indices = np.digitize(delta_nonzero, bin_edges)

# Initialize an array to store the mean collapse fraction for each bin
delta_collapse_per_bin = np.zeros(len(bin_edges) - 1)
shape_collapse_per_bin = np.zeros(len(bin_edges) - 1)
loc_collapse_per_bin = np.zeros(len(bin_edges) - 1)
scale_collapse_per_bin = np.zeros(len(bin_edges) - 1)
N_collapse_per_bin = np.zeros(len(bin_edges) - 1)
min_collapse_per_bin = np.zeros(len(bin_edges) - 1)
max_collapse_per_bin = np.zeros(len(bin_edges) - 1)
avg_collapse_per_bin = np.zeros(len(bin_edges) - 1)

# Calculate the mean collapse fraction for each bin
for i in range(1, len(bin_edges)):
    # Select the collapse fractions & deltas that fall into the current bin
    in_bin_fcoll = fcoll_nonzero[bin_indices == i]
    in_bin_delta = delta_nonzero[bin_indices == i]
    if len(in_bin_fcoll)>10:
        plt.figure()
        plt.hist(in_bin_fcoll, bins=min(int(len(in_bin_fcoll)/2) +2,200), histtype='step', density=True, label='Data')
        shape, loc, scale = lognorm.fit(in_bin_fcoll)
        x = np.linspace(min(in_bin_fcoll), max(in_bin_fcoll), 1000)
        pdf_fitted = lognorm.pdf(x, shape, loc=loc, scale=scale)
        plt.plot(x, pdf_fitted, 'r-', label='Fitted lognormal', linestyle='dotted', alpha=0.5, label='fitted Log')
        plt.legend()
        plt.xlabel(r'$f_{coll}$')
        plt.title(r'$\delta=$ {:.3f}'.format((bin_edges[i-1]+bin_edges[i])/2))
        plt.semilogx()
        plt.show()
        plt.close()
