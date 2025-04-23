# Calculating the total fcoll in CubeP3M sim 

# source /p/project/lgreion/david/hestia_bin/myenv/bin/activate

import os
import sys
import scipy.stats 
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline, interp1d


min_resolved_mass = 10**(9)
fcoll_path = '/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/subgrid/LMACH/results/'
density_files = '/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/dens_fields/DM_count_dens/'

ratios = []
redshift_list= []
Ncube=256**3

zlist = np.loadtxt('/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/subgrid/LMACH/redshift_matching/matched_z.txt')

zred_HESTIA = zlist[:,0]
zred_highRes = zlist[:,1]
# zred=7.570
print(zred_HESTIA)
for j in range(2,len(zred_HESTIA)):
    print('Doing z={:.3f}'.format(zred_HESTIA[j]))
    # Code to extract infor from the high res simulation 4096^3 coarsened to 256^3
    sz = '%.3f' %zred_highRes[j]
    fcoll = []
    dens = []

    file_path = f'/p/project/lgreion/david/HESTIA_multirealizations/subgrid_modelling/LMACH_binning/binning/{sz}halo_mass256.bin'

    # Load data as a numpy array
    data = np.loadtxt(file_path)

    # Extract density values and fcoll values
    dens = data[:, 0].astype(float)                    # First column as density
    fcoll = np.sum(data[:, 1:-1].astype(float), axis=1) 

    print('Done reading')

    fcoll = np.array(fcoll)
    total_fcoll = np.sum(fcoll)
    print('Total fcoll: ', total_fcoll)
    print('log_10: ', np.log10(total_fcoll))
    tot_mass = np.sum(np.array(dens))
    print('Total mass in sim: ', tot_mass)
    print('Total collapse fraction at z={} is {}'.format(sz, total_fcoll/tot_mass))

    print('Results for subgrid implemented scatter:')
    fcoll = np.load(fcoll_path+'fcoll_z{:.3f}.npy'.format(zred_HESTIA[j]))
    fcoll = fcoll.reshape(256,256,256)
    dens = np.load(density_files+'DM_dens_cgs_{:.3f}_1024_256.npy'.format(zred_HESTIA[j]))
    box_size = 100 # in comoving Mpc/h
    Ngrid=256

    pc=  3.086e18 #1 pc in cm
    Mpc = 1e6*pc
    h = 0.677

    Msol =  1.989e33 #grams
    Ngrid_volume = (((box_size/h)**3) / (Ngrid**3)) * (Mpc**3) # volume of subgrid in g/cm**3 

    Mass = dens* fcoll * Ngrid_volume / Msol  
    total_implemented_mass =  np.sum(Mass)
    print('Total implemented fcoll: ', total_implemented_mass)  
    print('log_10: ', np.log10(total_implemented_mass))

    print('Ratio fcoll_impl: hestia/cubep3m: ', total_implemented_mass/total_fcoll)
    ratios = np.append(ratios, total_implemented_mass/total_fcoll)
    redshift_list = np.append(redshift_list,zred_HESTIA[j])

plt.plot(redshift_list, ratios)
plt.xlabel('z')
plt.ylabel('(total_implemented_mass)/(tot LMACH mass in CubeP3M)')
plt.title('Ratio in the implemented mass (LMACHs)')
plt.savefig('ratios_in_mass_LMACHS.png')

data_to_save = np.column_stack((redshift_list, ratios))  # Stack columns together

# Save to a file, with each row having a redshift and ratio
np.savetxt('ratios_in_mass_LMACHs.txt', data_to_save, header='Redshift Ratios', fmt='%.6f', delimiter='\t')
