# This shoudl be run once you have maneaged to generate the ratios using the provided code

import numpy as np
import pandas as pd


density_files = '/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/dens_fields/DM_count_dens/'

scale_ratios = '/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/subgrid/LMACH/normalise_src/ratios_in_mass_LMACHs.txt' # the ratio of fcoll_implemented/fcoll_CubeP3M
fcoll_path = '/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/subgrid/LMACH/results/'


scale_ratios = np.loadtxt(scale_ratios)
redshifts = scale_ratios[:,0]


box_size = 100 # in comoving Mpc/h
Ngrid=256

pc=  3.086e18 #1 pc in cm
Mpc = 1e6*pc
h = 0.677

Msol =  1.989e33 #grams

Ngrid_volume = (((box_size/h)**3) / (Ngrid**3)) * (Mpc**3) # volume of subgrid in g/cm**3 

m=0
for i in redshifts:
    
    print('Doing z={:.3f}'.format(i))
    fcoll = np.load(fcoll_path+'fcoll_z{:.3f}.npy'.format(i))
    fcoll = fcoll.reshape(256,256,256)
    dens = np.load(density_files+'DM_dens_cgs_{:.3f}_1024_256.npy'.format(i))
    Mass = dens* fcoll * Ngrid_volume / Msol    

    Mass = Mass/scale_ratios[m,1] #scaling the mass to the high res sim

    x, y, z = np.indices(Mass.shape)

    # Step 3: Flatten the arrays
    x_flat = x.flatten()
    y_flat = y.flatten()
    z_flat = z.flatten()
    mass_flat = Mass.flatten()
    
    src = pd.DataFrame({
    'x': x_flat + 1,  # Add 1 to convert from 0-based index to 1-based if necessary
    'y': y_flat + 1,
    'z': z_flat + 1,
    'HMACH': mass_flat  # Assign the flattened mass values to LMACH
    })

    with open('./src_LMACHsubgrid_{:.3f}.txt'.format(i), 'w') as file:
            # Write the number of grid cells
            file.write(f"{len(src)}\n")
            src.to_csv(file, sep=' ', index=False, header=False)

    m+=1
