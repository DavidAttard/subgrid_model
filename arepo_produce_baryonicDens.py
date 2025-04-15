import numpy as np
import glob 
import re

def read_cic_binary(filename):
    """ Reads CIC binary output file """
    
    with open(filename, 'rb') as f:
        # Read header
        trail_byte = np.fromfile(f, dtype=np.int32, count=1)  # Should be 12 for the 3 integers (4 bytes each)
        # print('trail_byte', trail_byte)
        ncells = np.fromfile(f, dtype=np.int32, count=1)[0]
        # print('ncells: ', ncells)
        
        trail_byte = np.fromfile(f, dtype=np.int32, count=2)  # Should be 12 for the 3 integers (4 bytes each)
        # print('trail_byte', trail_byte)
        ntot = np.fromfile(f, dtype=np.int64, count=1)[0]
        # print('ntot: ', ntot)
        
        trail_byte = np.fromfile(f, dtype=np.int32, count=2)
        # print('trail_byte', trail_byte)

        real_boxsize = np.fromfile(f, dtype=np.float32, count=1)[0]
        # print('real_boxsize (Mpc): ', real_boxsize)

        trail_byte = np.fromfile(f, dtype=np.int32, count=2)
        # print('trail_byte', trail_byte)

        real_mpart = np.fromfile(f, dtype=np.float32, count=1)[0]
        # print('real_mpart: ', real_mpart)

        # print(f"ncells: {ncells}, ntot: {ntot}, boxsize: {real_boxsize}, mpart: {real_mpart}")
        
        trail_byte = np.fromfile(f, dtype=np.int32, count=2)
        # print('trail_byte', trail_byte)
        
        # # Read data in slabs (slices along one axis)
        # cic_grid = np.zeros((ncells, ncells, ncells), dtype=np.float32)
        
        # for slab in range(ncells):
        #     cic_grid[:, :, slab] = np.fromfile(f, dtype=np.float32, count=ncells*ncells).reshape((ncells, ncells))

        toto = np.fromfile(f, dtype=np.float64, count=ncells*ncells*ncells)
        cic_grid = toto.reshape((ncells, ncells, ncells), order='F')

    return ncells, ntot, real_boxsize, real_mpart, cic_grid



path_to_cube = '/scratch/da500/HESTIA_multiRealisations_proj/coarsened_density_cic/09_18/'
path_to_snap = '/p/scratch/lgreion/david/runs/zoom_pyc2ray/coarser_densities/dens/'
save_path = '/scratch/da500/HESTIA_multiRealisations_proj/coarsened_density_cic/09_18/DM_density_files/'
realisation = '0918'
AHF_path = '/store/clues/HESTIA/RE_SIMS/4096/DM_ONLY/09_18/AHF_output/'
start_snap = 14
end_snap = 45

Ngrid = 256 # size of the mesh

# Cosmology
h = 0.677
OmegaB = 0.048
OmegaM = 0.318 
box_size = 100 #Mpc/h

pc=  3.086e18 #1 pc in cm
Mpc = 1e6*pc
G_grav = 6.6732e-8

H0 = 100.0*h
H0cgs = H0*1e5/Mpc

for i in range(start_snap,end_snap+1): 

    filename = '{}hestia_{}_{:03d}.BCIC.{}'.format(path_to_cube, realisation, i, Ngrid)
    print(filename)
    ncells, ntot, real_boxsize, real_mpart, cic_grid = read_cic_binary(filename)

    cic_grid = cic_grid/np.sum(cic_grid)

    rho_crit_0 = 3.0*H0cgs*H0cgs/(8.0*np.pi*G_grav)

    # Matter density field in (comoving) g/cm^3
    baryonic_dens = cic_grid*rho_crit_0*OmegaB*(Ngrid)**3 

    # Pattern to match files ending with "_halos"
    AHF_directory = '{}/HESTIA_100Mpc_4096_09_18.{:03d}.*.AHF_halos'.format(AHF_path, i)

    # Get list of files ending with "_halos"
    halos_files = glob.glob(AHF_directory)
    # Print full file paths
    for AHF_file in halos_files:
        print(AHF_file)

    # Regular expression pattern to match the number between 'z' and '.AHF'
    pattern = r'z(\d+\.\d+)'

    # Use re.search to find the match
    match = re.search(pattern, AHF_file)

    # Extract the number from the match
    if match:
        redshift = match.group(1)
        print(redshift)
        z  = float(redshift)
    
    np.save(save_path+'dens_cgs_{:.3f}.npy'.format(z), baryonic_dens)
    print(z)


print('Theoretical value: ',rho_crit_0*OmegaB )
print('Avg value in sim: ', np.average(baryonic_dens))

    

    
