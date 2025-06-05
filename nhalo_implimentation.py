# source /p/project/lgreion/david/hestia_bin/myenv/bin/activate

# Note you have to change the fname file to the file containing the density+Nhalo data
# Also change the mock density field file name to the file containing the density field for 
# which you are generating the mock halo for
# Also change Ncube to the the number of cells in your grid

# replace 24.597 to zred_HESTIA[j]
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams.update({'font.size': 15})

min_resolved_mass = 10**8
frac_diff = []
redshift_diff = []

# Function to ensure each bin has at least min_count values
def ensure_min_bin_count(bin_edges, data, min_count=10):
    # Initialize variables
    new_bin_edges = [bin_edges[0]]  # Start with the first bin edge
    count_in_bin = 0  # Count of elements in the current bin

    for i in range(1, len(bin_edges)):
        # Count elements in the current bin
        count_in_bin += np.sum((data >= bin_edges[i-1]) & (data < bin_edges[i]))

        if count_in_bin >= min_count:  # If bin has enough values
            new_bin_edges.append(bin_edges[i])  # Close the bin
            count_in_bin = 0  # Reset count for the next bin

    # If any remaining data is left in the last bin, include the final edge
    if new_bin_edges[-1] != bin_edges[-1]:
        new_bin_edges.append(bin_edges[-1])

    return np.array(new_bin_edges)

def adaptive_bins_with_initial_empty(data, start, end, bin_width, min_count, max_bin_width):
    """
    Keep initial empty bins with fixed bin width. Once data is found, adaptively widen bins
    up to max_bin_width to meet min_count.
    """
    bin_edges = [start]
    current = start
    found_data = False

    while current < end:
        next_edge = current + bin_width
        count = np.sum((data >= current) & (data < next_edge))

        if not found_data:
            bin_edges.append(next_edge)
            if count > 0:
                found_data = True
                current = next_edge
            else:
                current = next_edge
        else:
            width = bin_width
            while width <= max_bin_width and current + width <= end:
                count = np.sum((data >= current) & (data < current + width))
                if count >= min_count:
                    bin_edges.append(current + width)
                    current += width
                    break
                width += bin_width
            else:
                # Force creation of a bin with `bin_width` if none meet the min_count
                if current + bin_width <= end:
                    bin_edges.append(current + bin_width)
                    current += bin_width
                else:
                    bin_edges.append(end)
                    break

    return np.array(bin_edges)
  


Ncube=14**3

redshifts = np.loadtxt('./redshift_list.txt')
print(redshifts[1:-30:4])
for redshift in redshifts[1:-45:3]:

    # Code to extract infor from the high res simulation 4096^3 coarsened to 256^3
    # fname = f"/p/project/lgreion/david/subgrid_MultiGrid/LMACH/binning/merged_coarsened_halo_num_256_64_z{sz}.dat" 
    # fname = f"/its/research/prace/cubepm_090630_6_1728_6.3Mpc/postprocessing/{sz}halo_num14.bin"
    fname = f"/its/research/prace/cubepm_090630_6_1728_6.3Mpc/postprocessing/{redshift:.03f}halo_num14.bin"
    # fname = f"/p/project/lgreion/david/subgrid_MultiGrid/LMACH/binning/merged_coarsened_256_64_z{sz}.dat"

    # This reads column 0 into dens, column 1 (last) into halo numbers in one pass.
    dens, halo_num1, halo_num2, halo_num3, halo_num4, halo_num5, halo_num6, halo_num7, halo_num8 = np.loadtxt(
        fname, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8), unpack=True
    )

    halo_num = halo_num1 + halo_num2 + halo_num3 + halo_num4 + halo_num5 + halo_num6 + halo_num7 + halo_num8

    # dens, halo_num = np.loadtxt(fname, usecols=(0, 1), unpack=True)

    overdense = dens / dens.mean()
    delta = overdense-1
    delta_nonzero = delta[halo_num!=0]
    halo_num_nonzero = halo_num[halo_num!=0]
    min_halo_num = 0

    mock_density_field = f"/its/research/prace/cubepm_090630_6_1728_6.3Mpc/postprocessing/{redshift:.03f}halo_num14.bin"
    # This reads column 0 into dens, column 1 (last) into halo numbers in one pass.
    dens_mock = np.loadtxt( fname, usecols=(0), unpack=True)

    overdense_mock = dens_mock / dens_mock.mean()
    delta_mock = overdense_mock-1
    overdenisty_hestia = delta_mock
    Nhalo = np.zeros_like(overdenisty_hestia)


    bin_edges = adaptive_bins_with_initial_empty(
        delta_nonzero,
        start=-1,
        end=np.max(delta_nonzero),
        bin_width=0.01,
        min_count=100,
        max_bin_width=0.2
    )

    # Digitize the overdensity values to determine which bin they fall into
    bin_indices = np.digitize(delta_nonzero, bin_edges)
    bin_indices_inc_0 = np.digitize(delta, bin_edges)

    print("bin_indices_inc_0", bin_indices_inc_0)


    # Initialize an array to store the mean collapse fraction for each bin
    delta_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    shape_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    loc_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    scale_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    N_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    min_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    max_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    avg_collapse_per_bin = np.zeros(len(bin_edges) - 1)
    prob_0 = np.zeros(len(bin_edges) - 1)
    tot_Nhalos = np.zeros(len(bin_edges) - 1)
    N_collapse_per_bin_including_zero = np.zeros(len(bin_edges) - 1)
    counts_total_array = np.zeros(len(bin_edges) - 1)
    # Calculate the mean collapse fraction for each bin
    for i in range(1, len(bin_edges)):
        # Select the collapse fractions & deltas that fall into the current bin
        in_bin_halo_num = halo_num_nonzero[bin_indices == i]
        in_bin_halo_num_inc_0 = halo_num[bin_indices_inc_0 == i]
        in_bin_delta = delta_nonzero[bin_indices == i]
        # Compute the mean collapse fraction for this bin
        if len(in_bin_halo_num) > 0:
            # parameter fitting for lognormal
            shape_collapse_per_bin[i-1], loc_collapse_per_bin[i-1], scale_collapse_per_bin[i-1] = lognorm.fit(in_bin_halo_num)  
            # avg delta value of bins
            delta_collapse_per_bin[i-1] = np.mean(in_bin_delta)
            
            # number of cells in the bin with a non-zero nhalo
            N_collapse_per_bin[i-1] = len(in_bin_halo_num)

            N_collapse_per_bin_including_zero[i-1] = len(in_bin_halo_num_inc_0)

            # min and max values of number of halos that are found within all cells of a particular bin
            min_collapse_per_bin[i-1] = min(in_bin_halo_num)
            max_collapse_per_bin[i-1] = max(in_bin_halo_num)

            # avg value of nhalos found in all cells within the bin 
            avg_collapse_per_bin[i-1] = np.mean(in_bin_halo_num)


            counts_total = N_collapse_per_bin_including_zero[i-1]
            counts_non_zero = N_collapse_per_bin[i-1]
            prob_0[i-1] = (counts_total - counts_non_zero) / counts_total
            counts_total_array[i-1] = counts_total 
            
            tot_Nhalos[i-1] = np.sum(in_bin_halo_num)
            plt.figure(figsize=(6, 4))
            x_vals = np.linspace(min(in_bin_halo_num), max(in_bin_halo_num), 1000)

            mu = np.mean(np.log(in_bin_halo_num))
            sigma = np.std(np.log(in_bin_halo_num))
            shape_collapse_per_bin[i-1] = sigma
            scale_collapse_per_bin[i-1] = np.exp(mu)
            loc = loc_collapse_per_bin
            # log_pdf_loc0 = lognorm.pdf(x_vals, shape, loc=0, scale=scale)

            fitted_pdf = lognorm.pdf(x_vals, shape_collapse_per_bin[i-1],
                                        loc=loc_collapse_per_bin[i-1],
                                        scale=scale_collapse_per_bin[i-1])
            
            
            # Normalize histogram to match PDF scale
            plt.hist(in_bin_halo_num, bins=20, density = True, alpha=1, color='black', label='Data', histtype='step')
            plt.plot(x_vals, fitted_pdf, 'r-', lw=2, label='Lognorm fit')
            plt.semilogy()
            bin_range = f"[{bin_edges[i-1]:.2f}, {bin_edges[i]:.2f})"
            plt.title(rf'Bin {i-1}: $\delta$ in [{bin_edges[i-1]:.2f}, {bin_edges[i]:.2f}), nonzero cells in bin {N_collapse_per_bin[i-1]}', fontsize ='12')
            plt.xlabel('Number of Halos')
            plt.ylabel('Probability Density')
            plt.legend()
            plt.tight_layout()
            if N_collapse_per_bin[i-1] > 10:
                # plt.savefig(f'lognorm_fits_dynamic/z_bin_{i-1}.png')
                plt.close()

    # Find the index of the first non-zero element
    first_nonzero_idx = np.argmax(prob_0 != 0)
    # Set all values before that index to 1
    prob_0[:first_nonzero_idx] = 1

    data = np.column_stack((delta_collapse_per_bin, shape_collapse_per_bin, loc_collapse_per_bin, scale_collapse_per_bin, N_collapse_per_bin, min_collapse_per_bin, max_collapse_per_bin,avg_collapse_per_bin, prob_0, tot_Nhalos, counts_total_array))

    fout=open('./logfiles/logfile_fits_{:.3f}.txt'.format(redshift), "w")
    np.savetxt(fout,data,fmt="%.3e",header='delta, shape, loc, scale, N halo in bins, minimum halo_num, max halo_num, mean halo_num, prob_0, tot_Nhalos, counts_total_inbin_inc_zero')

    plt.figure()
    plt.plot(np.log10(delta_collapse_per_bin[counts_total_array > 0] + 1), prob_0[counts_total_array > 0])

    plt.xlabel('$\log_{10}(\delta + 1)$')
    plt.ylabel(r'Probability $N_{halo, 5:8}=0$')
    plt.savefig('./empty_halo_num/nhalo_prob_z{:.3f}.png'.format(redshift))

    # Implementing the Nhalos in the mock halo catalogue 
    bin_indices = np.digitize(overdenisty_hestia, bin_edges)
    unique_bin_indices = np.unique(bin_indices)

    for bin in unique_bin_indices:
        bin_index = bin - 1  # Correct the offset

        if bin_index < 0 or bin_index >= len(prob_0):
            continue  # Safely skip invalid indices
        indices = np.where(bin_indices == bin)[0]

        # print('Total cells in bin in CubeP3M: ', N_collapse_per_bin_including_zero[bin-1])
        # print('Total cells in bin in HESTIA: ', len(indices))

        if N_collapse_per_bin_including_zero[bin-1] == 0:
            continue
            
        # The error fraction calculates how the difference in the number of cells that exists in a given bin     
        error_fraction = (N_collapse_per_bin_including_zero[bin-1] - len(indices)) / N_collapse_per_bin_including_zero[bin-1]
        expected_nhalo_implemented  =  tot_Nhalos[bin-1]*(1-error_fraction)

        # print('error_fraction = ', error_fraction)
        # print('expected_nhalo_implemented = ', expected_nhalo_implemented)
        # print('Nhalos in current bin in CubeP3m: ', tot_Nhalos[bin-1])

        # Remove those cells which will be given a 0 value
        num_zero = int(prob_0[bin-1] * len(indices))
        to_zero = np.random.choice(indices, size=num_zero, replace=False)
        Nhalo[to_zero] = 0
        
        # Contains everything in indices except the values that were in to_zero
        indices = np.setdiff1d(indices, to_zero)


        # If statistic is generated using more tha 10 data points       
        if N_collapse_per_bin[bin-1]>10:
            shape = shape_collapse_per_bin[bin-1]
            loc = loc_collapse_per_bin[bin-1]
            scale = scale_collapse_per_bin[bin-1]

            if (scale > 0) & (shape>0): 
                #generate random Nhalo values
                Nhalo_possible_values = np.rint(lognorm.rvs(shape, loc=0, scale=scale, size=len(overdenisty_hestia[indices])*200))
            
                if len(Nhalo_possible_values[(Nhalo_possible_values>0) & (Nhalo_possible_values <  max_collapse_per_bin[bin-1])]) > len(overdenisty_hestia[indices]):
                    counter = 0

                    while np.absolute(np.sum(Nhalo[indices])- expected_nhalo_implemented)/expected_nhalo_implemented > 0.02: 
                        counter +=1
                        # print('max_collapse_per_bin[bin-1]: ', max_collapse_per_bin[bin-1])
                        Nhalo_possible_values = Nhalo_possible_values[(Nhalo_possible_values>0) & (Nhalo_possible_values <  max_collapse_per_bin[bin-1])]
                        # print('max(Nhalo_possible_values): ', np.max(Nhalo_possible_values))

                        Nhalo[indices] = np.random.choice(Nhalo_possible_values, size=len(overdenisty_hestia[indices]), replace=False)
                        # print('Difference in Nhalos implemented in the bin and expected implemented: ', np.absolute(np.sum(Nhalo[indices])- expected_nhalo_implemented))
                        if counter > 50:
                            break

                else:
                # Replace any Nhalo values that are larger than what is allowed to be the average value 
                    avg_indices = np.where( (Nhalo_possible_values<0) & (Nhalo_possible_values >  max_collapse_per_bin[bin-1]))[0]
                    Nhalo_possible_values[avg_indices] = np.rint(avg_collapse_per_bin[bin-1])

                    # Select randomly the required number of Nhalos to be filled
                    Nhalo[indices] = np.random.choice(Nhalo_possible_values, size=len(overdenisty_hestia[indices]), replace=False)
                

            else: 
                # print('Using avg: ')
                Nhalo[indices] = np.rint(avg_collapse_per_bin[bin-1])

        else: 
            # print('Using avg...too little stats:')
            Nhalo[indices] = np.rint(avg_collapse_per_bin[bin-1])

        # print('Nhalos implemented in HESTIA in current bin: ', np.sum(Nhalo[indices]))
        # print('Cells to fill in  HESTIA: ', len(Nhalo[indices]))
        # print('Cells filled in CubeP3M: ', N_collapse_per_bin[bin-1])

    print('Total N halos implemented (HESTIA): ', np.sum(Nhalo)) 
    grid1 = 256**3   
    print('Total N halos in High-res (CubeP3M): ', np.sum(halo_num[:grid1]))  
    print('Fractional diff: ', np.sum(Nhalo)/np.sum(halo_num[:grid1]))
    frac_diff = np.append(frac_diff, np.sum(Nhalo)/np.sum(halo_num[:grid1]))
    redshift_diff = np.append(redshift_diff, redshift)

    # Log-transform the data and filter for the first heatmap
    x1 = np.log10(overdenisty_hestia[Nhalo > 0]+1)
    y1 = np.log10(Nhalo[Nhalo > 0])

    # Create a 2D histogram for the first heatmap
    heatmap1, xedges1, yedges1 = np.histogram2d(x1, y1, bins=100)


    # Log-transform the data and filter for the second heatmap
    x2 = np.log10(overdense[halo_num>0])
    y2 = np.log10(halo_num[halo_num>0])

    # Create a 2D histogram for the second heatmap
    heatmap2, xedges2, yedges2 = np.histogram2d(x2, y2, bins=100)

    combined_max = max(heatmap1.max(), heatmap2.max())
    combined_min = max(1, min(heatmap1[heatmap1 > 0].min(), heatmap2[heatmap2 > 0].min()))
    norm = LogNorm(vmin=combined_min, vmax=combined_max)

    # Set up the figure without sharey
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))  # No sharey here

    # Plot heatmap 1
    im1 = ax1.imshow(heatmap1.T, origin='lower', cmap='viridis', aspect='auto',
                    extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]], norm=norm)
    ax1.set_title('Implemented Mock catalogue')
    ax1.set_xlabel('$\log_{10}(\delta + 1)$')
    ax1.set_ylabel(r'log$_{10}$(N$_{halo,5:8})$')

    # Plot heatmap 2
    im2 = ax2.imshow(heatmap2.T, origin='lower', cmap='viridis', aspect='auto',
                    extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]], norm=norm)
    ax2.set_title('High Res simulation')
    ax2.set_xlabel('$\log_{10}(\delta + 1)$')
    # ax2.set_ylabel(r'N$_{halo}$')

    # Manually align y-axis limits
    shared_ymin = min(yedges1[0], yedges2[0])
    shared_ymax = max(yedges1[-1], yedges2[-1])
    ax1.set_ylim(shared_ymin, shared_ymax)
    ax2.set_ylim(shared_ymin, shared_ymax)

    shared_xmin = min(xedges1[0], xedges2[0])
    shared_xmax = max(xedges1[-1], xedges2[-1])
    ax1.set_xlim(shared_xmin, shared_xmax)
    ax2.set_xlim(shared_xmin, shared_xmax)

    # Shared colorbar
    fig.colorbar(im2, ax=[ax1, ax2], label='Density')


    plt.savefig('./results/scatterImp_N5_8_14_z{:.3f}.png'.format(redshift))
    np.save('./results/halo_num_z{:.3f}'.format(redshift), Nhalo)


    
fractional_diff = np.column_stack((frac_diff, redshift_diff))
plt.figure()
plt.plot(redshift_diff, frac_diff)
plt.xlabel(r'$z$')
plt.ylabel(r'$N_{\text{total, mock}}/ N_{\text{total, CubeP3M}}$')
plt.tight_layout()
plt.savefig('fractional_diff.png')
fout=open('fractional_diff.txt', "w")
np.savetxt(fout,fractional_diff,fmt="%.3e")
