import os
import sys
import scipy.stats 
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline, interp1d


min_resolved_mass = 10**9

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


Ncube=256**3

zlist = np.loadtxt('/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/subgrid/HMACH_9_9.4/redshift_matching/matched_z.txt')

zred_HESTIA = zlist[:,0]
zred_highRes = zlist[:,1]
# zred=7.570
print(zred_HESTIA)
# for j in range(7,len(zred_HESTIA)):
for j in range(28,len(zred_HESTIA)):

    print('Doing z={:.3f}'.format(zred_HESTIA[j]))
    # Code to extract infor from the high res simulation 4096^3 coarsened to 256^3
    sz = '%.3f' %zred_highRes[j]
    txt_file = open('/p/project/lgreion/david/HESTIA_multirealizations/subgrid_modelling/HMACH_binning/binning/'+sz+'halo_mass256.bin')
    fcoll = []
    dens = []

    for line in txt_file:
        # we can process file line by line here, for simplicity I am taking count of lines
        elements = line.split()
        dens_value = float(elements[0])
        dens.append(dens_value)
        fcoll_value = float(elements[-1])
        fcoll.append(fcoll_value)
        
    txt_file.close()

    # The density values are stored in the 1st entry of every two lines
    #dens = dens[0::2]
    overdense = np.array(dens/np.average(dens))

    fcoll = np.array(fcoll)# Code to extract infor from the high res simulation 4096^3 coarsened to 256^3
    min_fcoll = min_resolved_mass/np.array(dens)

    # Loading the HESTIA data
    density_data=np.load('/p/project/lgreion/david/HESTIA_multirealizations/runs/Full_box/37_11/dens_fields/DM_count_dens/DM_dens_cgs_{:.3f}_1024_256.npy'.format(zred_HESTIA[j])).flatten()
    print('Loaded: DM_dens_cgs_{:.3f}_1024_256.npy'.format(zred_HESTIA[j]))

    # Find the average density
    rho_avg = np.sum(density_data)/Ncube
    print('rho_avg: ', rho_avg)
    density_field = density_data/rho_avg

    # Assume these are your input arrays
    delta = overdense-1
    delta_nonzero = delta[fcoll!=0]
    fcoll_nonzero = fcoll[fcoll!=0]

    # initial_bin_edges = np.arange(-1, 15 + 0.1, 0.1)
    # initial_bin_edges = np.append(initial_bin_edges, [20,25,30,35,40,50,60,70,80,90,100,120,140,160,180,200,240,260,300])  # Wider bins for larger values

    # Initial bin edges from -1 to 15 with width 0.1
    initial_bin_edges = np.arange(-1, 15 + 0.1, 0.1)
    bins = list(initial_bin_edges)

    # Define variables for binning delta > 15
    start_bin = 15
    current_bin_width = 0.1

    # Sort delta_nonzero values above 15 for efficient processing
    delta_above_15 = np.sort(delta_nonzero[delta_nonzero > start_bin])

    # Loop to create bins for delta > 15
    while len(delta_above_15) > 10:
        # Calculate the end of the current bin
        next_bin_edge = start_bin + current_bin_width
        
        # Count the entries in the current bin
        count_in_bin = np.sum((delta_above_15 >= start_bin) & (delta_above_15 < next_bin_edge))
        
        if count_in_bin >= 10:
            # If there are at least 10 entries, accept the current bin width
            bins.append(next_bin_edge)
            # Remove entries that fall into the current bin
            delta_above_15 = delta_above_15[delta_above_15 >= next_bin_edge]
            # Reset start_bin and bin width
            start_bin = next_bin_edge
            current_bin_width = 5  # Reset bin width to 0.1 for next bin
            print(bins)
            print(len(delta_above_15))
        else:
            # If fewer than 10 entries, increase the bin width and try again
            current_bin_width += 5

    # Convert bins to numpy array for easy use
    bin_edges = np.array(bins)


    # Digitize the data into initial bins
    bin_indices = np.digitize(delta_nonzero, bins=bin_edges) - 1

    print(bin_edges)
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
        # Compute the mean collapse fraction for this bin
        if len(in_bin_fcoll) > 0:
            
            shape_collapse_per_bin[i-1], loc_collapse_per_bin[i-1], scale_collapse_per_bin[i-1] = lognorm.fit(in_bin_fcoll)  
            delta_collapse_per_bin[i-1] = np.mean(in_bin_delta)
            N_collapse_per_bin[i-1] = len(in_bin_fcoll)
            min_collapse_per_bin[i-1] = min(in_bin_fcoll)
            max_collapse_per_bin[i-1] = max(in_bin_fcoll)
            avg_collapse_per_bin[i-1] = np.mean(in_bin_fcoll)
        
    data = np.column_stack((delta_collapse_per_bin, shape_collapse_per_bin, loc_collapse_per_bin, scale_collapse_per_bin, N_collapse_per_bin, min_collapse_per_bin, max_collapse_per_bin,avg_collapse_per_bin))

    if len(data) == 0:        

        # Digitize the data into initial bins
        bin_indices = np.digitize(delta_nonzero, bins=initial_bin_edges) - 1


        # Call the function to adjust the bin edges based on the data
        bin_edges = ensure_min_bin_count(initial_bin_edges, delta_nonzero, min_count=10)

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

            # Compute the mean collapse fraction for this bin
            if len(in_bin_fcoll) > 0:

                shape_collapse_per_bin[i-1], loc_collapse_per_bin[i-1], scale_collapse_per_bin[i-1] = lognorm.fit(in_bin_fcoll)  
                delta_collapse_per_bin[i-1] = np.mean(in_bin_delta)
                N_collapse_per_bin[i-1] = len(in_bin_fcoll)
                min_collapse_per_bin[i-1] = min(in_bin_fcoll)
                max_collapse_per_bin[i-1] = max(in_bin_fcoll)
                avg_collapse_per_bin[i-1] = np.mean(in_bin_fcoll)
        
    data = np.column_stack((delta_collapse_per_bin, shape_collapse_per_bin, loc_collapse_per_bin, scale_collapse_per_bin, N_collapse_per_bin, min_collapse_per_bin, max_collapse_per_bin,avg_collapse_per_bin))

        
    print()
    fout=open('./logfiles/logfile_fits_{:.3f}.txt'.format(zred_HESTIA[j]), "w")
    np.savetxt(fout,data,fmt="%.3e",header='delta, shape, loc, scale, N fcoll in bins, minimum fcoll, max fcoll, mean fcoll')


    # filename = './new_logfile.txt'
    # #filename = '7.570_lognormdistr_114Mpcbox256'
    # filein = open(filename, 'r')
    # pdf = np.genfromtxt(filein, float)

    pdf = data

    pdf=pdf[pdf[:,0]!=0]
    pdf0 = pdf[1:][np.argsort(pdf[1:,0])] ## Make sure delta is increasing

    # First interpolate ln(mu) based on 1+delta, get rid of bins with 0
    overdensity_shape = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,1][pdf0[:,4]>10],k=1)
    shape = overdensity_shape(density_field)

    overdense_loc = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,2][pdf0[:,4]>10],k=1)
    loc = overdense_loc(density_field)

    overdense_scale = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,3][pdf0[:,4]>10],k=1)
    scale = overdense_scale(density_field)
    
    delta = delta[:10000]
    min_fcoll = min_fcoll[:10000]
    delta_sorted_indices = np.argsort(delta)  # Get the indices that would sort delta
    delta_sorted = delta[delta_sorted_indices]  # Sort delta
    min_fcoll_sorted = min_fcoll[delta_sorted_indices]  # Sort min_fcoll to match delta's order

    # Downsample the sorted delta and min_fcoll for the interpolation
    delta_sampled = delta_sorted + 1  # Add 1 if necessary for your specific data
    min_fcoll_sampled = min_fcoll_sorted

    # Remove duplicates to ensure delta_sampled is strictly increasing
    delta_unique, unique_indices = np.unique(delta_sampled, return_index=True)
    min_fcoll_unique = min_fcoll_sampled[unique_indices]

    # Create the interpolating spline
    delta_min = InterpolatedUnivariateSpline(delta_unique, min_fcoll_unique, k=1)
    # delta_min = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>20]+1,pdf0[:,5][pdf0[:,4]>20], k=1)
    min_points = delta_min(density_field)

    delta_max = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,6][pdf0[:,4]>10], k=1)
    max_points = delta_max(density_field)

    delta_mean = InterpolatedUnivariateSpline(pdf0[:,0]+1,pdf0[:,7], k=1)
    mean_points = delta_mean(density_field)

    fig = plt.figure(figsize=(12,14))
    ax1 = fig.add_subplot(611)
    ax2 = fig.add_subplot(612)
    ax3 = fig.add_subplot(613)
    ax4 = fig.add_subplot(614)
    ax5 = fig.add_subplot(615)
    ax6 = fig.add_subplot(616)


    # Downsample every 500th point for plotting
    density_field_sampled = density_field[::500]
    shape_sampled = shape[::500]
    loc_sampled = loc[::500]
    scale_sampled = scale[::500]
    min_points_sampled = min_points[::500]
    max_points_sampled = max_points[::500]
    mean_points_sampled = mean_points[::500]

    # Plot each subplot using the sampled data

    # ax1
    print('Max log_10(delta+1) in cubeP3M = ', np.max(np.log10(pdf0[:, 0] + 1)))
    ax1.plot(np.log10(pdf0[:, 0] + 1), pdf0[:, 1], "k-", label='Data')
    ax1.scatter(np.log10(density_field_sampled), shape_sampled, color='red', label='interpolation')

    # ax2
    ax2.plot(np.log10(pdf0[:, 0] + 1), pdf0[:, 2], "k-")
    ax2.scatter(np.log10(density_field_sampled), loc_sampled, color='red')

    # ax3
    ax3.plot(np.log10(pdf0[:, 0] + 1), pdf0[:, 3], "k-")
    ax3.scatter(np.log10(density_field_sampled), scale_sampled, color='red')

    # ax4
    ax4.plot(np.log10(density_field_sampled), min_points_sampled, "r.")
    ax4.plot(np.log10(pdf[:, 0] + 1), pdf[:, 5], "k-", lw=2.5)

    # ax5
    ax5.plot(np.log10(density_field_sampled), max_points_sampled, "r.")
    ax5.plot(np.log10(pdf[:, 0] + 1), pdf[:, 6], "k-", lw=2.5)

    # ax6
    ax6.plot(np.log10(density_field_sampled), mean_points_sampled, "r.")
    ax6.plot(np.log10(pdf[:, 0] + 1), pdf[:, 7], "k-", lw=2.5)



    ax1.set_title("shape")
    ax2.set_title("loc")
    ax3.set_title("scale")
    ax4.set_title("min f_coll")
    ax5.set_title("max f_coll")
    ax6.set_title("mean f_coll")


    ax6.set_xlabel(r'log$_{10}(1+\delta)$', fontsize=15)
    ax1.legend()
    plt.tight_layout()
    min_points[np.isnan(min_points)] = 0
    plt.savefig('./logfiles/interpolation_z{:.3f}.png'.format(zred_HESTIA[j]))

    print('Log-file generated')

    Nhalo = np.zeros_like(density_field)
    # Create a mask for valid scale values (scale > 0)
    mask = (scale > 0) & (shape>0)

    # Apply the lognormal sampling only where the mask is True
    if np.any(mask):
        Nhalo[mask] = lognorm.rvs(shape[mask], loc=loc[mask], scale=scale[mask])

    complex_condition = (scale > 0) & (shape>0) & (min_points<max_points)  & ((Nhalo>max_points) | (Nhalo<min_points))
    complex_condition_index = np.where(complex_condition==True)
    complex_condition_index= np.array(complex_condition_index[0])

    simple_mask = (scale > 0) & (shape>0) & (Nhalo<min_points) & ~complex_condition
    simple_mask_index = np.where(simple_mask==True)
    simple_mask_index= np.array(simple_mask_index[0])

    i =0
    for k in complex_condition_index:
        while (Nhalo[k]<min_points[k]) or (Nhalo[k]>max_points[k]):
            possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=1000) 
            possible_values = possible_values[(possible_values>min_points[k]) & (possible_values<max_points[k])]
            if len(possible_values)==0:
                possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=1000) 
                possible_values = possible_values[(possible_values>min_points[k]) & (possible_values<max_points[k])]
                if len(possible_values)==0:
                    Nhalo[k] = mean_points[k]
                    break
            
            Nhalo[k] = np.random.choice(possible_values)
            
    
    counter =0        
    for k in simple_mask_index:
        while (Nhalo[k]<min_points[k]):
            possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=1000) 
            possible_values = possible_values[(possible_values>min_points[k])]
            if len(possible_values)==0:
                possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=1000) 
                possible_values = possible_values[(possible_values>min_points[k])]
                if len(possible_values)==0:
                    Nhalo[k] = mean_points[k]
                    break
            if len(possible_values)!=0:
                Nhalo[k] = np.random.choice(possible_values)

        # Nhalo[k] = np.random.choice(possible_values)

    if np.max(Nhalo)>1:
        indices = np.where(Nhalo > 1)[0]
        print(indices)
        if indices.size > 0:
            for m in indices:  # Iterate directly over the indices
                Nhalo[m] = mean_points[m]

    # Due to the lack of statistics in the outer bins we want to populate these with the avg f_coll values
    max_overdensity = np.argmax(density_field)
    print('max_overdensity in HESTIA (log10(delta+1)): ', np.log10(density_field[max_overdensity]))
    print('Nhalo at max_overdens: ', Nhalo[max_overdensity])
    if Nhalo[max_overdensity] == 0:
        non_empty_density = density_field[Nhalo != 0]
        largest_density_non_empty = np.max(non_empty_density)
        
        density_field_populate = np.where(density_field>largest_density_non_empty)[0]
        print(density_field_populate)
        print(len(density_field_populate))
        for k in density_field_populate:
            Nhalo[k]=mean_points[k]
    print('Nhalo at max_overdens: ', Nhalo[max_overdensity])


    indices = np.where(Nhalo < min_points)[0]
    if len(indices)>0:
        for k in indices:
            Nhalo[k] = mean_points[k]
            if Nhalo[k]<min_points[k]:
                Nhalo[k]=min_points[k]

    print('implementing the fcoll=0')
    
    binning=np.array(bins)
    prob_zero_fcoll = {"%0.2f" %bin:[] for bin in binning} #probability of having a zero coll fraction for a given delta


    for binindex in range(len(binning)-1):
        fcoll_current_bin = fcoll[((overdense-1)>binning[binindex]) & ((overdense-1)<binning[binindex+1])]
        nhalo_zero_fcoll = fcoll_current_bin[fcoll_current_bin==0]
        if len(fcoll_current_bin)==0:
            prob_zero_fcoll["%0.2f" %binning[binindex]].append(1)
        else:
            prob_zero_fcoll["%0.2f" %binning[binindex]].append(len(nhalo_zero_fcoll)/len(fcoll_current_bin))

    # Assuming binning and prob_zero_fcoll are your variables
    x_values = []
    y_values = []

    # Extracting x and y values
    for binindex in range(len(binning)):
        bin_key = "%0.2f" % binning[binindex]
        if bin_key in prob_zero_fcoll:
            if len(prob_zero_fcoll[bin_key])==1:
                x_values.append(binning[binindex])
                y_values.append(prob_zero_fcoll[bin_key][0])
    
    plt.figure()
    plt.plot(np.log10(np.array(x_values)+1), np.array(y_values))
    plt.xlabel('$\log_{10}(\delta + 1)$')
    plt.ylabel(r'Probability $f_{coll}=0$')
    plt.savefig('./empty_fcoll/zero_fcoll_prob_z{:.3f}.png'.format(zred_HESTIA[j]))

    overdens_hestia = density_field-1

    Nhalo_new = np.copy(Nhalo)

    for binindex in range(len(binning)-1):
        bin_key = "%0.2f" % binning[binindex]
        Nhalo_current_bin = Nhalo_new[(overdens_hestia>binning[binindex]) & (overdens_hestia<binning[binindex+1])]
        
        if len(prob_zero_fcoll[bin_key])==1:
            prob_zero = prob_zero_fcoll[bin_key][0]

            # Calculate the number of entries to set to zero
            num_to_zero = int(prob_zero * len(Nhalo_current_bin))

            # Randomly select indices to set to zero
            indices_to_zero = np.random.choice(len(Nhalo_current_bin), size=num_to_zero, replace=False)

            # Set the selected entries to zero
            Nhalo_current_bin[indices_to_zero] = 0

            Nhalo_new[(overdens_hestia>binning[binindex]) & (overdens_hestia<binning[binindex+1])] = Nhalo_current_bin
    print('Nhalo at max_overdens: ', Nhalo_new[max_overdensity])

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    # Log-transform the data and filter for the first heatmap
    x1 = np.log10(density_field[Nhalo_new > 0])
    y1 = np.log10(Nhalo_new[Nhalo_new > 0])

    # Create a 2D histogram for the first heatmap
    heatmap1, xedges1, yedges1 = np.histogram2d(x1, y1, bins=100)


    # Log-transform the data and filter for the second heatmap
    x2 = np.log10(overdense[fcoll>0])
    y2 = np.log10(fcoll[fcoll>0])

    # Create a 2D histogram for the second heatmap
    heatmap2, xedges2, yedges2 = np.histogram2d(x2, y2, bins=100)

    # Create a figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot the first heatmap
    im1 = ax1.imshow(heatmap1.T, origin='lower', cmap='seismic', aspect='auto',
                    extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]], norm='log')



    ax1.set_title(r'HESTIA ${}^3$'.format(256))
    ax1.set_xlabel('$\log_{10}(\delta + 1)$')
    ax1.set_ylabel(r'$\log_{10}(f_{coll})$')
    fig.colorbar(im1, ax=ax1, label='Density')

    # Plot the second heatmap
    im2 = ax2.imshow(heatmap2.T, origin='lower', cmap='seismic', aspect='auto',
                    extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]], norm='log')

    # contours = ax1.contour(heatmap2.T, levels=5, origin='lower', 
    #                     extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]], 
    #                     colors='black', label='High-res sim contours', alpha =0.5)

    # # Correct the extent for ax2 to match heatmap2's xedges2 and yedges2
    # contours = ax2.contour(heatmap2.T, levels=5, origin='lower', 
    #                     extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]], 
    #                     colors='black', label='High-res sim contours', alpha =0.5)

    # contour_proxy1 = plt.Line2D([0], [0], color='black', label='High-res sim contours')

    # # Add legend with the contour line label
    # ax1.legend(handles=[contour_proxy1])
    # ax2.legend(handles=[contour_proxy1])


    ax2.set_title('High Res simulation')
    ax1.set_xlabel('$\log_{10}(\delta + 1)$')
    ax2.set_xlabel('$\log_{10}(\delta + 1)$')

    ax1.set_ylabel(r'$\log_{10}(f_{coll})$')
    ax2.set_ylabel(r'$\log_{10}(f_{coll})$')

    fig.colorbar(im2, ax=ax2, label='Density')


    ax1.set_ylim(np.log10(min(float(np.min(Nhalo_new[Nhalo_new > 0])), float(np.min(fcoll[fcoll>0])))), np.log10(max(float(np.max(Nhalo_new[Nhalo_new > 0])), float(np.max(fcoll[fcoll>0])))))
    ax2.set_ylim(np.log10(min(float(np.min(Nhalo_new[Nhalo_new > 0])), float(np.min(fcoll[fcoll>0])))), np.log10(max(float(np.max(Nhalo_new[Nhalo_new > 0])), float(np.max(fcoll[fcoll>0])))))

    ax1.set_xlim(np.log10(min(float(np.min(density_field[Nhalo_new > 0])), float(np.min(overdense[fcoll>0])))), np.log10(max(float(np.max(density_field[Nhalo_new > 0])), float(np.max(overdense[fcoll>0])))))
    ax2.set_xlim(np.log10(min(float(np.min(density_field[Nhalo_new > 0])), float(np.min(overdense[fcoll>0])))), np.log10(max(float(np.max(density_field[Nhalo_new > 0])), float(np.max(overdense[fcoll>0])))))
    fig.suptitle(r'Implemeted $f_{coll}$ for halos $10^9 <M_{halo}[M_{\odot}]< 10^{9.4}$ on $256^3$ grid', fontsize=16)

    plt.savefig('./results/scatterImp_M_9to9_4_1024_256_z{:.3f}.png'.format(zred_HESTIA[j]))
    np.save('./results/fcoll_z{:.3f}'.format(zred_HESTIA[j]), Nhalo_new)


    
    
