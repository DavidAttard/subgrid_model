import os
import sys
import scipy.stats 
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline, interp1d


Ncube=256**3

zlist = np.loadtxt('./redshift_matching/matched_z.txt')

zred_HESTIA = zlist[:,0]
zred_highRes = zlist[:,1]
# zred=7.570

for j in range(len(zred_HESTIA)):
    print('Doing z={:.3f}'.format(zred_HESTIA[j]))
    # Code to extract infor from the high res simulation 4096^3 coarsened to 256^3
    sz = '%.3f' %zred_highRes[j]
    txt_file = open('./binning/'+sz+'halo_mass256.bin')
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


    # Loading the HESTIA data
    density_data=np.load('./HESTIA_256/DM_density/dens_cgs_{:.3f}.npy'.format(zred_HESTIA[j])).flatten()
    print('Loaded: dens_cgs_{:.3f}.npy'.format(zred_HESTIA[j]))

    # Find the average density
    rho_avg = np.sum(density_data)/Ncube
    print('rho_avg: ', rho_avg)
    density_field = density_data/rho_avg

    # Assume these are your input arrays
    delta = overdense-1
    delta_nonzero = delta[fcoll!=0]
    fcoll_nonzero = fcoll[fcoll!=0]

    # Since we are dealing with log-normal distributions we are interested in log(x)
    log_fcoll = np.log(fcoll_nonzero)

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
        # Assume these are your input arrays
        delta = overdense-1
        delta_nonzero = delta[fcoll!=0]
        fcoll_nonzero = fcoll[fcoll!=0]

        # Since we are dealing with log-normal distributions we are interested in log(x)
        log_fcoll = np.log(fcoll_nonzero)

        # Define the bins for overdensities, starting from -1 to 10 with a bin width of 0.1
        bin_edges = np.arange(-1, 10 + 0.1, 0.1)
        # For Larger overdensities we have wider bins for sufficient statistics
        bin_edges = np.append(bin_edges, 15)
        bin_edges = np.append(bin_edges, 30)
        bin_edges = np.append(bin_edges, 60)
        bin_edges = np.append(bin_edges, 80)
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
    pdf0 = pdf[0:][np.argsort(pdf[1:,0])] ## Make sure delta is increasing

    # First interpolate ln(mu) based on 1+delta, get rid of bins with 0
    overdensity_shape = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,1][pdf0[:,4]>10],k=1)
    shape = overdensity_shape(density_field)

    overdense_loc = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,2][pdf0[:,4]>10],k=1)
    loc = overdense_loc(density_field)

    overdense_scale = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,3][pdf0[:,4]>10],k=1)
    scale = overdense_scale(density_field)


    delta_min = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,5][pdf0[:,4]>10], k=1)
    min_points = delta_min(density_field)

    delta_max = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,6][pdf0[:,4]>10], k=1)
    max_points = delta_max(density_field)

    delta_mean = InterpolatedUnivariateSpline(pdf0[:,0][pdf0[:,4]>10]+1,pdf0[:,7][pdf0[:,4]>10], k=1)
    mean_points = delta_mean(density_field)

    fig = plt.figure(figsize=(12,14))
    ax1 = fig.add_subplot(611)
    ax2 = fig.add_subplot(612)
    ax3 = fig.add_subplot(613)
    ax4 = fig.add_subplot(614)
    ax5 = fig.add_subplot(615)
    ax6 = fig.add_subplot(616)


    ax1.plot(np.log10(pdf0[:,0]+1),pdf0[:,1],"k-", label = 'Data')
    ax1.scatter(np.log10(density_field[::1000]),shape[::1000],color = 'red', label = 'interpolation')

    ax2.plot(np.log10(pdf0[:,0]+1),pdf0[:,2],"k-")
    ax2.scatter(np.log10(density_field[::1000]),loc[::1000],color='red')


    ax3.plot(np.log10(pdf0[:,0]+1),pdf0[:,3],"k-")
    ax3.scatter(np.log10(density_field[::1000]),scale[::1000],color='red')


    ax4.plot(np.log10(density_field[::1000]),min_points[::1000],"r.")
    ax4.plot(np.log10(pdf[:,0]+1),pdf[:,5],"k-", lw= 2.5)

    ax5.plot(np.log10(density_field[::1000]),max_points[::1000],"r.")
    ax5.plot(np.log10(pdf[:,0]+1),pdf[:,6],"k-", lw= 2.5)

    ax6.plot(np.log10(density_field[::1000]),mean_points[::1000],"r.")
    ax6.plot(np.log10(pdf[:,0]+1),pdf[:,7],"k-", lw= 2.5)


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

    Nhalo = np.zeros_like(density_field)
    # Create a mask for valid scale values (scale > 0)
    mask = (scale > 0) & (shape>0)

    # Apply the lognormal sampling only where the mask is True
    if np.any(mask):
        Nhalo[mask] = lognorm.rvs(shape[mask], loc=loc[mask], scale=scale[mask])

    complex_condition = (scale > 0) & (shape>0) & (min_points<max_points) & (min_points<mean_points) & (mean_points<max_points)  & ((Nhalo>max_points) | (Nhalo<min_points))
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
                possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=1000000) 
                possible_values = possible_values[(possible_values>min_points[k]) & (possible_values<max_points[k])]
            Nhalo[k] = np.random.choice(possible_values)
            
    for k in simple_mask_index:
        while (Nhalo[k]<min_points[k]):
            possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=1000) 
            possible_values = possible_values[(possible_values>min_points[k])]
            if len(possible_values)==0:
                possible_values = lognorm.rvs(shape[k], loc=loc[k], scale=scale[k], size=100000) 
                possible_values = possible_values[(possible_values>min_points[k])]
            Nhalo[k] = np.random.choice(possible_values)

    step = 0.1
    binning=np.arange(-1,10.0+step,step)
    binning=np.append(binning,15)
    binning = np.append(binning, 100)

    prob_zero_fcoll = {"%0.2f" %bin:[] for bin in binning} #probability of having a zero coll fraction for a given delta


    for binindex in range(len(binning)-1):
        fcoll_current_bin = fcoll[((overdense-1)>binning[binindex]) & ((overdense-1)<binning[binindex+1])]
        nhalo_zero_fcoll = fcoll_current_bin[fcoll_current_bin==0]
        if len(fcoll_current_bin)==0:
            prob_zero_fcoll["%0.2f" %binning[binindex]].append(0)
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
    ax1.set_title('HESTIA 256**3')
    ax1.set_xlabel('log10(density_field)')
    ax1.set_ylabel('log10(Nhalo_new)')
    fig.colorbar(im1, ax=ax1, label='Density')

    # Plot the second heatmap
    im2 = ax2.imshow(heatmap2.T, origin='lower', cmap='seismic', aspect='auto',
                    extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]], norm='log')
    ax2.plot(np.log10(pdf0[:,0]+1),np.log10(pdf0[:,6]),color = 'yellow', alpha=1, lw=3)
    ax2.plot(np.log10(pdf0[:,0]+1),np.log10(pdf0[:,5]),color = 'pink', alpha=1, lw=3)

    ax1.plot(np.log10(pdf0[:,0]+1),np.log10(pdf0[:,6]),color = 'yellow', alpha=1, lw=3)
    ax1.plot(np.log10(pdf0[:,0]+1),np.log10(pdf0[:,5]),color = 'pink', alpha=1, lw=3)


    ax1.scatter(np.log10(density_field[::1000]),np.log10(max_points[::1000]))
    ax1.plot(np.log10(pdf[:,0]+1),np.log10(pdf[:,6]),"yellow", lw= 2.5)


    ax2.set_title('High Res simulation')
    ax2.set_xlabel(r'log$_{10}(/delta+1)$')
    ax2.set_ylabel(r'log$_{10}(f_{coll})$')
    fig.colorbar(im2, ax=ax2, label='Density')

    ax1.set_ylim(-3,1)
    ax2.set_ylim(-3,1)

    ax1.set_xlim(-0.5,2)
    ax1.set_xlim(-0.5,2)

    # Show the plots
    plt.savefig('./results/fcoll_implementation_z{:.3f}.png'.format(zred_HESTIA[j]))
    np.save('./results/fcoll_z{:.3f}'.format(zred_HESTIA[j]), Nhalo_new)
