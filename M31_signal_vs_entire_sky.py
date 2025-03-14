from scipy.interpolate import interp1d
import seaborn as sns
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--test-patches', action='store_true', help='testing by showing patches on sphere as they get smaller')
parser.add_argument('--show-dispersion', action='store_true', help='Also give dispersion plot of Rotation Measure within Halo of Andromeda')
parser.add_argument('--annuli-anal', action='store_true', help='Conducting annulus analysis for histograms')
parser.add_argument('--rm-vs-proj-dist', action='store_true', help='Honours output in plotting RM against projected distance from M31 (as well as assemsemnt of the entire RM-sky)')
args = parser.parse_args()

from main import (
#Importing alias'
np, u, SkyCoord, plt, stats, WCS, warnings, fits,

#importing variables
d_bin_centers as d_bin_centers_m31,
bin_means as bin_means_m31,
bin_std as bin_std_m31,
bin_med as bin_med_m31,
rm, rm_err, eq_pos,
m31_sep_Rvir, rm_m31,
bin_num as bin_num_from_main,
bin_std_past_rvir, L_m31, cutoff,

#Thes eare statisics for the entire sepration from m31 center to all 
bin_means_past_rvir, bin_meds_past_rvir,

#importing functions
get_projected_d_old, confining_circle
)

#Hanlding unnneccesary clutter of printing from warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter("ignore", RuntimeWarning) #Suppresses stats SmallSampleWarning when using stats.sem()


def curr_dir_path():
    """Returns the folder path of the currently running script."""
    return os.path.dirname(os.path.abspath(__file__)) + "/"

def interpolate(data, num_points):
    #Function to quadratically interpolate data (increasing number of points for smoothness)
    x = np.arange(len(data))
    
    data_clean = np.nan_to_num(data, nan=np.nanmean(data)) #Replacing NaNs with the mean of the data
    f_interpolate = interp1d(x, data_clean, kind='quadratic')
    x_new = np.linspace(x.min(), x.max(), num_points)
    return f_interpolate(x_new)
    
def get_real_rm_data(): #get rm_data from catalog.dat
    return rm, rm_err

#to be used for scaling size of scatter based on their rm value
sizes = lambda vals: 30 #np.abs(vals) * 2.3
        
#If overlap = False then minimum separation distance between points must be given
#this is used for the real data points too
#NOTE: The global vairable patch_size is used in here as min_distance_deg
def get_random_points(num=30000, ra_range=(0, 360), dec_range=(-1, 1), 
                      for_test=False, overlap=True, **kw):
    
    ra_range = [ra_range[0], ra_range[1]] #In degrees
    if dec_range != (-1, 1):  #Handling custom Dec range
        dec_range = np.deg2rad(dec_range) #Conversion to rads
        dec_range = np.sin(dec_range) #conversion to sin(dec)
    if overlap:  # If overlapping is allowed
        if for_test:
            L = []
            # Generating random RA and Dec points
            ra_vals = np.random.uniform(*ra_range, num) * u.deg  # RA points in degrees
            
            # Ensuring that arcsin is within [-1, 1]
            dec_vals_clipped = np.clip(np.random.uniform(*dec_range, num), -1, 1)
            dec_vals = np.rad2deg(np.arcsin(dec_vals_clipped)) * u.deg  # Dec points
            L.append(ra_vals)
            L.append(dec_vals)

            # Dummy RM values and errors for testing
            L.append(np.random.choice([-1, 1], size=num))  # Random RM signs (-1 or 1)
            L.append(np.random.uniform(-0.01, 0.01, num))  # Random small errors

        else:
            valid_ra, valid_dec = [], []
            radius = kw["radius"] if "radius" in kw else 0
            count = 0
            min_distance_deg = kw["min_distance_deg"] if "min_distance_deg" in kw else 0
            while len(valid_ra) < num:
                
                remaining = num - len(valid_ra)  # How many more points to generate
                ra = np.random.uniform(*ra_range, remaining)  # Random RA points
                
                #Esuring dec values are valid for arcsin and the ENTIRE patch is above -40 dec
                sin_dec_min = np.sin(np.radians(-40)) + np.radians(min_distance_deg)
                dec_vals_clipped = np.clip(np.random.uniform(*dec_range, remaining), sin_dec_min, 1)
                dec = np.arcsin(dec_vals_clipped)  # Random Dec points in radians

                # #Making sure the circle patches dont include the empty Dec region
                # #Their centres are alsready never in this region but their areas might be, so this mucst be handeled
                # if n:=sum(dec - np.radians(min_distance_deg) < np.sin(-40)) > 0: continue#Hence num value remains the same in next loop
                
                valid_mask = (np.rad2deg(dec) + radius <= 90) & (np.rad2deg(dec) - radius >= -90)
                # if not (n:= np.sum(valid_mask)): print(n) ; import sys; sys.exit()

                dec = np.rad2deg(dec)  # Random Dec points
                valid_ra.extend(ra[valid_mask])
                valid_dec.extend(dec[valid_mask])
                count += 1
            L = np.array(valid_ra[:num]), np.array(valid_dec[:num])
        
        return L
    
    else: #If overlapping of pathces is prohibited
        min_distance_deg=kw["min_distance_deg"]
        min_distance_rad = np.radians(min_distance_deg) #dealing with radians
        
        #Lists to store accepted RA and Dec values
        ra_values = []
        dec_values = []
        
        while len(ra_values) < num: #Generating points one by one to make sure there is no overlapping of patches
            
            #Generating random point
            ra_new = np.random.uniform(*ra_range) 
            
            # Clip the Dec point to stay within the specified range
            #"loc" is he mean and "scale" contols the spread about that mean (stdev)
            dec_new = np.clip(np.random.normal(loc=0, scale=30), *dec_range)
            
            if (ra_new + min_distance_deg > 360 or ra_new - min_distance_deg < 0) or (dec_new + min_distance_deg > 90 or dec_new - min_distance_deg < -90):
                continue

            if ra_values: #If not empty
                #Convert current points and new point to radians
                ra_values_rad = np.radians(ra_values)
                dec_values_rad = np.radians(dec_values)
                ra_new_rad = np.radians(ra_new)
                dec_new_rad = np.radians(dec_new)
                
                #Compute angular distance using spherical law of cosines
                cos_distances = (np.cos(dec_values_rad) * np.cos(dec_new_rad) +
                             np.sin(dec_values_rad) * np.sin(dec_new_rad) *
                             np.cos(ra_values_rad - ra_new_rad))
                
                #Convert cosines to angular distances in radians
                distances = np.arccos(cos_distances)
                
                #If any distance is smaller than the minimum, reject this point
                #Also making sure it doesnt surpass declination upper limit since this is a point for the corner of the patch
                #And this square patch will have a height and width of min_distance_deg
                if np.any(distances < min_distance_rad):
                    continue #try plottting another point
            
            
            #If all distances are valid, add the point to the list
            ra_values.append(ra_new)
            dec_values.append(dec_new)
        
        return ra_values, dec_values

def get_wcs(filename):
    hdu = fits.open(filename)[0]
    return WCS(hdu.header)
    
def collection_of_points_from_WCS_sphere():#Placing patches on sphere with real RM values

    for patch_ra, patch_dec in list(zip(patch_ra_points, patch_dec_points)):
        
        (patch_region_coord_pos,
          rm_val_filament_pos,
          patch_region_coord_neg,
          rm_val_filament_neg
          ) = confining_circle(ax1, #needed to use WCS tools from astropy 
                                ra=patch_ra, #degrees
                                dec=patch_dec, #degrees
                                radius=patch_size, #radius of the circle
                                angle=7, #used for addition onto polar angle (may be unnecessary to be honest)
                                polar_angle=-10, #degrees 
                                positions=eq_pos, #all points on WCS axis
                                pos_mask=rm_s>0, #All positive points on WCS axis
                                neg_mask = rm_s<0, #All negative points on WCS axis
                                rm_s=rm_s, 
                                rm_errs=rm_errs,
                                return_data=True, return_err=False,
                                plot=(False, #plotting patch 
                                      False) #plotting scatter within patch
                                )
        
        RM_values = list(rm_val_filament_pos)
        RM_values.extend(rm_val_filament_neg)
        RM_values_per_patch.append(RM_values)
        
        RM_coords = SkyCoord([patch_region_coord_pos, patch_region_coord_neg])
        
        #Essential grouping of skycoords to Skycoord list
        RM_coords =  np.vstack((RM_coords.ra.deg, RM_coords.dec.deg)).T
        
        try:
            RM_coords = SkyCoord(ra=RM_coords[:, 0], dec=RM_coords[:, 1], unit=(u.deg, u.deg), frame='icrs')
        except IndexError: #If RM_coord is empty fill with nan value.
            RM_coords = SkyCoord(ra=np.nan, dec=np.nan, unit=(u.deg, u.deg), frame='icrs')
    
        RM_coords_per_patch.append(RM_coords)
        
def get_mean_and_med_stats(sep_vals, rm_vals, bin_num):
    rm_vals = np.abs(rm_vals) # |RM|
    
    #Calculating mean and median of RM values within patch
    bin_means, bin_edges, binnumber = stats.binned_statistic(sep_vals, rm_vals, statistic='mean', bins=bin_num)
    bin_med, bin_edges, binnumber = stats.binned_statistic(sep_vals, rm_vals, statistic='median', bins=bin_num)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2
    
    #Converting to projected distance away from center of m31(in kpc)
    d_bin_centers = get_projected_d_old(bin_centers).value
    
    #Standard error of the mean for error bars
    bin_std, bin_edges, binnumber = stats.binned_statistic(
        sep_vals, rm_vals,
        statistic=stats.sem,  # Standard Error of Mean
        bins=bin_num)

    return d_bin_centers, bin_means, bin_med, bin_std


def annuli_analysis(save_plot=False, stack_indiv_patch=False): 
    print("Annuli Histogram analysis and plotting have begun")
    
    def construct_and_plot_annuli(distance_flat, rm_vals_flat_mean, rm_vals_flat_median, save_plot=save_plot):
        distance_flat = np.asarray(distance_flat)
        rm_vals_flat_mean = np.asarray(rm_vals_flat_mean)
        rm_vals_flat_median = np.asarray(rm_vals_flat_median)

        #<bin_num_from_main> number of bins for annular regions
        bin_edges = np.linspace(distance_flat.min(), distance_flat.max(), bin_num_from_main+1) #bin number based 
        bin_indices = np.digitize(distance_flat, bins=bin_edges) #Assigning each RM to a bin
        
        annul_area = np.pi * (bin_edges[1]**2 - bin_edges[0]**2) #Areas are essentially discrete and are the same.
        # Storing RM values per annulus
        rm_per_annulus_mean = {i: [] for i in range(1, len(bin_edges)+1)}  # Dictionary to store RM values for each bin (mean)
        rm_per_annulus_median = {i: [] for i in range(1, len(bin_edges)+1)}  # Dictionary to store RM values for each bin (median)

        for i, bin_idx in enumerate(bin_indices):
            rm_per_annulus_mean[bin_idx].append(rm_vals_flat_mean[i]/annul_area)  # Append RM values for mean to corresponding bin
            rm_per_annulus_median[bin_idx].append(rm_vals_flat_median[i]/annul_area)  # Append RM values for median to corresponding bin

        # Converting lists to numpy arrays for easier handling
        rm_per_annulus_mean = {k: np.array(v) for k, v in rm_per_annulus_mean.items() if len(v) > 0}
        rm_per_annulus_median = {k: np.array(v) for k, v in rm_per_annulus_median.items() if len(v) > 0}

        if stack_indiv_patch: 
            annuli_to_plot = np.arange(0, round(cutoff.value)+1) # Adjusting annuli ranges based on bin edges (in degrees)
        else: 
            annuli_to_plot = np.arange(0, bin_num_from_main+1)  # Adjusting annuli ranges based on bin edges (in kpc)

        anul_dist_type = "deg" if stack_indiv_patch else "kpc"
        b_m_1 = bin_means_m31 if stack_indiv_patch else bin_means_past_rvir
        b_m_2 = bin_med_m31 if stack_indiv_patch else bin_meds_past_rvir
        std = bin_std_m31 if stack_indiv_patch else bin_std_past_rvir #assuming same standard deviation

        # Loop through the annuli
        for bin_idx in annuli_to_plot:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Create 1x2 subplot grid
            
            x_axis_label = r" (rad/m$^{2}$/$\xi$)"
            histbin = 50 #number of beans per histogram's RM-axis (x-axis)
            if bin_idx in rm_per_annulus_mean and bin_idx in rm_per_annulus_median:

                # Plot for "Mean" subplot (left side)
                counts, _, _ = axes[0].hist(rm_per_annulus_mean[bin_idx], bins=histbin, alpha=0.5)
                axes[0].set_title("Mean")
                axes[0].set_xlabel("RM per radial area " +  x_axis_label)
                axes[0].set_ylabel("Counts")
                axes[0].set_ylim(0, np.max(counts) * 1.1)
                
                # Plot for "Median" subplot (right side)
                counts, _, _ = axes[1].hist(rm_per_annulus_median[bin_idx], bins=histbin, alpha=0.5)
                axes[1].set_title("Median")
                axes[1].set_xlabel("RM per radial area " +  x_axis_label)
                # axes[1].set_ylabel("Counts")
                axes[1].set_ylim(0, np.max(counts) * 1.1)

                # For M31 relative annulus (mean RM)
                axes[0].axvline(x=b_m_1[bin_idx-1]/annul_area, 
                                label=f"m31 = {b_m_1[bin_idx-1]/annul_area:.4f}", 
                                color='k', linestyle='--', linewidth=.8)
                axes[1].axvline(x=b_m_2[bin_idx-1]/annul_area, 
                                label=f"m31 = {b_m_2[bin_idx-1]/annul_area:.4f}", 
                                color='k', linestyle='--', linewidth=.8)

                # Filling the region around the mean RM value within 1 sigma (mean RM)
                axes[0].fill_betweenx(
                    y=np.linspace(0, 700, 100),
                    x1=(b_m_1[bin_idx-1] - std[bin_idx-1])/annul_area,
                    x2=(b_m_1[bin_idx-1] + std[bin_idx-1])/annul_area,
                    color='k', alpha=0.2, edgecolor="none",
                    label=r"$1\sigma \approx {:.2f}$".format(std[bin_idx-1]/annul_area)
                )
                # Filling the region around the median RM value within 1 sigma (median RM)
                axes[1].fill_betweenx(
                    y=np.linspace(0, 700, 100),
                    x1=(b_m_2[bin_idx-1] - std[bin_idx-1])/annul_area,
                    x2=(b_m_2[bin_idx-1] + std[bin_idx-1])/annul_area,
                    color='k', alpha=0.2, edgecolor="none",
                    label=r"$1\sigma \approx {:.2f}$".format(std[bin_idx-1]/annul_area)
                )

                axes[0].legend(); axes[1].legend()
                fig.suptitle(f"Range: {bin_edges[bin_idx-1]:.2f} - {bin_edges[bin_idx]:.2f} {anul_dist_type}" + r" ($\xi$" + f" ={annul_area:.2f}" + f"{anul_dist_type}" r"$^2$)")
                
                # Save or display the figure
                if save_plot:
                    path = curr_dir_path() + "Results/"
                    plt.savefig(f"{path}annuli_plot_{bin_idx}.png", dpi=600, bbox_inches="tight")
                    plt.clf()  # clearing the figure (not deleting it)
                else:
                    plt.show()

                # break #Testing out one plot

        if save_plot: 
            plt.close()  # Deleting the figure to clear memory
            print(f"All images saved to {path}")
    
    # Stack all patches together without any mean analysis
    if stack_indiv_patch:
        
        # Converting Radial separation from relative patch to projected distance to be used for BG correction
        projected_distances = [
            [get_projected_d_old(val) for val in sublist.value] 
            for sublist in RM_coords_sep
        ]
        
        RM_values_per_patch_corr = [
            [indiv_bg_corr(RM_val, bin_cent=proj_d_val, absol=False) for RM_val, proj_d_val in zip(RM_patch, proj_d_patch)]
            for RM_patch, proj_d_patch in zip(RM_values_per_patch, projected_distances)
        ]

        # Flattening the lists for easier computation
        flat_sep_vals = np.concatenate([patch.value for patch in RM_coords_sep])  # Separation distances (degrees)
        flat_rm_vals_mean = np.concatenate([RM_values_per_patch_corr_mean for RM_values_per_patch_corr_mean in RM_values_per_patch])  # Corresponding RM values for mean
        flat_rm_vals_median = np.concatenate([RM_values_per_patch_corr_median for RM_values_per_patch_corr_median in RM_values_per_patch])  # Corresponding RM values for median
        construct_and_plot_annuli(flat_sep_vals, flat_rm_vals_mean, flat_rm_vals_median)

    # Using binned values for MEAN and MEDIAN (As initially discussed with DJ and Russ)
    else: 
        D_Bin_centers = np.concatenate([bin_centers for bin_centers in all_d_bin_centers])
        Avg_Means = np.concatenate([avg_mean for avg_mean in all_means_1]) 
        Avg_Medians = np.concatenate([avg_med for avg_med in all_medians_1]) 
        construct_and_plot_annuli(D_Bin_centers, Avg_Means, Avg_Medians)

        
def plot_indidividual_patch_stats(ax, d_bin_centers, bin_mean, bin_med, bin_std):
    # Plot the mean and median with error bars
    ax.errorbar(d_bin_centers, np.absolute(bin_mean), yerr=bin_std, fmt='k.-', alpha=.4)
    ax.errorbar(d_bin_centers, np.absolute(bin_med), yerr=bin_std, fmt='g.-', alpha=.4, capsize=2)

def plot_m31_stats(ax):
    ax.errorbar(d_bin_centers_m31, np.absolute(bin_means_m31), yerr = bin_std_m31, 
                color="black", fmt = '.-', alpha=.6)#, label ="$\mu_{M31}$")
    ax.errorbar(d_bin_centers_m31, np.absolute(bin_med_m31), yerr=bin_std_m31, 
                color='orange', fmt='.-', capsize=2, markeredgecolor="k", alpha=.6)

def plot_m31_dispersion(bin_num, save_plot=False):

    plt.figure(figsize = (10, 6))

    #Dispersion of RM values in each bin (Standard Error of the Means)
    bin_std, _, _ = stats.binned_statistic(m31_sep_Rvir, rm_m31, statistic=stats.sem, bins = bin_num)
    plt.plot(d_bin_centers, bin_std, "ko")
    plt.xlabel('Projected Separation from M31[kpc]')
    plt.ylabel('$\sigma_{RM} [\mathrm{rad} \ \mathrm{m}^{-2}]$', rotation='horizontal', labelpad=60, fontsize=15)

    x_values = np.linspace(0, 296, 1000)

    plt.grid(True)
    plt.tight_layout()

    # #Plotting the curve_fit
    # coefficients = np.polyfit(d_bin_centers, bin_std, 3)
    # fit_line = np.poly1d(coefficients)
    # plt.plot(d_bin_centers.value, fit_line(d_bin_centers.value), 
    #             color = 'orange', linestyle = "--")

    plt.title("Dispersion of RM values in each bin")

    #Mentioned on page 841 of https://doi.org/10.1093/mnras/stad2811
    plt.axhline(xmin=0, xmax = np.max(x_values), y=6, linestyle="--",
                label="Observed $\sigma_{RM}$ indepenedent of Galactic Latitude")

    plt.legend(fontsize = 12, loc = 'upper center', bbox_to_anchor = (0.5, 1.2),
                framealpha = 0, ncols = (2,2))
    
    if save_plot: 
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}" + "Dispersion_M31.png", dpi=600, bbox_inches="tight")#Saving as image
        plt.clf() #clearing the figure (not deleting it)
    else:
        plt.close() #Deleting the figure to clear memory
        print(f"Dispersion image saved in " + f"{path}")

# In[A little extra testing]:        
#This code below produces 2 subplots on a figure.
#First plot is sphere with unformly random points with patches of equal area placed on the surface of the sphere
#Second is a 3d scatter plot representing how many points were collected by each patchcollected (with a given ra and dec)
def test_patches_on_sphere():
    random_ra, random_dec, random_mag, random_err = get_random_points(num=30000, for_test=True)
    
    positive_mask = random_mag > 0
    negative_mask = random_mag < 0
    
    #Converting the list of tuples to SkyCoord objects
    random_SC_coords = SkyCoord(ra=random_ra, dec=random_dec, unit=(u.deg, u.deg), frame='icrs')
    
    counter = 0
    # Testing whether patches of equal area collect equal points on spherical surface
    for step in np.arange(50, 0, -1):
        #Creating figure and axis with WCS projection
        fig = plt.figure(figsize=(16,5))
        ax = fig.add_subplot(131, projection=get_wcs("LGSNLGSR.SQLGBB.FITS"), 
                            slices=('x', 'y', 0, 0))
        ax.grid(True, alpha=1, color = "k")
        ax.set_aspect('equal') # Making Sphere purely spherical-looking
        # ax.grid(True, alpha=0.5, color='k')
        
        try:
            ax.scatter(random_ra, random_dec, s=1, marker=".", 
            transform=ax.get_transform('world'))
        except ValueError:
            print(random_ra.shape, random_dec.shape) 
        
        ax.set_xlabel('Right Ascension (RA)')
        ax.set_ylabel('Declination (Dec)')
        ax.set_title('Random Sky Coverage on WCS Projection')
        
        # #ZOOMING OUT (in degrees)
        # ra_min, ra_max = 80, -20
        # dec_min, dec_max = 16, 55
        
        # # Convert these to pixel coordinates using WCS
        # bottom_left = wcs.world_to_pixel_values(ra_min, dec_min, 0, 0)
        # top_right = wcs.world_to_pixel_values(ra_max, dec_max, 0, 0)
        
        # #Extracting pixel coordinates
        # x_min, y_min, _, _ = bottom_left
        # x_max, y_max, _, _ = top_right
        
        # Setting axis limits
        # ax.set_xlim(x_min, x_max)
        # # ax.set_ylim(y_min, y_max)
        
        dec = 90
        storage=[]
        for ra in range(0, 360+step, step):
            # if ra + step > 360 or ra - step < 0:
            #     storage.append((np.nan, ra, dec))
            #     continue
            
            step_ra = step #remember steps taken by ra alone
    
            for dec in range(-90, 90+step, step):
    
                # #Make sure it doesnt surpass declination upper limit
                # if dec + step > 90 or dec - step < -90 or np.abs(dec) + step > 90:
                #     storage.append((np.nan, ra, dec))
                #     continue
            
                (patch_region_coord_pos,
                rm_val_filament_pos,
                patch_region_coord_neg,
                rm_val_filament_neg
                ) = confining_circle(ax, 
                                    ra=ra, #degrees
                                    dec=dec, #degrees
                                    radius=step, #radius of the circle
                                    angle=7, #used for addition onto polar angle (may be unnecessary to be honest)
                                    polar_angle=-10, #degrees 
                                    positions=random_SC_coords, #all points on WCS axis
                                    pos_mask=positive_mask, #All positive points on WCS axis
                                    neg_mask = negative_mask, #All negative points on WCS axis
                                    rm_s=random_mag, 
                                    rm_errs=random_err,
                                    return_data=True, return_err=False,
                                    plot=(True, #plotting patch 
                                            True) #plotting scatter within patch
                                    )
            
                #Store the dec value after the counts (counts, dec) tuple
                #Notice that this is for the same patch that has both ra and dec widths
                storage.append((len(rm_val_filament_neg)+len(rm_val_filament_pos), ra, dec))
            break
            step = step_ra #get back on steps taken in ra
    
        storage = np.array(storage)  # Ensure it's a numpy array
        valid_data = storage[~np.isnan(storage).any(axis=1)]  # Filter out NaN values
    
        COUNTS = valid_data[:, 0]  # First column: counts
        RA = valid_data[:, 1]       # Second column: RA
        DEC = valid_data[:, 2]      # Third column: Dec
        
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133)
        
        n_bins = 100
        
        # Bin edges
        ra_bins = np.linspace(0, 360, n_bins + 1)  # RA bins from 0 to 360 degrees
        dec_bins = np.linspace(-90, 90, n_bins + 1)  # DEC bins from -90 to 90 degrees
        
        # Sum counts per RA bin
        ra_hist, _ = np.histogram(RA, bins=ra_bins, weights=COUNTS)
        
        # Sum counts per DEC bin
        dec_hist, _ = np.histogram(DEC, bins=dec_bins, weights=COUNTS)
        
        # Histogram for RA bins
        ax1.bar(ra_bins[:-1], ra_hist, width=0.4*step, edgecolor='black', align='edge', color='blue')
        ax1.set_title('Sum of Counts per RA Bin')
        ax1.set_xlabel('Right Ascension (degrees)')
        ax1.set_ylabel('Total Counts (RA)')
        ax1.set_aspect('auto')
        
        # Histogram for DEC bins
        ax2.bar(dec_bins[:-1], dec_hist, width=0.4*step, edgecolor='black', align='edge', color='green')
        ax2.set_title('Sum of Counts per DEC Bin')
        ax2.set_xlabel('Declination (degrees)')
        ax2.set_ylabel('Total Counts (DEC)')
        ax2.set_aspect('auto')
        
        # ax.view_init(elev=10, azim=20)  # Adjust azim to rotate around the z-axis
        plt.suptitle(f"Number {counter}")
        plt.tight_layout()
        plt.show()
        counter += 1

def Shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_Avg_medians, int_Avg_medians_std, int_D_bin_centers):
    def shade_ms_mimic(data,#average mean/ medians (interpolated)
                        std,#standard deviation (interpolated)
                        bin_centers, #Bin centers
                        cmap, name):

        y_upper= np.max(data + std)  #100

        X, Y = np.meshgrid(bin_centers, np.linspace(np.min(data-std), np.max(data+std), 100))

        #Computing dispersion-density based on the fill_between regions
        Z = np.exp(-((Y - data[:, None])**2 / (2 * std[:, None]**2)))

        fig, ax = plt.subplots(figsize=(10, 5))
        sns.heatmap(Z, cmap=cmap, alpha=0.6, cbar=True, xticklabels=False, yticklabels=False)
        y = np.arange(0,y_upper+10, 10)
        ax.set_yticks(y); ax.set_yticklabels(y)
        ax.set_xlabel(r'R$_{projected}$ [kpc]', fontsize=12)
        ax.set_ylabel('|RM|', fontsize=12)
        ax.set_ylim(0, y_upper)
        plt.title('Dispersion of RM Variations ('+name+')')

    shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_D_bin_centers, name="Means", cmap="viridis")
    shade_ms_mimic(int_Avg_medians, int_Avg_medians_std, int_D_bin_centers, cmap="viridis", name="Medians")

    plt.show()

#Individiual backgorund correction
def indiv_bg_corr(arr, bin_cent, absol=True):
    """
    bin_cent - Center of bin of projected distance or just projected distance relative to the RM values given in 'arr'
    arr - RM arrays . RM separation distances in degrees
    """

    # If arr and bin_cent are lists of arrays, process each separately
    if isinstance(arr, list) and isinstance(bin_cent, list):
        return [indiv_bg_corr(a, b, absol) for a, b in zip(arr, bin_cent)]
    
    # Ensure bin_cent is an array and force it to have units
    bin_cent = np.asarray(bin_cent) * u.kpc if not hasattr(bin_cent, "unit") else bin_cent
    bin_cent = bin_cent.to_value()  # Convert to unitless
    
    arr = np.asarray(arr)  # Ensure arr is also an array (important for indexing later)
    
    arr_bg = np.where(bin_cent > 300, arr, 0)  # Fill all values within virial radius with 0

    # Ensure arr_bg is at least 1D
    arr_bg = np.atleast_1d(arr_bg)

    arr_bg_no_0 = arr_bg[np.flatnonzero(arr_bg)]  # Remove zeros to compute mean of background
    BG = arr_bg_no_0.mean() if arr_bg_no_0.size > 0 else 0

    return np.abs(arr - BG) if absol else arr - BG

#Note that it might take too long to fit patches that dont overlap each other 
#If the number of patches are too many and/or the size of the patches are too big
patch_size = 30 #in degrees (same as M31 Virial Radius)

"""IMPORTANT"""
number_of_patches = int(1e4) #Creating laaaarge nubmer of patches (choose smaller vlue if you only want to see output features)


BINS = 50
rm_s, rm_errs = get_real_rm_data()

patch_ra_points, patch_dec_points = get_random_points(num=number_of_patches, 
                                                  ra_range=(min(eq_pos.ra.deg), 
                                                            max(eq_pos.ra.deg)),
                                                  dec_range=(min(eq_pos.dec.deg),
                                                              max(eq_pos.dec.deg)),
                                                  overlap=True, #Allowing circles to overlap
                                                  radius = patch_size, #in degrees
                                                  min_distance_deg = patch_size #minimum separation of points in degrees
                                                  )

Patch_pos = SkyCoord(ra=patch_ra_points, dec=patch_dec_points, unit=(u.deg,u.deg), frame='icrs')
print("Saved random points as SkyCoord")

# RM_coords_per_patch will be filled with values when "collection_of_points_from_WCS_sphere()" is called
# Same goes for RM_values_per_patch
RM_values_per_patch, RM_coords_per_patch = [], []

#Creating axis that will be used when doing WCS axis point collections
fig1 = plt.figure(figsize=(6, 5))
ax1 = fig1.add_subplot(111, projection=get_wcs("LGSNLGSR.SQLGBB.FITS"), slices=('x', 'y', 0, 0))
fig1.clf() #No need to show figure. Only axis is needed

print("(No longer plotting) but collecting of points from each circular patch has begun...")
collection_of_points_from_WCS_sphere() #IMPORTANT
print("Collection of points by each circular patch on WCS sphere is complete")

print("Getting separation of RM from center of relative patch")
#Get separation of RM from center of relative patch.
RM_coords_sep = [rm_coords.separation(patch_pos) 
                 for rm_coords, patch_pos in 
                 list(zip(RM_coords_per_patch, Patch_pos))]

all_d_bin_centers=[] #For x-axis
all_means = []
all_medians = []
all_bin_stds = []

fig2 = plt.figure(figsize=(10, 5))
ax2 = fig2.add_subplot(111)
print("Mean and Median calculations have begun")
for i in range(len(RM_coords_sep)): #Searching through each patch
    
    if not (RM_coords_per_patch[i].shape == ()):  #Checking if its not empty or filled with NaNs
        if not (RM_coords_per_patch[i].shape[0] < 15):  # Ensure sufficient number of points for stats module to work
            d_bin_centers, bin_mean, bin_med, bin_std = get_mean_and_med_stats(RM_coords_sep[i], RM_values_per_patch[i], bin_num=BINS)
            
            all_d_bin_centers.append(d_bin_centers) #For x axis
            all_means.append(bin_mean) #For y-axis
            all_medians.append(bin_med) #For y-axis
            all_bin_stds.append(bin_std)
            
            # #Background correction for the individual plots
            # bin_mean_1 = indiv_bg_corr(bin_mean, d_bin_centers)
            # bin_med_1 = indiv_bg_corr(bin_med, d_bin_centers)

            # #This has been commented out to remove clatter
            # #The they are all being collected and will be averaged to make a final one
            # plot_indidividual_patch_stats(ax2, d_bin_centers, bin_mean_1, bin_med_1, bin_std)
print("Mean and Median calculations have ended")

if __name__ == "__main__": #continue (this makes it easier to excecute "M31_signal_density.py" file)

    #getting mean of background
    D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                    max([max(centers) for centers in all_d_bin_centers]), 
                                    num=BINS)
    
    all_means_bg = np.where(D_bin_centers > 300, all_means, 0) #Fill all mean values within virial radius with 0
    all_medians_bg = np.where(D_bin_centers > 300, all_medians, 0) #Fill all median values within virial radius with 0 
    all_means_bg = [a[a != 0] for a in all_means_bg] #removing all the zeros to only focus on mean of backgrounds 
    all_medians_bg = [a[a != 0] for a in all_medians_bg]  #removing all the zeros to only focus on median of backgrounds
    BG_mean = np.mean(all_means_bg)
    BG_median = np.mean(all_medians_bg)

    Avg_means = np.mean(all_means, axis=0) - BG_mean
    Avg_medians = np.mean(all_medians, axis=0) - BG_median
    Avg_means_std = np.std(all_means, axis=0)
    Avg_medians_std= np.std(all_medians, axis=0)

    #Interpolating all data to have 100 poitns for more smoothness
    #Default type of smoothing is via quadratic interpolation
    Tup = ()
    for data in [Avg_means, Avg_medians, Avg_means_std, Avg_medians_std]:
        Tup += (interpolate(data, num_points=100),)
    int_Avg_means, int_Avg_medians, int_Avg_means_std, int_Avg_medians_std = Tup

    # Find a common range of D_bin_centers for non-interpolated data (BINS =30)
    D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                    max([max(centers) for centers in all_d_bin_centers]), 
                                    num=BINS #BINS #Number of points for interpolation (smoothens it out)
                                            #Must be same as bin_num parameter in function "get_mean_and_med_stats"
                                    )
    
    # Find a common range of D_bin_centers for interpolation
    int_D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                    max([max(centers) for centers in all_d_bin_centers]), 
                                    num=100 #BINS #Number of points for interpolation (smoothens it out)
                                            #Msut be same as bin_num parameter in function "get_mean_and_med_stats"
                                    )
    
    plot_m31_stats(ax2) #Plots the data from intial starterkit (So nothing new here)

    ax2.errorbar(D_bin_centers, np.absolute(Avg_means), yerr = Avg_means_std, fmt = 'b.-')#,label ="$\mu_{\mu patch}$"
    ax2.errorbar(D_bin_centers, np.absolute(Avg_medians), yerr = Avg_medians_std, 
                color= 'green', fmt='.-', capsize = 2)#,label ='$Median_{\mu patch}$')
    
    # #plotting average (Interpolated for smoother effect)
    # ax2.plot(D_bin_centers, Avg_means, color='blue', label='$\mu_{\mu patch}$')
    # ax2.plot(D_bin_centers, Avg_medians, color='green', label='$Median_{\mu patch}$')

    # After all plots, fill the area between the highest and lowest lines
    ax2.fill_between(int_D_bin_centers, 
                    int_Avg_means - int_Avg_means_std , 
                    int_Avg_means + int_Avg_means_std, 
                    color='blue', alpha=0.2)#, label='$\mu_{patch}$' + ' Coverage')
    ax2.fill_between(int_D_bin_centers, 
                    int_Avg_medians - int_Avg_medians_std, 
                    int_Avg_medians + int_Avg_medians_std, 
                    color='green', alpha=0.4)#, label='$Median_{patch}$' + ' Coverage')
    ax2.set_xlabel(r'R$_{projected}$ [kpc]',fontsize=12)
    ax2.set_ylabel('|RM|', fontsize=12)
    ax2.set_xlim(0,300)
    ax2.set_ylim(0,)

    # ax2.fill_between(D_bin_centers, Avg_means_std , Avg_means_std, color='blue', alpha=0.2, label='Mean standard deviation of the patches')
    # ax2.fill_between(D_bin_centers, Avg_medians_std, Avg_medians_std, color='green', alpha=0.4, label='Median standard deviation of the patches')

    # ax2.legend(fontsize = 12, loc = 'upper center', bbox_to_anchor = (.5, 1.2), 
    #             framealpha = 0, ncols = (2,4))
    plt.tight_layout()
    #plt.show()
    
    if args.rm_vs_proj_dist:
        path = curr_dir_path() + "Results/"
        fig2.savefig(f"{path}M31_signal_vs_entire_sky_{number_of_patches}_patches.png", dpi=600, bbox_inches="tight") #saving the image
        print(f"M31_signal_vs_entire_sky_{number_of_patches}_patches.png has been successfully saved in Results directory")
    
    plt.close(fig2) #Deleting the Figure
    
    # #Data close to Gundo's shade_ms plots to spot any storng outliers
    # Shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_Avg_medians, int_Avg_medians_std, int_D_bin_centers)

    if args.test_patches: #show patches on sphere as they get smaller
        test_patches_on_sphere()

    if args.show_dispersion:
        plot_m31_dispersion(bin_num_from_main)

    all_means_1 = indiv_bg_corr(all_means, all_d_bin_centers, absol=False)
    all_medians_1 = indiv_bg_corr(all_medians, all_d_bin_centers, absol=False)
    #MASTERS addition to identifying significance in M31's halo compared to sky via annulus analysis
    if args.annuli_anal: annuli_analysis(save_plot=True)