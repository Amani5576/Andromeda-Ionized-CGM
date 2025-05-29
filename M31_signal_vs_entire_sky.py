from scipy.interpolate import interp1d
import seaborn as sns
import os
import matplotlib.animation as animation
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from numpy.polynomial.polynomial import Polynomial
import pickle
from concurrent.futures import ThreadPoolExecutor
from matplotlib.ticker import MaxNLocator
import pandas as pd
from tabulate import tabulate

from main import (
#Importing alias'
np, u, SkyCoord, plt, stats, WCS, warnings, fits, Patch,

#importing variables
d_bin_centers as d_bin_centers_m31,
bin_means as bin_means_m31,
bin_med as bin_med_m31,
bin_std as bin_std_m31,
bin_edges_mean_m31,
bin_edges_median_m31,
rm, #RM values of the entrie sky
rm_err, 
eq_pos,#RM coordinates of the entrie sky in ICRS coordinates
m31_sep_Rvir, rm_m31, err_m31,
m31_sep_bg, rm_bg, err_bg,
bin_num as bin_num_from_main,
bin_std_past_rvir, L_m31, cutoff,
bin_num as bin_num_main,
PA_rm, PA_bg,
R_vir, #Virial radius of M31 in astropy.units of kpc
d_m31, #Distance to Andormeda in kpc
L_m31, #Rvir in degrees
cutoff, #Cutoff value for background region limit in degrees
rm_pos_gal_lat, #Galactic Latitudes of RM's within virial radius
rm_pos_gal_lat_bg, #Galactic Latitudes of RM's for M31's Background
#These are statisics for the entire sepration from m31 center to all 
bin_means_past_rvir, bin_meds_past_rvir,
m31_lat, #Galactic Latitude of M31 in degrees
b_upper_lim, #Limit of |b|>5 degrees
args,

#importing functions
get_projected_d, get_sep_angle,
confining_circle, get_wcs, 
BG_correction, curr_dir_path, 
create_annuli_binning_structure,
apply_plot_attributes,
get_discrete_colors,
plot_m31_stats,
curr_dir_path,
get_discrete_colors,
)

#Hanlding unnneccesary clutter of printing from warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter("ignore", RuntimeWarning) #Suppresses stats SmallSampleWarning when using stats.sem()

def load_pickle(filename, message):
    with open(filename, "rb") as f:
        data = pickle.load(f)
    print(message)
    return tuple(data.values()) 

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
def get_random_points(num=30_000, ra_range=(0, 360), dec_range=(-1, 1), 
                      for_test=False, overlap=True, **kw):
    
    num_fake = num
    ra_range = list(ra_range) #Should be in degrees (without the units)
    if dec_range != (-1, 1):  #Handling custom Dec range
        dec_range = np.deg2rad(dec_range) #Conversion to rads
        dec_range = np.sin(dec_range) #conversion to sin(dec)
    if overlap:  # If overlapping is allowed
        if for_test:
            L = []
            # Generating random RA and Dec points
            ra = np.random.uniform(*ra_range, num) * u.deg  # RA points in degrees
            
            # Ensuring that arcsin is within [-1, 1]
            dec_vals = np.random.uniform(*dec_range, num)
            dec_vals_clipped = np.clip(dec_vals, -1, 1)
            dec = np.rad2deg(np.arcsin(dec_vals_clipped)) * u.deg  # Dec points
            L.append(ra)
            L.append(dec)

            # Dummy RM values and errors for testing
            L.append(np.random.choice([-1, 1], size=num))  # Random RM signs (-1 or 1)
            L.append(np.random.uniform(-0.01, 0.01, num))  # Random small errors

        else:
            valid_dec = []
            radius = kw.get("radius", 10) #Radius of Random patch to be in DEGREES
            count = 0
            min_distance_rad = np.deg2rad(kw.get("min_distance_deg", 0))
            while len(valid_dec) < num: #Keep on filling up the vlaid RA and DEC lists with values
                
                dec_vals = np.random.uniform(*dec_range, size = num_fake)
                dec_rad = np.arcsin(dec_vals)  # Random Dec points in radians

                # Mask to exclude patches whose edge falls below -40 deg
                valid_mask = dec_rad >= np.radians(-40 + radius)

                # #Ensuring dec values stay within valid celestial bounds
                # valid_mask &= (np.rad2deg(dec) + radius <= 90)
                # valid_mask &= (np.rad2deg(dec) - radius >= -90)

                dec = np.rad2deg(dec_rad)  # Random Dec points
                valid_dec = dec[valid_mask]

                num_fake +=200
                count += 1
                
                # print(f"{num_fake=}")
                # print(f"{len(valid_dec)=}")
                # print(f"{count=}")
                # print(f"_________________")
            
            ra = np.random.uniform(*ra_range, num_fake)  # Random RA points
            L = np.array(ra[:num]), np.array(valid_dec[:num])
            # print(f"{min(valid_dec)=}")
            # print(f"{max(valid_dec)=}")
            # import sys; sys.exit()
            print(f"10000 Random centre points made. Discarded {num_fake-num} points")
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
        
def get_mean_and_med_stats(x_vals, y_vals, bin_num, x_is_r_proj=True, absol=True, **kw):
    
    """
    x vals: Polar Angle or Separation Distance measurements
    y_vals: Usually going to be the Rotation Measures
    bin_num: (Probably same as variable "bin_num_from_main")
    x_is_r_proj: Set to True if x_vals is Separation Distance measurements
    """
    if absol: y_vals = np.abs(y_vals)

    # print(y_vals); 
    # print(f"{len(x_vals)=}")
    # print(f"{len(y_vals)=}")
    # print(f"{bin_num=}")
    
    #Calculating mean and median of RM values within patch
    bin_means, bin_edges_mean, binnumber = stats.binned_statistic(x_vals, y_vals, statistic='mean', bins=bin_num)
    bin_med, bin_edges_med, binnumber = stats.binned_statistic(x_vals, y_vals, statistic='median', bins=bin_num)

    #making bin_edges same as m31's bin edges:
    bin_edges_mean = bin_edges_mean_m31
    bin_edges_med = bin_edges_median_m31

    if x_is_r_proj: 
        bin_width_mean, bin_width_med = (bin_edges_mean[1] - bin_edges_mean[0]) , (bin_edges_med[1] - bin_edges_med[0])
        bin_centers_mean, bin_centers_med = bin_edges_mean[1:] - bin_width_mean / 2 , bin_edges_med[1:] - bin_width_med / 2
    else: #making all bins the same for PA vs RM
        bw_deg = kw["bin_width_degrees"]
        bin_width_mean, bin_width_med = bw_deg, bw_deg #Bin size is 15 degrees
        bin_centers_mean, bin_centers_med = np.arange(0,360+bw_deg, bw_deg)[1:]- bin_width_mean / 2, np.arange(0,360+bw_deg, bw_deg)[1:]- bin_width_med / 2

    #Standard error of the mean for error bars
    bin_std, bin_edges, binnumber = stats.binned_statistic(
        x_vals, y_vals,
        statistic=stats.sem,  # Standard Error of Mean
        bins=bin_num)

    if x_is_r_proj: #previously only used edges from median.
        #Converting to projected distance away from center of m31(in kpc)
        d_bin_centers = get_projected_d(bin_centers_mean*u.deg, d_m31).value
        return d_bin_centers, bin_means, bin_med, bin_std, bin_edges
    
    else:
        return bin_centers_mean, bin_centers_med, bin_means, bin_med, bin_std, bin_edges


def annuli_analysis_m31(rm_m31=rm_m31, err_m31=err_m31, save_plot=False):
    print("M31 Halo Annuli Histogram analysis and plotting have begun")

    def construct_and_plot_annuli(distance_flat, rm_vals_flat, rm_errs_flat, 
                                distance_flat_bg, rm_vals_flat_bg, rm_errs_flat_bg, 
                                  save_plot=save_plot):

        # Converting all to numpy lists if needbe
        L = (distance_flat, rm_vals_flat, rm_errs_flat,
            distance_flat_bg, rm_vals_flat_bg, rm_errs_flat_bg)
        (distance_flat, rm_vals_flat, rm_errs_flat,
            distance_flat_bg, rm_vals_flat_bg, rm_errs_flat_bg) = map(np.asarray, L)

        # #Adding in background region (Which has also been background subtracted)
        # distance_flat = np.concatenate([distance_flat,distance_flat_bg])
        # rm_vals_flat = np.concatenate([rm_vals_flat, rm_vals_flat_bg])
        # rm_errs_flat = np.concatenate([rm_errs_flat, rm_errs_flat_bg])
        

        def add_colorbar(ax, data, colormap, label, num_ticks=10):
            """
            Adds a colorbar to the given axis.
            
            Parameters:
            - ax: Matplotlib axis to attach the colorbar to.
            - data: Data for normalization (e.g., distances or annulus indices).
            - colormap: The colormap to use.
            - label: Label for the colorbar.
            - num_ticks: Number of ticks to generate for the colorbar.
            """
            norm = mcolors.Normalize(vmin=np.min(data), vmax=np.max(data))
            sm = cm.ScalarMappable(cmap=colormap, norm=norm)
            sm.set_array([])  # Required for colorbar

            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label(label, rotation=-90, labelpad=20)

            # Generate evenly spaced ticks
            custom_ticks = sorted(np.floor(np.linspace(np.min(data), np.max(data), num_ticks)))
            cbar.set_ticks(custom_ticks)
            cbar.ax.tick_params(labelsize=7)

            return cbar

        def plot_seaborn_hist(ax, rm_per_annulus, histbin, cmap_name, d_flat, annul_area):
            """
            Plots a stacked histogram using Seaborn with a colorbar indicating radial distance.

            Parameters:
            - ax: Matplotlib axis to plot on.
            - rm_per_annulus: Dictionary containing RM values per annulus.
            - histbin: Number of bins for histogram.
            - cmap_name: Colormap name for seaborn histogram.
            - d_flat: Array of radial distances.
            - annul_area: Area of discrete annulus (Should be of equal area) (ITS NOT !)
            """

            # ax.set_title(r"Discrete Annulus Area $\xi$" + f" = {annul_area:.2f} " + f"{annul_dist_type}" r"$^2$")
            
            #Flattenning rm data
            all_rm_values = np.concatenate(list(rm_per_annulus.values()))

            #Array of annulus indices for hue
            annulus_indices = np.concatenate([[idx] * len(v) for idx, v in enumerate(rm_per_annulus.values())])

            #Normalizing weights by the annulus area
            weights = np.concatenate([np.full(len(v), 1 / annul_area) for v in rm_per_annulus.values()])

            #Colormap object
            colormap = plt.get_cmap(cmap_name, len(rm_per_annulus))

            # Plot histogram
            sns.histplot(
                x = all_rm_values,
                hue=annulus_indices,
                bins=histbin,
                palette=colormap,
                element="step",      # Outline-only histogram (no fill)
                color="none",        # Ensure bars are fully transparent
                linewidth=0.04,      # Set edge thickness
                edgecolor="k",
                multiple="layer",
                weights=weights,
                ax=ax
            )

            add_colorbar(ax, data=d_flat, colormap=colormap, label="Radial Distance [kpc]", num_ticks=10)
            add_colorbar(ax, data=annulus_indices, colormap=colormap, label="Annulus Index", num_ticks=10)
            
            if ax.get_legend(): ax.get_legend().remove()
            ax.minorticks_on()

        def set_figure():
            fig, ax = plt.subplots(1, 1, figsize=(12, 6))  #1x1 subplot grid
            ax.set_xlim(-150, 150)
            ax.set_xlabel("RM [rad m$^{-2}$]")
            ax.set_ylabel("Counts / " + r"$\xi$" + f" [{annul_dist_type}" + r"$^{-2}$]")
            return fig, ax 

        #Binning edges based on distances
        bin_edges = np.linspace(distance_flat.min(), distance_flat.max(), bin_num_from_main)
        bin_edges_deg = get_sep_angle(bin_edges, d_m31)

        # Updated function calls with correct arguments
        [rm_per_annulus, _], annul_area, bin_edges = create_annuli_binning_structure(distance_flat, [rm_vals_flat, rm_errs_flat], bin_num_from_main+1)
        [rm_per_annulus_bg, _], annul_area_bg, _ = create_annuli_binning_structure(distance_flat_bg, [rm_vals_flat_bg, rm_errs_flat_bg], bin_num_from_main+1)

        annuli_to_plot = np.arange(0, bin_num_from_main + 1)
        annul_dist_type = "kpc"

        b_width = 8  # Width of histogram bars
        if args.overplot:
            fig, ax = set_figure()

        Counts = []
        histbin = 60  # Number of bins for a single histogram

        for bin_idx in annuli_to_plot:
            
            if not args.overplot: 
                fig, ax = set_figure() #Then plot individually

            if bin_idx in rm_per_annulus:
                if args.overplot: #Will be run only once then broken of the loop
                    if args.seaborn_hist:
                        plot_seaborn_hist(
                            ax=ax,
                            rm_per_annulus=rm_per_annulus, #Halo + BG
                            histbin=histbin,
                            cmap_name="viridis_r",
                            d_flat=distance_flat,
                            annul_area=annul_area[0] #Assuming all area is the same
                        )


                        # plot_seaborn_hist(
                        #     ax=ax,
                        #     rm_per_annulus=rm_per_annulus_bg, #BG ontop of previous seaborn hist
                        #     histbin=histbin,
                        #     cmap_name="magma",
                        #     d_flat=distance_flat_bg,
                        #     annul_area=annul_area_bg
                        # )

                        break #Seaborn can plot everything all at once
                    else: #USING BASIC MATPLOTLIB

                        # Flatten data
                        all_rm_values = np.concatenate(list(rm_per_annulus.values()))
                        all_rm_values_bg = np.concatenate(list(rm_per_annulus_bg.values()))
                        # weights = np.concatenate([np.full(len(v), 1 / annul_area) for v in rm_per_annulus.values()])

                        #Making one single annulus for BG region
                        all_rm_values_bg_flat = all_rm_values_bg
                        all_rm_values_bg_flat.flatten()

                        #Defining bins
                        bins = np.linspace(np.min(all_rm_values), np.max(all_rm_values), histbin-1)
                        # bins_bg = np.linspace(np.min(all_rm_values_bg), np.max(all_rm_values_bg), histbin-1)

                        #Colormap normalization
                        cmap = cm.get_cmap("tab20_r", len(rm_per_annulus))
                        norm = mcolors.Normalize(vmin=0, vmax=len(rm_per_annulus) - 1)
                        # hatch_patterns = ["//", "\\\\","--"]# "xx"]#, "--", "oo", "OO", "**", "..", "||", "++"]

                        def convert_kpc2_to_deg2(bin_edges, d_m31):
                            """
                            Convert annular bin areas from kpc² to deg² using steradian conversion.

                            Parameters:
                            -----------
                            bin_edges : array-like
                                The bin edges in kpc.
                            d_m31 : astropy.Quantity
                                Distance to M31 in kpc (must have length units).

                            Returns:
                            --------
                            area_deg2 : astropy.Quantity
                                The area of each annular bin in square degrees.
                            """
                            try:
                                bin_edges = np.asarray(bin_edges)

                                # Ensure bin_edges has an even number of elements
                                if len(bin_edges) % 2 != 0:
                                    bin_edges = np.insert(bin_edges, 0, 0)  # Insert 0 at the start to maintain even count

                                # Ensure d_m31 has proper length units
                                if not isinstance(d_m31, u.Quantity):
                                    raise ValueError("d_m31 must be an astropy Quantity with length units.")

                                # Compute angular separations (convert to radians)
                                theta_1 = (np.arctan(bin_edges[:-1] / d_m31.to(u.kpc).value) * u.deg).to(u.rad)  # Inner radii
                                theta_2 = (np.arctan(bin_edges[1:] / d_m31.to(u.kpc).value) * u.deg).to(u.rad)   # Outer radii

                                # # Compute area in square degrees #version 1
                                # area_deg2 = (theta_2**2 - theta_1**2) * (180 / np.pi)**2 #Already in degrees squared

                                # # Compute area in square degrees #version 2
                                # area_deg2 = 2*np.pi* ( np.cos(theta_1) - np.cos(theta_2)) * (180 / np.pi)**2 #Already in degrees squared

                                # Compute area in square degrees #version 2
                                area_deg2 = 4*np.pi* ( np.sin(theta_2/2)**2 - np.sin(theta_1/2)**2) * (180 / np.pi)**2 #Already in degrees squared

                                return area_deg2

                            except AttributeError as e:
                                print(f"AttributeError: {e}")
                                print("Ensure that bin_edges is a numerical array and d_m31 has proper units (e.g., kpc).")
                                return None
                            except ValueError as e:
                                print(f"ValueError: {e}")
                                return None

                        #Annulus areas within CGM in degrees
                        annul_area_deg = convert_kpc2_to_deg2(bin_edges, d_m31) 
                        
                        # # Ensure rm_per_annulus keys match the expected number of bins
                        # if len(rm_per_annulus) > len(annul_area_deg):
                        #     print("Warning: rm_per_annulus has extra bins, adjusting...")
                        #     rm_per_annulus = {k: v for k, v in list(rm_per_annulus.items())[:len(annul_area_deg)]}

                        # Plot histograms with "step" type for outlines
                        for idx, (label, values) in enumerate(rm_per_annulus.items()):

                            color = cmap(norm(idx))  # Get color from viridis colormap
                            # hatch = hatch_patterns[idx % len(hatch_patterns)]  # Cycle through hatching styles

                            #HISTOGRAM OF RM COLLECTED PER ANNULI WITHIN CGM OF M31
                            plt.hist(values, 
                                    bins=bins, 
                                    weights=np.full(len(values), 1 / annul_area_deg[idx-1]),
                                    histtype="step", 
                                    color=color,
                                    facecolor="none",
                                    # hatch=hatch,
                                    label=f"{label}")
                            
                            # plt.title(r"Halo of M31 | Annuli Area $\xi$" + f" = {annul_area_deg.value:.2f} " + r"$[\circ]^{2}$")

                        # annul_area_entire_bg = np.pi * (get_projected_d(cutoff, d_m31).value **2 - R_vir.value **2)
                        annul_area_entire_bg_deg = 4*np.pi* ( np.cos(cutoff.to(u.rad)) - np.cos((L_m31*u.deg).to(u.rad))) #* (180/np.pi)**2 #Annuli of Backgorund in degrees^2

                        # # HISTOGRAM OF RM COLLECTED FOR ONE ANNULI WITHIN M31'S CGM'S BG
                        # plt.hist(all_rm_values_bg_flat, 
                        #         bins=bins_bg, 
                        #         weights=np.full(len(all_rm_values_bg_flat), 1 / annul_area_entire_bg),
                        #         histtype="step", 
                        #         linestyle="--",
                        #         facecolor="none",
                        #         hatch="value",
                        #         label="Background (Area = " f"{annul_area_entire_bg:.2f} " + r"$kpc^{2}$)",
                        #         color="k")
                            
                        plt.legend(loc="upper center", ncol=7, framealpha=0, bbox_to_anchor = (0.5, 1.25))
                        
                        if save_plot: #Save individual plots
                            path = curr_dir_path() + "Results/"
                            plt.savefig(f"{path}m31_annuli_plot_{bin_idx}.png", dpi=600, bbox_inches="tight")
                            plt.close()
                    
                        # Create side-by-side subplots
                        fig, ax = plt.subplots(1, 3, figsize=(12, 5))

                        # Plot 1: Total counts of RM values per annulus index
                        # rm_per_annulus["BG"] = all_rm_values_bg_flat #Including the background
                        ax[0].plot(rm_per_annulus.keys(), [len(v) for v in rm_per_annulus.values()], 
                                   marker=".", linestyle="", color="b")
                        ax[0].set_xlabel("Annulus Index")
                        ax[0].set_ylabel("Total Count of RM Values")
                        ax[0].grid(True)

                        # annul_area_deg = np.concatenate([annul_area_deg, [annul_area_entire_bg_deg]])
                        # bin_edges_deg = np.concatenate([bin_edges_deg, [cutoff]])
                        # bin_edges = np.concatenate([bin_edges, [get_projected_d(cutoff, d_m31).value]])

                        # Plot 2: Annulus area vs. radius
                        print(f"{len(bin_edges_deg)=}")
                        print(f"{len(annul_area_deg)=}")

                        # print(f"{rm_per_annulus=}")
                        ax[1].plot(bin_edges_deg[1:], annul_area_deg, marker=".", linestyle="", color="r")
                        # ax[1].plot(bin_edges[2:], annul_area_deg, marker="*", linestyle="--", color="g", label="upper bound")
                        ax[1].set_xlabel("Radius (degrees)")
                        ax[1].set_ylabel(r"Annulus Area [$deg^2$]")
                        ax[1].grid(True)

                        # Plot 3:  Total counts of RM values vs. radius
                        ax[2].plot(bin_edges_deg, [len(v) for v in rm_per_annulus.values()], marker=".", linestyle="", color="g")
                        ax[2].set_xlabel("Radius (degrees)")
                        ax[2].set_ylabel("Total Count of RM Values")
                        ax[2].grid(True)

                        # Adjust layout
                        plt.tight_layout()

                        # Save figure
                        plt.savefig(f"{curr_dir_path()}Results/combined_annuli_plots.png", dpi=600, bbox_inches="tight")
                        plt.close()

                        break

                scale = 1/annul_area
                counts, _, patches = ax.hist(rm_per_annulus[bin_idx], bins=histbin, alpha=0.1, color="k")

                for p in patches:
                    p.set_height(p.get_height() *scale)
                    p.set_width(b_width)
                    ax.set_ylim(0, np.max(counts) * 1.05) #ylimit for mean subplot

                mark_m31_sem_vals_on_annulus_hist(ax, bin_means_past_rvir, bin_idx, bin_std_past_rvir, label_prefix="m31")
                
                fig.suptitle(f"Annulus ({bin_edges[bin_idx-1]:.2f} < r < {bin_edges[bin_idx]:.2f}) {annul_dist_type}")

                if save_plot: #Save individual plots
                    path = curr_dir_path() + "Results/"
                    plt.savefig(f"{path}m31_annuli_plot_{bin_idx}.png", dpi=600, bbox_inches="tight")
                    plt.clf()
                else:
                    plt.show()

        
        if save_plot:
            if not args.m31_annuli_anal:
                plt.close()  #Deleting figure
                print(f"All images saved to {path}")
            else:
                # ax.set_ylim(0, max(Counts) * 1.1)
                plt.xlim(-350, 200)
                path = curr_dir_path() + "Results/"
                plt.savefig(f"{path}m31_annuli_overplot.png", dpi=600, bbox_inches="tight")
                print(f"M31 Overplot has been saved to {path}")
        else:
            ax.set_ylim(0, max(Counts) * 1.1)
            plt.show()

    #Running annuli analysis just for M31 Halo
    m31_distances = tuple(map(get_projected_d, m31_sep_Rvir, [d_m31]*len(m31_sep_Rvir))) #from seprated distance in degrees - within CGM of M31 - to kpc
    m31_distances = list(map(lambda m31_d: m31_d.value, m31_distances))
    m31_distances_bg = tuple(map(get_projected_d, m31_sep_bg, [d_m31]*len(m31_sep_Rvir))) #from seprated distance in degrees - of the CGM Background of M31 - to kpc
    m31_distances_bg = list(map(lambda m31_d: m31_d.value, m31_distances_bg))

    construct_and_plot_annuli(m31_distances, rm_m31, err_m31,
                            m31_distances_bg, rm_bg, err_bg)
    
def mark_m31_sem_vals_on_annulus_hist(ax, b_m, bin_idx, std, label_prefix="m31"):
    """Marks M31-related values and fills +-1 sigma region on a histogram axis."""

    ax.axvline(x=b_m[bin_idx-1], 
               label=f"{label_prefix} = {b_m[bin_idx-1]:.2f}", 
               color='k', linestyle='--', linewidth=.8)

    ax.fill_betweenx(y=np.linspace(0, 100, 100),
                     x1=(b_m[bin_idx-1] - std[bin_idx-1]),
                     x2=(b_m[bin_idx-1] + std[bin_idx-1]),
                     color='k', alpha=0.2, edgecolor="none",
                     label=r"$\sigma \approx {:.2f}$".format(std[bin_idx-1]))

    ax.legend()

def annuli_analysis_random(all_means_corr, all_medians_corr, save_plot=False, stack_indiv_patch=False): 

    """
    This function conducts the annulus analysis of all the random patches in the sky.
    For assessment done for M31 see the function "annuli_analysis_m31".

    stack_indiv_patch - Taking all patches of the sky and stacking them on top of each other
                        to analyse the ultimate change in RM over a given annulus; along the stack
    """
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
            rm_per_annulus_mean[bin_idx].append(rm_vals_flat_mean[i])  # Append RM values for mean to corresponding bin
            rm_per_annulus_median[bin_idx].append(rm_vals_flat_median[i])  # Append RM values for median to corresponding bin

        #Converting lists to numpy arrays for easier handling
        rm_per_annulus_mean = {k: np.array(v) for k, v in rm_per_annulus_mean.items() if len(v) > 0}
        rm_per_annulus_median = {k: np.array(v) for k, v in rm_per_annulus_median.items() if len(v) > 0}

        if stack_indiv_patch: 
            annuli_to_plot = np.arange(0, round(cutoff.value)+1) # Adjusting annuli ranges based on bin edges (in degrees)
        else: 
            annuli_to_plot = np.arange(0, bin_num_from_main+1)  # Adjusting annuli ranges based on bin edges (in kpc)

        annul_dist_type = "deg" if stack_indiv_patch else "kpc"
        b_m_1 = bin_means_m31 if stack_indiv_patch else bin_means_past_rvir
        b_m_2 = bin_med_m31 if stack_indiv_patch else bin_meds_past_rvir
        std = bin_std_m31 if stack_indiv_patch else bin_std_past_rvir #assuming same standard deviation

        b_width = .8 #width of bars in the overplot
        if args.overplot:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))  #1x2 subplot grid
            Counts = {"mean":[], "med":[]} #To be able to retirve highest count value for histograms

        x_axis_label = r" [rad m$^{-2}$]"
        histbin = 50 #number of bins per histogram's RM-axis (x-axis)
        #Looping through annuli
        for bin_idx in annuli_to_plot:

            if not args.overplot: #Then plot individually
                fig, axes = plt.subplots(1, 2, figsize=(12, 6))  #1x2 subplot grid
                
            if bin_idx in rm_per_annulus_mean and bin_idx in rm_per_annulus_median:

                xlim, ylim = (-100, 200), (0, 17)
                #Plotting for "Mean" subplot (left side)
                counts, _, patches_mean = axes[0].hist(rm_per_annulus_mean[bin_idx], bins=histbin, alpha=0.1, color="k")
                axes[0].set_title("Mean")
                axes[0].set_xlabel("RM " +  x_axis_label)
                axes[0].set_ylabel("Counts" + r"/$\xi$" + f" [{annul_dist_type}"+ r"$^{-2}$]")

                #dividing Counts of Mean by Annulus Area
                for p in patches_mean: 
                    p.set_height(p.get_height() / annul_area) #Crucial for accurate diviiding counts by area
                    p.set_width(b_width) #Decreasing width of bars
                counts /= annul_area

                if args.overplot: 
                    Counts["mean"].append(np.max(counts))
                else: 
                    axes[0].set_ylim(0, np.max(counts) * 1.1)

                # axes[0].set_ylim(*ylim)
                axes[0].set_xlim(-75, 125)
                
                #Plotting for "Median" subplot (right side)
                counts, _, patches_med = axes[1].hist(rm_per_annulus_median[bin_idx], bins=histbin, alpha=0.1, color="k")
                axes[1].set_title("Median")
                axes[1].set_xlabel("RM " +  x_axis_label)
                # axes[1].set_ylabel("Counts")

                #dividing Counts of Median by Annulus Area
                for p in patches_med: 
                    p.set_height(p.get_height() / annul_area) #Crucial for accurate diviiding counts by area
                    p.set_width(b_width) #Decreasing width of bars
                counts /= annul_area

                if args.overplot: 
                    Counts["med"].append(np.max(counts))
                else: 
                    axes[1].set_ylim(0, np.max(counts) * 1.1)

                # axes[1].set_ylim(*ylim)
                axes[1].set_xlim(-75, 100)

                if not args.overplot: #Add vertical line (and std) to individual histogram
                    mark_m31_sem_vals_on_annulus_hist(axes[0], b_m_1, bin_idx, std, label_prefix="m31")
                    mark_m31_sem_vals_on_annulus_hist(axes[1], b_m_2, bin_idx, std, label_prefix="m31")

                
                if not args.overplot:
                    fig.suptitle(f" ({bin_edges[bin_idx-1]:.2f}" + r" $< r_{proj} <$ " + f"{bin_edges[bin_idx]:.2f}) {annul_dist_type}" + r"  |  Area $\xi$" + f" = {annul_area:.2f} " + f"{annul_dist_type}" r"$^2$")
                # else:
                    # fig.suptitle(r"Discrete Annulus Area $\xi$" + f" = {annul_area:.2f} " + f"{annul_dist_type}" r"$^2$")

                if not args.overplot: #Then plot individually
                    #Saving or displaying
                    if save_plot:
                        path = curr_dir_path() + "Results/"
                        plt.savefig(f"{path}annuli_plot_{bin_idx}.png", dpi=600, bbox_inches="tight")
                        plt.clf()  # clearing the figure (not deleting it)
                    else:
                        plt.show()

                # if bin_idx == 3: break #Testing out one (or a few) plots

        if args.overplot:

            #Give correct y-axis limits from maximum histogram
            axes[0].set_ylim(0, max(Counts["mean"])* 1.1)
            axes[1].set_ylim(0, max(Counts["med"])* 1.1)

            if save_plot:#Finally Saving the overplots
                path = curr_dir_path() + "Results/"
                plt.savefig(f"{path}annuli_overplot.png", dpi=600, bbox_inches="tight")
            else:
                plt.show() #Otherwise show the overplot

        if save_plot: #Since its not easy to make plots interactively on ilifu
            plt.close()  #Deleting the figure to clear memory
            print(f"All images saved to {path}")
    
        if args.annuli_video: #Saving a video if needbe
            image_files = sorted(glob.glob(f"{path}annuli_plot_*.png"))

            fig, ax = plt.subplots(figsize=(5.3,3.2))

            #Removing everything unnecessary from outer figure
            ax.set_xticks([])  
            ax.set_yticks([])  
            ax.set_frame_on(False) #Including boarders.

            img = plt.imshow(plt.imread(image_files[0]))

            def update(frame):
                img.set_array(plt.imread(image_files[frame]))
                return [img]

            print(f"Saving of annuli_video.mp4 in {path} ...")

            ani = animation.FuncAnimation(fig, update, frames=len(image_files), interval=500)
            ani.save(f"{path}annuli_video.mp4", fps=2, writer="ffmpeg", dpi=400)
            
            print(f"Video saved to {path}annuli_video.mp4")

    # Stack all patches together without any mean analysis
    if stack_indiv_patch:

        # Flattening the lists for easier computation (These are for random patches all over the sky)
        flat_sep_vals = np.concatenate([patch.value for patch in CGM_RM_coords_sep])  # Separation distances (degrees)
        flat_rm_vals_mean = np.concatenate([CGM_RM_coords_per_patch_mean for CGM_RM_coords_per_patch_mean in CGM_RM_coords_per_patch])  # Corresponding RM values for mean
        flat_rm_vals_median = np.concatenate([CGM_RM_coords_per_patch_median for CGM_RM_coords_per_patch_median in CGM_RM_coords_per_patch])  # Corresponding RM values for median
        construct_and_plot_annuli(flat_sep_vals, flat_rm_vals_mean, flat_rm_vals_median)

    # Using binned values for MEAN and MEDIAN (As initially discussed with DJ and Russ)
    else:
        D_Bin_centers = np.concatenate([bin_centers for bin_centers in all_d_bin_centers])
        Avg_Means = np.concatenate([avg_mean for avg_mean in all_means_corr]) 
        Avg_Medians = np.concatenate([avg_med for avg_med in all_medians_corr]) 
        construct_and_plot_annuli(D_Bin_centers, Avg_Means, Avg_Medians)

def plot_binned_gal_lat(B, RM, bin_edges, save_plot=False):
    """Raw RM against galactic Latitude - with pretty colorbar :)"""

    # Number of bins
    n_bins = bin_num_from_main + 1
    
    # Create bin edges as in Code 1
    radial_lims = (0, 300)
    bin_edges = np.linspace(*radial_lims, n_bins)
    bin_width = bin_edges[1] - bin_edges[0]
    bin_centers = bin_edges[:-1] + bin_width / 2  # Bin centers

    # Create a discrete version of the colormap
    base_cmap = plt.get_cmap("jet")
    cmap = base_cmap(np.linspace(0, 1, n_bins))  # n_bins distinct colors
    cmap_discrete = plt.matplotlib.colors.ListedColormap(cmap)

    # Normalize the colorbar using the bin edges
    norm_cb = plt.matplotlib.colors.BoundaryNorm(boundaries=bin_edges, ncolors=n_bins-1)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    # Plot background (black)
    sctt = ax.scatter(rm_pos_gal_lat_bg, rm_bg, marker="x", alpha=.4, s=10, c="k", label="Background", linewidths=.5)

    # Plot each bin with a distinct color
    for bin_idx in range(1, len(RM)+1):
        if bin_idx in RM:
            # Using bin_centers as the colormap index
            color_idx = np.argmin(np.abs(bin_centers - bin_centers[bin_idx-1]))
            ax.scatter(B[bin_idx], RM[bin_idx], marker=".", alpha=1, s=24, edgecolors='none', color=cmap_discrete(color_idx))

    # Colorbar using the same discrete colormap
    sm = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm_cb)
    sm.set_array([])

    # Set the ticks at bin_centers and make sure the colorbar has proper boundaries
    cbar = plt.colorbar(sm, ax=ax, boundaries=bin_edges)
    cbar.set_ticks(bin_centers)  # Ensure ticks align with bin_centers

    cbar.set_label("Radial Distance [kpc]", rotation=-90, labelpad=20)
    cbar.ax.tick_params(size=3, which='both', length=0)

    # Vertical line at M31's latitude
    ax.axvline(x=m31_lat.deg, linestyle="--", color="k", label=r"$b_{M31}$")

    ax.grid(True)
    ax.legend()

    apply_plot_attributes(xlim=(-51,10), ylim=(-300,200), 
                          xlabel="Galactic Latitude (b)" + r"[$^{\circ}$] ",
                          bbox=(0.5,1.1), framealpha=1, ncols=2)

    if save_plot:
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}Galactic_Lat_bin_plots_kpc.png", dpi=600, bbox_inches="tight")
        print(f"Galactic Latitude vs RM plot saved in {path}")

        plt.close()
    else:
        plt.show()

def fit_and_plot_cosine(x, y, guesses):
    """
    Fits a cosine function to the provided x and y data, then plots the results.

    Parameters:
    x : array-like
        The x values of the data.
    y : array-like
        The y values of the data.
    guesses : tuple
        Initial guesses for the fitting parameters (A, B, C, D).
    """
    def cos_wave(x, A, B, C, D):
        return A * np.cos(np.deg2rad(B) * x + np.deg2rad(C)) + D
    try:
        # Fit the curve
        popt, _ = curve_fit(cos_wave, xdata=x, ydata=y, p0=guesses)
        
        # Generate fitted data
        x_fit = np.linspace(min(x), max(x), 1000)
        y_fit = cos_wave(x_fit, *popt)
        
        # Label function
        label = lambda vals: rf"RM = ${vals[0]:.2f} \cdot \cos({vals[1]:.2f} \cdot \theta {'+' if vals[2]>0 else ''} {vals[2]:.2f}^{{\circ}}) {'+' if vals[3]>0 else ''} {vals[3]:.2f} \text{{ rad }}m^{{-2}}$"
        
        # Plotting the fit
        # plt.scatter(x, y, color='r', label='Data')
        plt.plot(x_fit, y_fit, color='k', linestyle='--', label=label(popt))
        plt.legend()

    except RuntimeError:
        print("Curve fitting failed")


def plot_binned_PA(PA, RM, bin_edges, save_plot=False):

    """
    Plots Rotation Measure (RM) as a function of polar angle/Polar_Angle (in the anticlockwise direction) according to binned distance from center of m31
    """

    PA_bin_width = 30 #in degrees
    radial_bin_width = int(300/bin_num_from_main) #in kpc

    #For distances
    bin_width = bin_edges[1] - bin_edges[0]
    bin_centers = bin_edges - bin_width / 2

    bin_centers_mean_GMMStats = []
    bin_centers_med_GMMStats = []
    bin_means_GMMStats = []
    bin_med_GMMStats = []
    bin_std_GMMStats = []
    bin_edges_GMMStats = []
    
    #from Dictionaries to list of np.ndarray siblists (if need be)
    RM_new = [rm for rm in RM.values()] if isinstance(RM, dict) else RM
    PA_new = [pa for pa in PA.values()] if isinstance(PA, dict) else PA
    
    for i in range(len(RM_new)): #Assuming ranges of RA and PA are equivalent
        
        if not (RM_new[i].shape == ()):  #Checking if its not empty or filled with NaNs
            #Note , this inner binning is for a clean 30 bins due to 360 degrees PA
            bin_cent_mean, bin_cent_med, bin_mean, bin_med, bin_std, bin_edges = get_mean_and_med_stats(PA_new[i],RM_new[i], bin_num=int(360/PA_bin_width) #360 divided by "bin_width_degrees"
                                                                                                        , x_is_r_proj=False, bin_width_degrees=PA_bin_width, absol=False)
            bin_centers_mean_GMMStats.append(bin_cent_mean) #For x axis
            bin_centers_med_GMMStats.append(bin_cent_med) #For x axis
            bin_means_GMMStats.append(bin_mean) #For y-axis
            bin_med_GMMStats.append(bin_med) #For y-axis
            bin_std_GMMStats.append(bin_std)
            bin_edges_GMMStats.append(bin_edges) #Havent been used for anything yet....

    # # This code below was hashed out since it doesnt make a difference in cleaning up the messy looking plot
    # #sorting in terms of Polar angle (x-axis) (respectively for mean and median)
    # sorted_idx_mean = [np.argsort(lst) for lst in bin_centers_mean_GMMStats]
    # sorted_idx_med = [np.argsort(lst) for lst in bin_centers_med_GMMStats]
    # bin_centers_mean_GMMStats = [lst[idx] for lst, idx in zip(bin_centers_mean_GMMStats,sorted_idx_mean)]
    # bin_centers_med_GMMStats = [lst[idx] for lst, idx in zip(bin_centers_med_GMMStats,sorted_idx_med)]
    # bin_means_GMMStats = [lst[idx] for lst, idx in zip(bin_means_GMMStats, sorted_idx_mean)]
    # bin_med_GMMStats = [lst[idx] for lst, idx in zip(bin_med_GMMStats, sorted_idx_med)]

    L = (bin_centers_mean_GMMStats, bin_centers_med_GMMStats, 
         bin_means_GMMStats, bin_med_GMMStats, 
         bin_std_GMMStats, bin_edges_GMMStats)
    
    (bin_centers_mean_GMMStats, bin_centers_med_GMMStats,
     bin_means_GMMStats, bin_med_GMMStats, 
     bin_std_GMMStats, bin_edges_GMMStats) = map(np.asarray, L)

    # plt.figure(figsize=(16,8))

    # # plt.errorbar(bin_centers_mean_GMMStats.flatten(), bin_means_GMMStats.flatten(), yerr = bin_std_GMMStats.flatten(), 
    # #             color="blue", fmt = '.', alpha=.6, label="Mean")

    # # # Flatten all arrays and sorting them based on "bin_centers_flat"
    bin_centers_flat = bin_centers_med_GMMStats.flatten()
    bin_med_flat = bin_med_GMMStats.flatten()
    bin_std_flat = bin_std_GMMStats.flatten()
    # # sort_idx = np.argsort(bin_centers_flat)
    # # bin_centers_sorted = bin_centers_flat[sort_idx]
    # # bin_med_sorted = bin_med_flat[sort_idx]
    # # bin_std_sorted = bin_std_flat[sort_idx]
    # #removing all non-values for the plot to have connections:
    # # valid_mask = ~np.isnan(bin_centers_sorted) & ~np.isnan(bin_med_sorted) & ~np.isnan(bin_std_sorted)
    no_nan_mask = ~np.isnan(bin_centers_flat) & ~np.isnan(bin_med_flat) & ~np.isnan(bin_std_flat)

    # # bin_centers_sorted = bin_centers_sorted[valid_mask]
    # # bin_med_sorted = bin_med_sorted[valid_mask]
    # # bin_std_sorted = bin_std_sorted[valid_mask]

    # # plt.errorbar(bin_centers_sorted, bin_med_sorted, yerr=bin_std_sorted, 
    # #             ecolor='r', marker="", linestyle="-", capsize=0, alpha=.5, elinewidth=1)#, label="Median")
    
    # plt.errorbar(bin_centers_med_GMMStats.flatten(), bin_med_GMMStats.flatten(), yerr=bin_std_GMMStats.flatten(), 
    #             ecolor='r', marker="", linestyle="-", capsize=0, alpha=.5, elinewidth=1)#, label="Median")
    # # print(f"{bin_centers_mean_GMMStats.flatten()}") #shows that there are a few nan vlaues for std
    # # plt.scatter(bin_centers_mean_GMMStats, bin_means_GMMStats, label="Mean")
    plt.plot(bin_centers_med_GMMStats.flatten(), bin_med_GMMStats.flatten(), linestyle="-")
    
    plt.title("Median")
    apply_plot_attributes(push_title_up=1.1, leg=False, xlim=(0,360))

    if save_plot:
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}Polar_Angle_bin_plots_median.png", dpi=600, bbox_inches="tight")
    else:
        plt.show()
    
    #This section plots the Folded and fitted PA_vs + RM plot
    #_____________________________________________________________________________
    plt.figure(figsize=(16,8))

    
    bin_centers_flat_folded = bin_centers_flat%180
    # print(bin_med_flat_folded); import sys; sys.exit()
    plt.scatter(bin_centers_flat_folded, bin_med_flat)
    # plt.errorbar(bin_centers_med_GMMStats.flatten()%180, bin_med_GMMStats.flatten(), yerr=bin_std_GMMStats.flatten(), 
    #             ecolor='r', marker="", linestyle="-", capsize=0, alpha=.5, elinewidth=1)#, label="Median")

    fit_and_plot_cosine(x=bin_centers_flat_folded[no_nan_mask], 
                        y=bin_med_flat[no_nan_mask], 
                        guesses = (50, 
                                   3, 
                                   10, 
                                   np.mean(bin_med_flat[no_nan_mask])))

    plt.legend()
    plt.title(f"Median (bw_PA={PA_bin_width}" + r"$^{\circ}$" + f", bw_radial_proj = {radial_bin_width}"+r" $kpc$)")
    apply_plot_attributes(push_title_up=1.1, leg=False)
    
    if save_plot:
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}Polar_Angle_bin_plots_median_bwPA_{PA_bin_width}_bwrad_{radial_bin_width}kpc.png", dpi=600, bbox_inches="tight")
    else:
        plt.show()
    
    #This section just plots the Raw RM against Polar_Angle - with pretty colorbar :)
    #____________________________________________________________________________

    colormap_bg_combos = [
    # ("brg", "k"),
    # ("CMRmap", "#73e831"),
    # ("coolwarm", "k"),
    # ("gist_ncar", "k"),
    # ("gist_rainbow", "k"),
    # ("gnuplot2", "#73e831"),
    # ("hsv", "k"),
    ("jet", "k"),
    # ("nipy_spectral", "k"),
    # ("tab20", "#73e831"),
    # ("tab20", "k")
]
    
    for cmap_name, bg_color in colormap_bg_combos:
        print(f"Plotting with colormap_bg_color: {cmap_name}_{bg_color}")
        
        # Number of bins
        n_bins = bin_num_from_main + 1

        # Create bin edges as in Code 1
        radial_lims = (0, 300)
        bin_edges = np.linspace(*radial_lims, n_bins)
        bin_width = bin_edges[1] - bin_edges[0]
        bin_centers = bin_edges[:-1] + bin_width / 2  # Bin centers
        # print(f'{len(bin_centers)=}')
        # print(f'{RM.keys()}'); import sys; sys.exit()

        cmap_discrete = get_discrete_colors(radial_lims, n_bins, cmap_name)

        #Normalizing the colorbar using the bin edges
        norm_cb = plt.matplotlib.colors.BoundaryNorm(boundaries=bin_edges, ncolors=n_bins-1)

        fig, ax = plt.subplots(1, 1, figsize=(16, 8))

        #Background Points
        ax.scatter(PA_rm_deg_bg, rm_bg, marker="x", alpha=.4, s=15, c=bg_color, label="Background", linewidths=.5)

        #Plotting each bin with a distinct color
        # counter = 0
        # print(f"{[k for k, _ in RM.items()]=}")
        # print(f"{len(RM)=}")
        for bin_idx in range(1, len(RM)+1):
            if bin_idx in RM:
                # Using bin_centers as the colormap index
                color_idx = np.argmin(np.abs(bin_centers - bin_centers[bin_idx-1]))  # Match index to bin_centers
                ax.scatter(PA[bin_idx], RM[bin_idx], marker=".", alpha=1, s=12, color=cmap_discrete(color_idx))
            #     counter += 1
            # print(counter)
        ax.set_xticks(np.arange(0, 361, 30))

        # Colorbar using the same discrete colormap
        sm = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm_cb)
        sm.set_array([])

        # Set the ticks at bin_centers and make sure the colorbar has proper boundaries
        cbar = plt.colorbar(sm, ax=ax, boundaries=bin_edges)
        cbar.set_ticks(bin_centers)  # Ensure ticks align with bin_centers

        cbar.set_label("Radial Distance [kpc]", rotation=-90, labelpad=20)
        cbar.ax.tick_params(size=3, which='both', length=0)

        # Apply plot attributes (title, axis limits, etc.)
        apply_plot_attributes(push_title_up=1.1, xlim=(0, 360), ylim=(-300, 300))

        ax.legend(bbox_to_anchor=(1.05, 1.1), loc="upper left", framealpha=0)


        if save_plot:
            path = curr_dir_path() + "Results/"
            filename = f"{path}Polar_Angle_bin_colormap_{cmap_name}_{bg_color}_bwPA_{PA_bin_width}_bwrad_{radial_bin_width}kpc.png"
            plt.savefig(filename, dpi=600, bbox_inches="tight")
            plt.close()
            
            # break
        # break
    if save_plot:
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}Polar_Angle_bin_plots_bwPA_{PA_bin_width}_bwrad_{radial_bin_width}kpc.png", dpi=600, bbox_inches="tight")
        print("Saving the values from RM vs PA to text files")
        # Saving Polar Angle vs RM data
        # Saving Polar Angle vs RM data
        with open(f"{path}Polar_Angle_vs_RM_data.txt", "w") as f:
            f.write("# PA_bin\tRadialBinCenter[kpc]\tPA[deg]\tRM[rad/m^2]\n")
            for bin_idx in PA:
                # Ensure bin_idx is within the length of bin_centers
                if bin_idx < len(bin_centers) + 1:
                    radial_bin_center = bin_centers[bin_idx -1]
                else:
                    # Falling back if out of range
                    # (If this appears in txt file recheck code until it doesnt)
                    radial_bin_center = -999 

                for pa_val, rm_val in zip(PA[bin_idx], RM[bin_idx]):
                    f.write(f"{bin_idx}\t{radial_bin_center:.2f}\t{pa_val:.4f}\t{rm_val:.4f}\n")

        #Saving background region too
        with open(f"{path}Polar_Angle_vs_RM_background.txt", "w") as f:
            f.write("# PA[deg]\tRM[rad/m^2]\n")
            for pa_bg, rm_bg_val in zip(PA_rm_deg_bg, rm_bg):
                f.write(f"{pa_bg:.4f}\t{rm_bg_val:.4f}\n")

        plt.close()
    else:
        plt.show()

    # #this section plots as above but for individual subplots (of RM vs Polar_Angle) per radial bin
    # #_________________________________________________________________________________________
    
    # for bin_idx in range(1,len(RM)+1):
    #     fig, ax = plt.subplots(figsize=(12, 6))

    #     if bin_idx in RM:
    #         ax.scatter(PA[bin_idx], RM[bin_idx], marker=".", alpha=0.7, s=12, color="k")

    #         apply_plot_attributes(ax, leg=False, xlim=None, ylim=None)  # Pass 'ax' to ensure correct application

    #         if save_plot:
    #             path = curr_dir_path() + "Results/"
    #             plt.savefig(f"{path}Polar_Angle_bin_plots_bwPA_{PA_bin_width}_bwrad_{radial_bin_width}kpc_{bin_idx}.png", dpi=600, bbox_inches="tight")
    #             plt.close(fig)  # Close the figure after saving
    #         else:
    #             plt.show()

    # Ensure all figures are closed after the loop
    plt.close("all")  


def plot_indidividual_patch_stats(ax, d_bin_centers, bin_mean, bin_med, bin_std):
    # Plot the mean and median with error bars
    ax.errorbar(d_bin_centers, np.absolute(bin_mean), yerr=bin_std, fmt='k.-', alpha=.4)
    ax.errorbar(d_bin_centers, np.absolute(bin_med), yerr=bin_std, fmt='g.-', alpha=.4, capsize=2)

def plot_dispersion(bin_num, sep_vals, rm, centers, save_plot=False, **kw):

    plt.figure(figsize = (10, 6))

    #Dispersion of RM values in each bin (Standard Error of the Means)
    bin_std, _, _ = stats.binned_statistic(sep_vals, rm, statistic=stats.sem, bins = bin_num)
    plt.plot(centers, bin_std, "o", label="M31", color="blue")
    if kw["bin_std_random"] is not None: plt.plot(centers, kw["bin_std_random"], "o", label=r"Random $R_{vir}$", color="r")
    plt.xlabel('Projected Separation [kpc]')
    plt.ylabel('$\sigma_{RM}$ \n $[\mathrm{rad} \ \mathrm{m}^{-2}]$', rotation='horizontal', labelpad=60, fontsize=15)

    plt.grid(True)
    plt.tight_layout()

    # #Plotting the curve_fit
    # coefficients = np.polyfit(d_bin_centers, bin_std, 3)
    # fit_line = np.poly1d(coefficients)
    # plt.plot(d_bin_centers.value, fit_line(d_bin_centers.value), 
    #             color = 'orange', linestyle = "--")
    
    plt.title("Dispersion of RM values in each bin")

    #Mentioned on page 841 of https://doi.org/10.1093/mnras/stad2811
    x_values = np.linspace(0, 296, 1000)
    plt.axhline(xmin=0, xmax = np.max(x_values), y=6, linestyle="--",
                label="(Böckmann et al. 2023)" #"Observed $\sigma_{RM}$ indepenedent of Galactic Latitude"
                )

    plt.legend(fontsize = 12, loc = 'upper center', bbox_to_anchor = (0.5, 1.2),
                framealpha = 0, ncols =3)
    
    if save_plot: 
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}" + f"Dispersion{'_Mean' if args.mean else ('_Median' if args.median else '')}.png", dpi=600, bbox_inches="tight")#Saving as image
        print(f"Dispersion image saved in " + f"{path}")
        plt.clf() #clearing the figure (not deleting it)
    else:
        plt.show()
        plt.close() #Deleting the figure to clear memory

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

#Individiual background correction (Old Method. New method -s better - Spline)
def indiv_bg_corr(arr, bin_cent, absol=True):
    """
    bin_cent - Center of bin of projected distance or just projected distance relative to the RM values given in 'arr'
    arr - RM arrays . RM separation distances in degrees
    """

    # If arr and bin_cent are lists of arrays, process each separately
    if isinstance(arr, list) and isinstance(bin_cent, list):
        return [indiv_bg_corr(a, b, absol) for a, b in zip(arr, bin_cent)]
    
    if hasattr(bin_cent, "unit"):
        bin_cent = bin_cent.value
    # print(f"{np.asarray(bin_cent).shape=}")
    arr = np.asarray(arr)  # Ensure arr is also an array (important for indexing later)
    
    try:
        arr_bg = np.where(bin_cent > 300, arr, 0)  # Fill all values within virial radius with 0
    except TypeError:  #'>' not supported between instances of 'list' and 'int'
        bin_cent = np.asarray(bin_cent)
        print(f"{bin_cent[0]=}"); import sys; sys.exit()
        arr_bg = np.where(np.array(bin_cent, dtype=float) > 300, arr, 0)

    # Ensure arr_bg is at least 1D
    arr_bg = np.atleast_1d(arr_bg)

    arr_bg_no_0 = arr_bg[np.flatnonzero(arr_bg)]  # Remove zeros to compute mean of background
    BG = arr_bg_no_0.mean() if arr_bg_no_0.size > 0 else 0

    return np.abs(arr - BG) if absol else arr - BG


def compare_distributions(random_samples,
                          region_sample,
                          name="",
                          min_samples=50,
                          qq_threshold=30,
                          perm_n=10000,
                          save_plot=False):
    """
    Compare two datasets with adaptive testing:
    1. KS test
    2. Shapiro-Wilk normality on each sample
    3. If region normal -> parametric: Q-Q & ECDF
       If random normal but region not -> permutation test
       If both non-normal -> Mann-Whitney U test
    4. Anderson-Darling test

    Parameters:
    - random_samples: list or 2D array of random sky values
    - region_sample: list or 1D array for specific region
    - name: label for printing and plots
    - min_samples: threshold for warning about KS power
    - qq_threshold: sample size below which Q-Q plotted
    - perm_n: number of permutations for permutation test
    - save_plot: whether to save Q-Q/ECDF plots
    """
    
    if save_plot: plt.close() #stop any empty figures
    # Flatten inputs
    if isinstance(random_samples, (list, np.ndarray)) and hasattr(random_samples, '__len__'):
        arr = np.asarray(random_samples)
        if arr.ndim > 1:
            rand = np.concatenate(arr)
        else:
            rand = arr.flatten()
    region = np.asarray(region_sample).flatten()

    n_rand, n_reg = len(rand), len(region)
    total_n = n_rand + n_reg
    print(f"\n[{name} patch size {number_of_patches}] Sample sizes: random={n_rand}, region={n_reg} (total={total_n})")
    if total_n < min_samples:
        print(f"Warning: total sample size {total_n} < {min_samples}; KS may be underpowered.")

    # 1. KS Test
    ks_stat, ks_p = stats.ks_2samp(rand, region, nan_policy='omit')
    print(f"\n▶ KS Test: statistic={ks_stat}, p-value={ks_p}")

    # 2. Shapiro-Wilk Normality
    print("\n▶ Shapiro-Wilk Test for Normality:")
    print("Parametric tests are suitable if Normality exists")
    def run_shapiro(data, label):
        if 3 <= len(data) <= 5000:
            W, p = stats.shapiro(data)
            print(f"  {label}: {W}, p={p} -> {'normal' if p>0.05 else 'not normal'}")
            return p > 0.05
        else:
            print(f"  {label}: size {len(data)} out of range for Shapiro")
            return False
    rand_norm = run_shapiro(rand, "Random Samples")
    reg_norm  = run_shapiro(region, "Region Sample")

    # 3. Conditional Testing
    if reg_norm:
        print("\nRegion is normal: proceeding with parametric visuals (Q-Q & ECDF).")
        # Q-Q plots
        if n_rand < qq_threshold or n_reg < qq_threshold:
            for data, label in ((rand, 'Random'), (region, 'Region')):
                stats.probplot(data, dist='norm', plot=plt)
                plt.title(f"Q-Q Plot: {name} ({label})")
                if save_plot:
                    plt.savefig(f"Results/qq_{name}_{label}.png", dpi=300)
                else:
                    plt.show()
        # ECDF plot
        from statsmodels.distributions.empirical_distribution import ECDF
        ecdf1, ecdf2 = ECDF(rand), ECDF(region)
        plt.plot(ecdf1.x, ecdf1.y, label='Random')
        plt.plot(ecdf2.x, ecdf2.y, label='Region')
        plt.title(f"ECDF: {name}")
        plt.legend(); plt.xlabel('Value'); plt.ylabel('ECDF')
        if save_plot:
            plt.savefig(f"Results/ecdf_{name}.png", dpi=300)
        else:
            plt.show()

    elif rand_norm and not reg_norm:
        print("\nRandom normal but region not")
        print("Q-Q plots and ECDF plots have been skipped.")
        print("Performing permutation test that resamples the data to assess significance of the observed difference between the two samples. \n(Doesn't rely on strict distributional assumptions)")
        # permutation test on difference in means
        obs_diff = np.abs(np.mean(rand) - np.mean(region))
        pool = np.concatenate([rand, region])
        count = 0
        for _ in range(perm_n):
            perm = np.random.permutation(pool)
            if np.abs(np.mean(perm[:n_rand]) - np.mean(perm[n_rand:])) >= obs_diff:
                count += 1
        p_perm = count / perm_n
        print(f"Permutation test p-value={p_perm}")

    else:
        print("\nBoth non-normal: using Mann-Whitney U (non-parametric).")
        mw_stat, mw_p = stats.mannwhitneyu(rand, region, alternative='two-sided')
        print(f"Mann-Whitney U: statistic={mw_stat}, p-value={mw_p}")

    # 4. Anderson-Darling
    print("\n▶ Anderson-Darling: \n To determine how closely region sample is drawn from random sample")
    ad_res = stats.anderson_ksamp([rand, region],
                                  method=stats.PermutationMethod() #Using full permutation distribution instead of default approximation to avoid flooring the p-value
                                  )
    ad_stat, ad_sig = ad_res.statistic, ad_res.significance_level/100
    print(f"statistic={ad_stat}, p-value={ad_sig}")
    print("="*50)

"""
Note that it might take too long to fit patches that dont overlap each other 
if the number of patches are too many and/or the size of the patches are too big
"""
patch_size = 30 #in degrees (same as M31 Virial Radius)

"""IMPORTANT"""
number_of_patches = args.patch_num #Creating laaaarge nubmer of patches (choose smaller vlue if you only want to see output features)
    
BINS = bin_num_from_main

if args.pickle: #Then overwrite then process and overwrite existing pickled data
    rm_s, rm_errs = get_real_rm_data()
    
    #Creating random center points for Random R_vir assessment
    patch_ra_points, patch_dec_points = get_random_points(num=number_of_patches, 
                                                    ra_range=(0,360),
                                                    #(min(eq_pos.ra.deg), 
                                                        #           max(eq_pos.ra.deg)),
                                                    dec_range=(-90#(-40+patch_size #Ensuring dec values are valid for arcsin and the ENTIRE patch is above -40 dec
                                                               ,90),
                                                    # (min(eq_pos.dec.deg),
                                                                # max(eq_pos.dec.deg)),
                                                    overlap=True, #Allowing circles to overlap
                                                    radius = patch_size, #in degrees
                                                    min_distance_deg = patch_size #minimum separation of points in degrees
                                                    )

    Patch_pos = SkyCoord(ra=patch_ra_points, dec=patch_dec_points, unit=(u.deg,u.deg), frame='icrs')
    print("Random points converted to SkyCoord object")

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
    
    #____________________________________________________________
    #Capturing Data per patch to udergo proper background subtraction
    #____________________________________________________________
    # CGM_RM_coords_per_patch = [#Gathering coordinates of RM landing in Random R_vir
    #      coords[mask]
    #      for coords, mask in zip(RM_coords_per_patch, [sep < L_m31 * u.deg for sep in RM_coords_sep])
    #  ]
    # BG_RM_coords_per_patch = [#Gathering coordinates of RM landing in Random R_vir's BG region
    #      coords[mask]
    #      for coords, mask in zip(RM_coords_per_patch, [sep > L_m31 * u.deg for sep in RM_coords_sep])
    #  ]
    # CGM_RM_values_per_patch = [#Gathering values of RM landing in Random R_vir
    #      np.asarray(values)[mask]
    #      for values, mask in zip(RM_values_per_patch, [sep.value < L_m31 for sep in RM_coords_sep])
    #  ]
    # BG_RM_values_per_patch = [#Gathering values of RM landing in Random R_vir's BG region
    #      np.asarray(values)[mask]
    #      for values, mask in zip(RM_values_per_patch, [sep.value > L_m31 for sep in RM_coords_sep])
    #  ]

    def filter_RM_coords(coords, sep, threshold, condition):
        """Helper function to filter RM coordinates based on a condition."""
        mask = condition(sep, threshold)
        return coords[mask]

    def filter_RM_values(values, sep, threshold, condition):
        """Helper function to filter RM values based on a condition."""
        mask = condition(sep.value, threshold)
        return np.asarray(values)[mask]

    def process_patch(coords_list, values_list, sep_list, L_m31):
        """Processes RM coordinates and values in parallel."""
        with ThreadPoolExecutor() as executor:
            CGM_RM_coords_per_patch = list(executor.map(
                filter_RM_coords, coords_list, sep_list, [L_m31 * u.deg] * len(sep_list), [np.less] * len(sep_list)
            ))
            BG_RM_coords_per_patch = list(executor.map(
                filter_RM_coords, coords_list, sep_list, [L_m31 * u.deg] * len(sep_list), [np.greater] * len(sep_list)
            ))
            CGM_RM_values_per_patch = list(executor.map(
                filter_RM_values, values_list, sep_list, [L_m31] * len(sep_list), [np.less] * len(sep_list)
            ))
            BG_RM_values_per_patch = list(executor.map(
                filter_RM_values, values_list, sep_list, [L_m31] * len(sep_list), [np.greater] * len(sep_list)
            ))
            CGM_RM_coords_sep = list(executor.map(#Filtering out separation lists as well (Smaller than L_m31)
                filter_RM_coords, sep_list, sep_list, [L_m31 * u.deg] * len(sep_list), [np.less] * len(sep_list)
            ))
            BG_RM_coords_sep = list(executor.map(#Filtering out separation lists as well (Greater than L_m31)
                filter_RM_coords, sep_list, sep_list, [L_m31 * u.deg] * len(sep_list), [np.greater] * len(sep_list)
            ))

        return (CGM_RM_coords_per_patch, BG_RM_coords_per_patch, 
                CGM_RM_values_per_patch, BG_RM_values_per_patch,
                CGM_RM_coords_sep, BG_RM_coords_sep)

# if args.pickle:
#     pickle_filename = "../Patch_pos.pkl"
#     data_to_pickle = {
#             "Patch_pos": Patch_pos
#         }
#     with open(pickle_filename, "wb") as f:
#         pickle.dump(data_to_pickle, f)
#     print("Random Patch Positions have been pickled individually") ; import sys; sys.exit()

if args.pickle:
    #Running in parallel
    (CGM_RM_coords_per_patch, BG_RM_coords_per_patch, 
    CGM_RM_values_per_patch, BG_RM_values_per_patch,
    CGM_RM_coords_sep, BG_RM_coords_sep) = process_patch(
        RM_coords_per_patch, RM_values_per_patch, RM_coords_sep, L_m31
    )
    print("Data for Random R_vir's have been structured for Proper BG substraction")
    pickle_filename = "../Structured_data_for_BG_Subtraction.pkl"
    data_to_pickle = {
            "RM_coords_sep": RM_coords_sep,
            "CGM_RM_coords_sep": CGM_RM_coords_sep,
            "BG_RM_coords_sep": BG_RM_coords_sep,
            "CGM_RM_coords_per_patch": CGM_RM_coords_per_patch,
            "CGM_RM_values_per_patch": CGM_RM_values_per_patch,
            "BG_RM_coords_per_patch": BG_RM_coords_per_patch,
            "BG_RM_values_per_patch": BG_RM_values_per_patch
        }
    with open(pickle_filename, "wb") as f:
        pickle.dump(data_to_pickle, f)
    print("Structured Data for Proper background subtraction has been saved successfully")

#Loading structured data for background subtraction
RM_coords_sep, CGM_RM_coords_sep, BG_RM_coords_sep, CGM_RM_coords_per_patch, CGM_RM_values_per_patch, BG_RM_coords_per_patch, BG_RM_values_per_patch = (
    load_pickle("../Structured_data_for_BG_Subtraction.pkl", 
                "Structured Data for Proper background subtraction has been reloaded from pickle file successfully")
)

if args.pickle:
    #____________________________________________________________
    # Undergoing Proper BG Correction for each Random Patch (using multithreading)
    #____________________________________________________________
    # CGM_RM_values_per_patch_corr = [BG_correction(rm_coords, rm_values, bg_coords, bg_values)
    #                            for rm_coords, rm_values,bg_coords, bg_values in
    #                            zip(CGM_RM_coords_per_patch, CGM_RM_values_per_patch,
    #                                BG_RM_coords_per_patch, BG_RM_values_per_patch)]
    with ThreadPoolExecutor() as executor:
        CGM_RM_values_per_patch_corr = list(executor.map(
            BG_correction, CGM_RM_coords_per_patch, CGM_RM_values_per_patch, BG_RM_coords_per_patch, BG_RM_values_per_patch
        ))
    print("Rotation Measures have been Corrected !! :)")

    pickle_filename = "../RM_Subtracted_with_spline.pkl"
    data_to_pickle = {
            "CGM_RM_values_per_patch_corr": CGM_RM_values_per_patch_corr,
        }
    with open(pickle_filename, "wb") as f:
        pickle.dump(data_to_pickle, f)
    print(" Properly Subtracted RM's have been dumped in pickle file successfully.")

#Loading properly subtracted RM values
(CGM_RM_values_per_patch_corr,) = (
    load_pickle("../RM_Subtracted_with_spline.pkl", 
                "Properly Subtracted RM's have been reloaded from pickle file successfully")
)

if args.pickle:
    all_d_bin_centers=[] #For x-axis
    all_means = []
    all_medians = []
    all_bin_stds = []
    all_bin_edges = ''

    fig2 = plt.figure(figsize=(10, 5))
    ax2 = fig2.add_subplot(111)
    print("Mean and Median calculations have begun")
    for i in range(len(RM_coords_sep)): #Searching through each patch
        
        if not (CGM_RM_coords_per_patch[i].shape == ()):  #Checking if its not empty or filled with NaNs
            if not (CGM_RM_coords_per_patch[i].shape[0] < 30):  # Ensure sufficient number of points for stats module to work
                d_bin_centers, bin_mean, bin_med, bin_std, bin_edges = get_mean_and_med_stats(CGM_RM_coords_sep[i], CGM_RM_values_per_patch_corr[i], bin_num=BINS)
                
                all_d_bin_centers.append(d_bin_centers) #For x axis
                all_means.append(bin_mean) #For y-axis
                all_medians.append(bin_med) #For y-axis
                all_bin_stds.append(bin_std)

                if i == 0 : #Only collect bin edges once
                    all_bin_edges = bin_edges

                # #Background correction for the individual plots
                # bin_mean_1 = indiv_bg_corr(bin_mean, d_bin_centers)
                # bin_med_1 = indiv_bg_corr(bin_med, d_bin_centers)

                #This has been commented out to remove clatter
                #The they are all being collected and will be averaged to make a final one
                #plot_indidividual_patch_stats(ax2, d_bin_centers, bin_mean_1, bin_med_1, bin_std)
            
    #Data to pickle
    data_to_pickle = {
        "ax2": ax2,
        "fig2": fig2,
        "RM_coords_sep": RM_coords_sep,
        "all_d_bin_centers": all_d_bin_centers,
        "all_means": all_means,
        "all_medians": all_medians,
        "all_bin_stds": all_bin_stds,
        "all_bin_edges": all_bin_edges,
        "rm_s": rm_s,
        "rm_errs": rm_errs,
        "patch_ra_points": patch_ra_points,
        "patch_dec_points": patch_dec_points,
        "Patch_pos": Patch_pos
    }

    #Save to pickle file
    with open("../RM_stats.pkl", "wb") as f:
        pickle.dump(data_to_pickle, f)
    print("Mean, Median calculations, and additional variables have been pickled successfully!")

#Loading RM statistics
(ax2, fig2, RM_coords_sep, all_d_bin_centers, all_means, all_medians, 
all_bin_stds, all_bin_edges, rm_s, rm_errs, patch_ra_points, 
patch_dec_points, Patch_pos) = load_pickle("../RM_stats.pkl", 
                "Pickled Mean & Median data loaded successfully!")

pickle_filename = os.path.join("..", "saved_data.pkl")

if args.pickle:
    D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                    max([max(centers) for centers in all_d_bin_centers]), 
                                    num=BINS)

    Avg_means = np.mean(all_means, axis=0) #Average Means for spline corrected Random R_vir
    Avg_medians = np.mean(all_medians, axis=0) #Average Medians for spline corrected Random R_vir
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

    # Variables to pickle
    data_to_save = {
        "D_bin_centers": D_bin_centers,
        "Avg_means": Avg_means,
        "Avg_medians": Avg_medians,
        "Avg_means_std": Avg_means_std,
        "Avg_medians_std": Avg_medians_std,
        "int_Avg_means": int_Avg_means,
        "int_Avg_medians": int_Avg_medians,
        "int_Avg_means_std": int_Avg_means_std,
        "int_Avg_medians_std": int_Avg_medians_std,
        "int_D_bin_centers": int_D_bin_centers
    }

    # Save data to pickle file
    with open(pickle_filename, "wb") as f:
        pickle.dump(data_to_save, f)

    print(f"Data for M31_vs_entire sky plotting saved to {pickle_filename}")

# Load the data immediately
with open(pickle_filename, "rb") as f:
    loaded_data = pickle.load(f)

(#Unpack loaded data for continued use
    D_bin_centers,
    Avg_means,
    Avg_medians,
    Avg_means_std,
    Avg_medians_std,
    int_Avg_means,
    int_Avg_medians,
    int_Avg_means_std,
    int_Avg_medians_std,
    int_D_bin_centers
) = (loaded_data[key] for key in loaded_data)
del loaded_data
print("Data reloaded for original M31_vs_entire_sky plot.")

if __name__ == "__main__": #continue (this makes it easier to excecute "M31_signal_density.py" file)

    if args.original_plot:
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
        if not args.save_plot:
            plt.show()
        else:
            path = curr_dir_path() + "Results/"
            fig2.savefig(f"{path}M31_signal_vs_entire_sky_{number_of_patches}_patches.png", dpi=600, bbox_inches="tight") #saving the image
            print(f"M31_signal_vs_entire_sky_{number_of_patches}_patches.png has been successfully saved in Results directory")
    
        plt.close(fig2) #Deleting the Figure
    
    # #Data close to Gundo's shade_ms plots to spot any storng outliers
    # Shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_Avg_medians, int_Avg_medians_std, int_D_bin_centers)
    if args.test_patches: #show patches on sphere as they get smaller
        test_patches_on_sphere()

    if args.CGM_dispersion:
        from main import d_bin_centers
        plot_dispersion(bin_num=bin_num_from_main, 
                        sep_vals=m31_sep_Rvir,
                        rm=rm_m31,
                        centers=d_bin_centers,
                        save_plot=args.save_plot,
                        bin_std_random=Avg_medians_std if args.median else (Avg_means_std if args.mean else None))

    if args.annuli_anal: 
        annuli_analysis_random(all_means, all_medians, save_plot=True)

    elif args.m31_annuli_anal:
        annuli_analysis_m31(save_plot=True)
    
    if args.rm_vs_PA or args.rm_vs_gal_lat:
        Sep_vals = np.concatenate([m31_sep_Rvir, m31_sep_bg])
        RM_and_BG = np.concatenate([rm_m31, rm_bg])
        m31_distances = tuple(map(get_projected_d, m31_sep_Rvir, [d_m31]*len(m31_sep_Rvir))) #from separated distance in degrees - within CGM of M31 - to kpc
        m31_distances = list(map(lambda m31_d: m31_d.value, m31_distances)) #Capturing only their value
    
    if args.rm_vs_PA: 
        PA_rm_deg = (PA_rm.deg -37.7) %360 #Converting from radians to degrees
        PA_rm_deg_bg = (PA_bg.deg -37.7) %360 #Converting from radians to degrees
        PA_rm_degs = np.concatenate([PA_rm_deg, PA_rm_deg_bg]) #Including the Background

        # print(f"{len(m31_distances)=}")
        # print(f"{len(PA_rm_deg)=}")
        # print(f"{len(rm_m31)=}")
        [distance_per_PA_bin, rm_per_PA_bin], _, bin_edges = create_annuli_binning_structure(bin_data=m31_distances, 
                                                                                             data=(PA_rm_deg, rm_m31), 
                                                                                             bin_num=bin_num_from_main+1, 
                                                                                             for_PA_or_b_plot=True)
        plot_binned_PA(distance_per_PA_bin, rm_per_PA_bin, bin_edges, save_plot=True)

    if args.rm_vs_gal_lat: 
        rm_pos_gal_lats = np.concatenate([rm_pos_gal_lat, rm_pos_gal_lat_bg]) #Inclusive of its background (if used...)
        # print(rm_pos_gal_lats); import sys; sys.exit()'
            
        [distance_per_b_bin, rm_per_b_bin], _, bin_edges = create_annuli_binning_structure(bin_data=m31_distances, 
                                                                                             data=(rm_pos_gal_lat, rm_m31), 
                                                                                             bin_num=bin_num_from_main+1, 
                                                                                             for_PA_or_b_plot=True)
        plot_binned_gal_lat(distance_per_b_bin, rm_per_b_bin, bin_edges, save_plot=True)
        
    if args.m31_ks_test:
        """
        Actually does more than a KS Test Depending on the data of the two smamples.
        Tests for Normalcy of the data and a series of conditions are added for further tests if KS cannot be conducted
        """
        if args.mean or not (args.mean or args.median): #Mean or both
            compare_distributions(Avg_means, bin_means_m31, name="Mean RM", save_plot=args.save_plot)
        if args.median or not (args.mean or args.median): #Mean or both
            compare_distributions(Avg_medians, bin_med_m31, name="Median RM", save_plot=args.save_plot)

    if args.print_ks_test:

        headers = ["Test Type", "10 R_vir", "100 R_vir", "1K R_vir", "10K R_vir"]
        index_labels = ["KS Test (stat, p)", "Shapiro - Random (W, p)", "Shapiro - Region (W, p)", "Permutation (p-value)", "Anderson-Darling (stat, p)"]

        # Mean RM test values
        mean_data = [
            ["1.0, 1.08e-05"] * 4,
            ["0.978, 0.953", "0.934, 0.485", "0.918, 0.340", "0.934, 0.485"],
            ["0.827, 0.031"] * 4,
            ["0.0"] * 4,
            ["9.857, 1e-06"] * 4
        ]

        # Median RM test values
        median_data = [
            ["1.0, p=1.08e-05"] * 4,
            ["0.923, 0.383", "0.986, 0.990", "0.988, 0.994", "0.923, 0.383"],
            ["0.868, 0.094"] * 4,
            ["-", "-", "-", "-"],
            ["9.857, 1e-06"] * 4
        ]

        # Create DataFrames
        mean_df = pd.DataFrame(mean_data, columns=headers[1:], index=index_labels)
        mean_df.insert(0, headers[0], index_labels)

        median_df = pd.DataFrame(median_data, columns=headers[1:], index=index_labels)
        median_df.insert(0, headers[0], index_labels)

        # Print tables nicely
        print("\nMean")
        print(tabulate(mean_df, headers='keys', tablefmt='fancy_grid', showindex=False))

        print("\nMedian")
        print(tabulate(median_df, headers='keys', tablefmt='fancy_grid', showindex=False))

    if args.rm_per_patch_hist or args.plot_random_cosdec:

        # Step 1: Convert patch center points to Galactic latitude (B)
        Random_cent_points_b = Patch_pos.galactic.b.deg

        # Step 2: Count RM values per patch
        rm_counts_per_patch = np.array([len(rm) for rm in CGM_RM_values_per_patch])

        if args.b_limit_small:
            filtered_latitudes = Random_cent_points_b
            filtered_rm_counts = rm_counts_per_patch
        else:#If using |b| > 5 even before spline correction
            mask = np.abs(Random_cent_points_b) >= b_upper_lim.value
            filtered_latitudes = Random_cent_points_b[mask]
            filtered_rm_counts = rm_counts_per_patch[mask]
        
        m31_lat_deg = m31_lat.deg
        def generate_lat_labels(rm_binned, cos_lat_bins, m31_lat_deg, m31_mark=True):
            """
            Generate cos(b) bin labels including count info and optional M31 b tag.
            """
            cos_m31 = abs(np.cos(np.deg2rad(m31_lat_deg)))
            labels = []
            for i in range(len(cos_lat_bins) - 1):
                count = len(rm_binned[i])
                label = f"{cos_lat_bins[i]:.2f} < cos(b) < {cos_lat_bins[i + 1]:.2f} [{count}]"
                if m31_mark:
                    if cos_lat_bins[i] <= cos_m31 < cos_lat_bins[i + 1]:
                        label += r"  ($b_{M31} \approx " + f"{m31_lat_deg:.2f}" + r"^\circ$)"
                labels.append(label)
            return labels
        
        # Step 3: Define cos(b) bins (absolute value, so from 0 to 1)
        bins = 10
        cos_lat_bins = np.linspace(0, 1, bins + 1)

        # Step 4: Initialize binned RM container
        rm_binned = [[] for _ in range(len(cos_lat_bins) - 1)]

        # Step 5: Separate positive and negative latitude RM counts into bins
        rm_binned_pos = [[] for _ in range(len(cos_lat_bins) - 1)]
        rm_binned_neg = [[] for _ in range(len(cos_lat_bins) - 1)]

        for i in range(len(filtered_rm_counts)):
            lat = filtered_latitudes[i]
            cos_lat = abs(np.cos(np.deg2rad(lat)))
            for j in range(len(cos_lat_bins) - 1):
                if cos_lat_bins[j] <= cos_lat < cos_lat_bins[j + 1]:
                    if lat >= 0:
                        rm_binned_pos[j].append(filtered_rm_counts[i])
                    else:
                        rm_binned_neg[j].append(filtered_rm_counts[i])
                    break
            
    if args.rm_per_patch_hist:
        """
        Creating Histogram plot of Number Random of Patches that have a certian number of collected RM points
        """

        #Step 6
        lat_labels_neg = generate_lat_labels(rm_binned_neg, cos_lat_bins, m31_lat_deg)
        lat_labels_pos = generate_lat_labels(rm_binned_pos, cos_lat_bins, m31_lat_deg, m31_mark=False)

        # Step 7: Get colors
        colors = get_discrete_colors((0, 1), len(rm_binned_pos), cmap_name="tab20").colors

        # Step 8: Plot in two subplots
        fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

        axs[0].hist(
            rm_binned_neg,
            bins=bins,
            stacked=True,
            color=colors,
            edgecolor='black',
            label=lat_labels_neg,
            linewidth=0.3
        )

        axs[0].set_title("Negative Galactic Latitude")
        axs[0].set_xlabel("RM Count")
        axs[0].set_ylabel("Frequency")
        axs[0].grid(True, axis="y", which="both")

        axs[1].hist(
            rm_binned_pos,
            bins=bins,
            stacked=True,
            color=colors,
            edgecolor='black',
            label=lat_labels_pos,
            linewidth=0.3
        )

        axs[1].set_title("Positive Galactic Latitude")
        axs[1].set_xlabel("RM Count")
        axs[1].grid(True, axis="y", which="both")

        for ax in axs:
            ax.legend(
                fontsize='small',
                title="Galactic latitude " + r"[$^{\circ}$] " + " (|b|>" + f"{b_upper_lim.value:.1f}" + r"$^{\circ}$)",
            )
        plt.tight_layout(rect=[0, 0.05, 1, 1])

        # Save
        path = curr_dir_path() + "Results/"
        plt.savefig(f"{path}Stacked_Histogram_Frequency_of_RMcounts.png", dpi=600, bbox_inches="tight")
        print(f"Stacked Histogram plot saved to {path}")

    if args.plot_random_cosdec:
        def plot_random_points_sin_dec(Patch_pos):

            # Step 1: Get Dec values in degrees from the random patch positions
            dec_vals = Patch_pos.dec.rad
            dec_vals = Patch_pos.transform_to('galactic').b.rad

            # Step 2: Calculate sin(dec) values
            sin_dec_vals = np.cos(dec_vals)

            # Step 3: Define bins for sin(dec), from -1 to 1
            bins = np.linspace(-1, 1, 20 + 1)  # 0.1 bins

            # Step 4: Separate sin(dec) into positive and negative parts
            sin_dec_neg = sin_dec_vals[dec_vals < 0]
            sin_dec_pos = sin_dec_vals[dec_vals >= 0]

            # Step 5: Plot
            fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

            axs[0].hist(sin_dec_neg, bins=bins, color='skyblue', edgecolor='black')
            axs[0].set_title("Negative Galactic Latitude")
            axs[0].set_xlabel(r"$\sin(b)$")
            axs[0].set_ylabel("Number of Random Points")
            axs[0].set_xlim(-1,0)
            axs[0].grid(True)

            axs[1].hist(sin_dec_pos, bins=bins, color='salmon', edgecolor='black')
            axs[1].set_title("Positive Galactic Latitude")
            axs[1].set_xlabel(r"$\sin(b)$")
            axs[0].set_xlim(-1,1)
            axs[1].grid(True)
            plt.show()

        plt.close()
        plot_random_points_sin_dec(Patch_pos)

    if args.cdf_anal: #Cumulative Density Function (Different from CDF of RM Denisty plots)
        
        def plot_cdf(rm_m31, rm_sky, label_suffix="Means"):

            plt.close() #Incase a random figure appears.
            """
            Plots the cumulative distribution of the absolute RM values for M31 and the sky.
            
            Parameters:
            - rm_m31: array-like, RM values for M31 region
            - rm_sky: array-like, RM values for the rest of the sky
            - label_suffix: str, optional label suffix for distinguishing plot legends
            """

            # Take absolute values
            rm_m31_abs = np.abs(rm_m31)
            rm_sky_abs = np.abs(rm_sky)

            # Sort values
            rm_m31_abs_sorted = np.sort(rm_m31_abs)
            rm_sky_abs_sorted = np.sort(rm_sky_abs)

            # Cumulative sum normalized
            cumulative_m31 = np.cumsum(rm_m31_abs_sorted) / np.sum(rm_m31_abs_sorted)
            cumulative_sky = np.cumsum(rm_sky_abs_sorted) / np.sum(rm_sky_abs_sorted)

            plt.figure(figsize=(10, 6))
            ls = 15
            plt.plot(rm_m31_abs_sorted, cumulative_m31, 
                    label=f"M31 ({label_suffix})", color="blue")
            plt.plot(rm_sky_abs_sorted, cumulative_sky, 
                    label=f"Rest of the Sky - {label_suffix}", color="red", linestyle="--")

            plt.xlabel("|RM| [rad m$^{-2}$]", fontsize=ls)
            plt.ylabel("Cumulative Fraction", fontsize=ls)
            plt.legend(fontsize=12)
            plt.xlim(0,)
            plt.ylim(0,)
            plt.grid(True, alpha=0.5)
            plt.show()

        plot_cdf(bin_means_m31, Avg_means, label_suffix="Means")
        plot_cdf(bin_med_m31, Avg_medians, label_suffix="Medians")

        