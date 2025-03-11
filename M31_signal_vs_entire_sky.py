from scipy.interpolate import interp1d
import seaborn as sns

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

#importing functions
get_projected_d_old, confining_circle
)

#Hanlding unnneccesary clutter of printing form warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter("ignore", RuntimeWarning) #Suppresses stats SmallSampleWarning when using stats.sem()

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
        statistic=lambda rm_vals: stats.sem(rm_vals),  # Standard Error of Mean
        bins=bin_num)

    return d_bin_centers, bin_means, bin_med, bin_std

def plot_indidividual_patch_stats(ax, d_bin_centers, bin_mean, bin_med, bin_std):
    # Plot the mean and median with error bars
    ax.errorbar(d_bin_centers, np.absolute(bin_mean), yerr=bin_std, fmt='k.-', alpha=.4)
    ax.errorbar(d_bin_centers, np.absolute(bin_med), yerr=bin_std, fmt='g.-', alpha=.4, capsize=2)

def plot_m31_stats(ax):
    ax.errorbar(d_bin_centers_m31, np.absolute(bin_means_m31), yerr = bin_std_m31, 
                color="black", fmt = '.-', alpha=.6)#, label ="$\mu_{M31}$")
    ax.errorbar(d_bin_centers_m31, np.absolute(bin_med_m31), yerr=bin_std_m31, 
                color='orange', fmt='.-', capsize=2, markeredgecolor="k", alpha=.6)

def plot_m31_dispersion(bin_num):

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
    plt.show()
    
def Shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_Avg_medians, int_Avg_medians_std, int_D_bin_centers):
    y_upper=100
    def shade_ms_mimic(data,#average mean/ medians (interpolated)
                        std,#standard deviation (interpolated)
                        bin_centers, #Bin centers
                        cmap, name):
            
        X, Y = np.meshgrid(bin_centers, np.linspace(np.min(data-std), np.max(data+std), 100))

        #Computing dispersion-density based on the fill_between regions
        Z = np.exp(-((Y - data[:, None])**2 / (2 * std[:, None]**2)))

        fig, ax = plt.subplots(figsize=(10, 5))
        sns.heatmap(Z, cmap=cmap, alpha=0.6, cbar=True, xticklabels=False, yticklabels=False)
        y = np.arange(0,y_upper+10, 10)
        ax.set_yticks(y); ax.set_yticklabels(y)
        ax.set_xlabel(r'R$_{projected}$ [kpc]', fontsize=12)
        ax.set_ylabel('|RM|', fontsize=12)
        ax.set_ylim(0, y_upper) #np.max(data + std))
        plt.title('Dispersion of RM Variations ('+name+')')

    shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_D_bin_centers, name="Means", cmap="viridis")
    shade_ms_mimic(int_Avg_medians, int_Avg_medians_std, int_D_bin_centers, cmap="viridis", name="Medians")
    
    plt.show()

#Note that it might take too long to fit patches that dont overlap each other 
#If the number of patches are too many and/or the size of the patches are too big
patch_size = 30 #in degrees (same as M31 Virial Radius)

"""IMPORTANT"""
number_of_patches = int(8e3) #Creating laaaarge nubmer of patches (choose smaller vlue if you only want to see output features)


BINS = 30
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

# RM_coords_per_patch will be filled with values when collection_of_points_from_WCS_sphere() is called
# Same goes for RM_values_per_patch
RM_values_per_patch, RM_coords_per_patch = [], []

#Creating axis that will be used when doing WCS axis point collections
fig1 = plt.figure(figsize=(6, 5))
ax1 = fig1.add_subplot(111, projection=get_wcs("LGSNLGSR.SQLGBB.FITS"), slices=('x', 'y', 0, 0))
fig1.clf() #No need to show figure. Only axis is needed

print("(No longer plotting) but collecting of points from each circular patch has begun...")
collection_of_points_from_WCS_sphere() #IMPORTANT
print("Collection of points by each circular patch on WCS sphere is complete")

fig2 = plt.figure(figsize=(10, 5))
ax2 = fig2.add_subplot(111)

print("Getting separation of RM from center of relative patch")
#Get separation of RM from center of relative patch.
RM_coords_sep = [rm_coords.separation(patch_pos) 
                 for rm_coords, patch_pos in 
                 list(zip(RM_coords_per_patch, Patch_pos))]

Max_med, Min_med, Max_mean, Min_mean, D_bin_centers = [None]*5

all_d_bin_centers=[] #For x-axis
all_means = []
all_medians = []
for i in range(len(RM_coords_sep)):
    
    if i == 0: print("Mean and Median calculations have begun")
    
    if not (RM_coords_per_patch[i].shape == ()):  #Checking if its not empty or filled with NaNs
        if not (RM_coords_per_patch[i].shape[0] < 15):  # Ensure sufficient number of points for stats module to work
            d_bin_centers, bin_mean, bin_med, bin_std = get_mean_and_med_stats(RM_coords_sep[i], RM_values_per_patch[i], bin_num=BINS)
            
            all_d_bin_centers.append(d_bin_centers) #For x axis
            all_means.append(bin_mean) #For y-axis
            all_medians.append(bin_med) #For y-axis
            
            #This has been commented out to remove clatter
            #The they are all being collected and will be averaged to make a final one
            # plot_indidividual_patch_stats(ax2, d_bin_centers, bin_mean, bin_med, bin_std)

#getting mean of background
D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                   max([max(centers) for centers in all_d_bin_centers]), 
                                   num=len(all_means[0]))


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
global D_bin_centers
D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                   max([max(centers) for centers in all_d_bin_centers]), 
                                   num=BINS #BINS #Number of points for interpolation (smoothens it out)
                                         #Must be same as bin_num parameter in function "get_mean_and_med_stats"
                                   )

# Find a common range of D_bin_centers for interpolation
int_D_bin_centers = np.linspace(min([min(centers) for centers in all_d_bin_centers]), 
                                   max([max(centers) for centers in all_d_bin_centers]), 
                                   num=100 #Number of points for interpolation (smoothens it out)
                                         #Msut be same as interploted num_points
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
ax2.set_xlim(0, 300)
ax2.set_ylim(0,)

# ax2.fill_between(D_bin_centers, Avg_means_std , Avg_means_std, color='blue', alpha=0.2, label='Mean standard deviation of the patches')
# ax2.fill_between(D_bin_centers, Avg_medians_std, Avg_medians_std, color='green', alpha=0.4, label='Median standard deviation of the patches')

ax2.legend(fontsize = 12, loc = 'upper center', bbox_to_anchor = (.5, 1.2), 
            framealpha = 0, ncols = (2,4))
plt.tight_layout()
plt.show()

#Data close to Gundo's shade_ms plots to spot any storng outliers
Shade_ms_mimic(int_Avg_means, int_Avg_means_std, int_Avg_medians, int_Avg_medians_std, int_D_bin_centers)

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
        
if (inpt := input("Want to show patches on spehre as they get smaller [Y,N]?").lower()) in ['y','yes']:
    test_patches_on_sphere()

if (inpt := input("Display Dispersion vs RM [Y,N]?").lower()) in ['y','yes']:
    plot_m31_dispersion(bin_num_from_main)
