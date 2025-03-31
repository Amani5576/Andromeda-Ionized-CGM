#Note that this script only declares variables and functions. No ouput graps are produced
from astropy.table import Table
from astropy import units as u
import numpy as np
import csv
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from scipy import stats
from astropy.coordinates import SkyCoord
from matplotlib.patches import Rectangle, Circle, Path
from scipy.optimize import curve_fit
from astropy.utils.exceptions import AstropyWarning
import warnings
from scipy.interpolate import RectBivariateSpline, griddata

# Suppress specific Astropy warnings
warnings.simplefilter('ignore', AstropyWarning)

#to be used for scaling size of scatter based on their rm value
sizes = lambda vals: 30 #np.abs(vals) * 2.3

R_vir = 300.*u.kpc #Virial Radius of Andromeda in kiloparsecs (kpc)
cutoff = 30.*u.deg #(in deg) Limit for taking to account background correction
d_m31 = 780.*u.kpc #Distance to Andromeda in kiloparsecs (kpc)

#Read in catalog.dat file
t = Table.read("catalog.dat", format = "ascii")

t.rename_column("col9", "l")
t.rename_column("col10", "b")
t.rename_column("col17", "RM")
t.rename_column("col18", "e_RM")

t.keep_columns(["l", "b", "RM", "e_RM"])

RM_unit = u.rad/u.m**2

"""IMPORTANT"""
sigma_detect_limit = 0

sig_mask = ~(sigma_detect_limit*np.abs(t["e_RM"])>np.abs(t["RM"]))
t["l"].unit = u.deg
t["b"].unit = u.deg
t["RM"].unit = RM_unit
t["e_RM"].unit = RM_unit
rm_raw = t["RM"]
rm_err_raw = t["e_RM"]
RM_lat = t["l"][sig_mask]
RM_lon = t["b"][sig_mask]
rm = rm_raw[sig_mask] #rm must be larger than e_RM
rm_err = rm_err_raw[sig_mask]

m31_maj = 73.*u.arcmin #Major axis diameter from POSS in arcmin
m31_min = 45.*u.arcmin #Minor axis diameter from POSS in arcmin
m31_pa = 37.7*u.deg #PA CCW from North from de Vaucouleurs et al. 1958.  

position = SkyCoord(RM_lat, RM_lon, unit = "degree", frame = "galactic")

eq_pos = position.transform_to('icrs') #Equitorial positions

m31_pos = SkyCoord(ra = "00:42:44.35", dec = "+41:16:08.6", unit = (u.hourangle, u.deg), frame = 'icrs')

new_frame_of_reference = m31_pos.skyoffset_frame()

#Convert coordinates to M31 reference frame
rm_m31_coord = eq_pos.transform_to(new_frame_of_reference)

# M31 separations
m31_sep = eq_pos.separation(m31_pos)
m31_theta = eq_pos.position_angle(m31_pos) #Polar angle

#Distance of Cloud 6 from M31 (just because)
cloud6_pos = SkyCoord("01:08:29.6 +37:45:00", unit = (u.hourangle, u.deg), frame = 'icrs')
#print(cloud6_pos)

m33_pos = SkyCoord("23.462042 30.660222", unit = (u.deg, u.deg), frame = 'icrs')

m33_m31coord = m33_pos.transform_to(new_frame_of_reference)

m33_sep = m33_m31coord.separation(m31_pos) #Angular distance from m31.
m33_theta = m33_m31coord.position_angle(m31_pos)

          #Virial radius,#distance to Andromeda
L_m31 = (np.arctan(R_vir/d_m31)*u.rad.to(u.deg)).value #no longer uses small angle approximation.

#Creating conditions.         EnsuresRM measurement for bg is within 30 degrees of M31
bg_condition = (m31_sep.deg > L_m31) & (m31_sep.deg < cutoff.value)
m31_condition = m31_sep.deg <= L_m31

#Applying conditions to filter RM positions for CGM of M31 and its BG alone.
bg_pos = rm_m31_coord[bg_condition]
bg_pos_icrs = bg_pos.transform_to("icrs") 
rm_pos = rm_m31_coord[m31_condition]
rm_pos_icrs = rm_pos.transform_to("icrs") 

#Done in year 2025 for plotting RM vs Galactic Azimuth
rm_pos_gal_lat = position.b.deg[m31_condition]
rm_pos_gal_lat_bg = position.b.deg[bg_condition]


# #Applying conditions to filter RM values and their errors
rm_bg = rm[bg_condition] #Record background RM
m31_sep_bg = (m31_sep.deg[bg_condition])*u.deg #Record background RM positions (angular)
err_bg = rm_err[bg_condition] #Record background RM error

rm_m31 = rm[m31_condition] #Record m31's RM
m31_sep_Rvir = (m31_sep.deg[m31_condition])*u.deg #Record m31's RM separated values.
err_m31 = rm_err[m31_condition] #Record m31's RM errors

def BG_correction(rm_coords, rm_values, bg_coords, bg_values):
    """
    Perform background subtraction on Rotation Measure (RM) values using 
    RectBivariateSpline interpolation.

    Parameters:
    -----------
    rm_coords : SkyCoord
        Sky coordinates of RM measurements to be background corrected.
    rm_values : np.ndarray
        RM values corresponding to `rm_coords`.
    bg_coords : SkyCoord
        Sky coordinates of background RM measurements.
    bg_values : np.ndarray
        RM values at background positions.

    Returns:
    --------
    np.ndarray
        Background-corrected RM values.
    """
    #Converting to degrees
    x_bg = bg_coords.ra.deg  #(N,)
    y_bg = bg_coords.dec.deg  #(N,)
    
    rm_x = rm_coords.ra.deg  #(M,)
    rm_y = rm_coords.dec.deg  #(M,)
    
    # Define a regular grid for interpolation
    x_grid = np.linspace(x_bg.min(), x_bg.max(), len(x_bg)) #(N,)
    y_grid = np.linspace(y_bg.min(), y_bg.max(), len(x_bg)) #(N,)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid) #each having dimensions (N,N)
    
    #Ensuring griddata inputs have correct dimensions
    bg_points = np.column_stack((x_bg, y_bg))  # (N,2) format required by griddata
    grid = np.column_stack((X_grid.ravel(), Y_grid.ravel()))  # (N*N,2)
    
    #Interpolating background values onto grid
    bg_grid = griddata(bg_points, bg_values, grid, method='linear')  # (N*N,)
    bg_grid = bg_grid.reshape(X_grid.shape)  # Reshape back to (N,N)
    
    # Handle NaNs (replace with nearest-neighbor interpolation)
    if np.isnan(bg_grid).any():
        bg_grid = griddata(points, bg_values, grid_points, method='nearest')
        bg_grid = bg_grid.reshape(X_grid.shape)

    #Fitting spline to interpolated background RM data
    fbeam = RectBivariateSpline(y_grid, x_grid, bg_grid)
    
    #Interpolating the background RM values at RM positions
    bg_values_interp = fbeam.ev(rm_y, rm_x)
    
    #Finally subtracting complex background (interpolated) from given rm coords.
    rm_corrected = rm_values - bg_values_interp
    
    return rm_corrected

#Conducting backgroudn subtraction.  
rm_m31 = BG_correction(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg)
rm_bg = BG_correction(bg_pos_icrs, rm_bg, bg_pos_icrs, rm_bg) #Doing same for background itself
#Validation/Sanity Check
# print(len(rm_m31), len(err_m31), len(np.power(err_m31, 2)), len(rm_bg))

bin_num = 22
#Calculate mean of RM values within Rvir of M31
bin_means, bin_edges, binnumber = stats.binned_statistic(m31_sep_Rvir, rm_m31, statistic = 'mean', bins = bin_num)

#Calculate mean of RM values within Rvir of M31 (inclusive of background RM for R_vir)
(bin_means_past_rvir, 
bin_edges_past_rvir, 
binnumber_past_rvir )= stats.binned_statistic(
    np.concatenate([m31_sep_Rvir,m31_sep_bg]), 
    np.concatenate([rm_m31,rm_bg]), 
    statistic = 'mean', bins = bin_num)

#Calculate mean of RM values within Rvir of M31 (inclusive of background RM for R_vir)
(bin_meds_past_rvir, 
bin_edges_past_rvir, 
binnumber_past_rvir )= stats.binned_statistic(
    np.concatenate([m31_sep_Rvir,m31_sep_bg]), 
    np.concatenate([rm_m31,rm_bg]), 
    statistic = 'median', bins = bin_num)

bin_med, bin_edges, binnumber = stats.binned_statistic(m31_sep_Rvir, rm_m31, statistic = 'median', bins = bin_num)

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2

#Convert angles to linear projected distances (No longer takes to account small angle approximation)
get_projected_d = lambda ang, d: (d * np.tan(ang.to(u.rad))).to(u.kpc)
get_projected_d_old = lambda pos: d_m31*np.arctan(np.radians(pos)) #specifically used for M31 (not RM)

get_sep_angle = lambda d_proj, d: (np.arctan(d_proj/d.value)*u.rad).to(u.deg)

#Convert projected radial distance (in kpc) and known distance to object (in kpc)
#to angular separation (in degrees)
get_angle_sep_from_distance = lambda R_proj, d: np.degrees(np.arctan(R_proj.to(u.kpc)/d.to(u.kpc)))
get_angle_sep_from_distance_reverse = lambda Ang, d: np.tan(Ang.to(u.rad)) * d.to(u.kpc)

d_bg = get_projected_d(m31_sep_bg, d_m31)
d_rm = get_projected_d(m31_sep_Rvir, d_m31) #Is only for within R_vir
d_m33 = get_projected_d(m33_sep, d_m31)
d_bin_centers = get_projected_d(bin_centers*u.deg, d_m31)

# SEMed = lambda sig, N: 1.253*sig/N #Using a bootstrap method of calculating standard error of median

#Using standard error of mean for error bars: 
bin_std, bin_edges, binnumber = stats.binned_statistic(m31_sep_Rvir, 
    rm_m31, 
    statistic = lambda rm_m31: stats.sem(rm_m31), #Standard Error ofMean
    bins = bin_num)  

#Using standard error of mean for error bars:  (inclusive of background RM for R_vir)
(bin_std_past_rvir, 
bin_edges_past_rvir, 
binnumber_past_rvir )= stats.binned_statistic(
    np.concatenate([m31_sep_Rvir,m31_sep_bg]), 
    np.concatenate([rm_m31,rm_bg]), 
    statistic = lambda rm_m31: stats.sem(rm_m31), #Standard Error ofMean
    bins = bin_num) 
    
#Calculating polar angles
shift = 180*u.deg
PA_bg = m31_theta[bg_condition] + shift
PA_rm = m31_theta[m31_condition] + shift
PA_m33 = m33_theta + shift

# =============================================================================
# M33 Radial distance from M31 --> d_m33
# M33 Angular distance from M31 --> m33_sep 
# Background RM Radial distance from M31 --> d_bg
# Background RM Angular distance from M31 --> m31_sep_bg
# R_vir RM Radial distance from M31 --> d_rm 
# R_vir RM Angular distance from M31 --> m31_sep_Rvir
# =============================================================================

def get_wcs(filename):
    hdu = fits.open(filename)[0]
    return WCS(hdu.header)

def update_projection(ax, row_num, col_num, projection, fig=None):
    """Perosnally modified but with assistance from:
    https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
    """
    """Updateprojection of a specific subplot in a grid of subplots.
    
    Parameters:
    - ax: Array of subplots.
    - row_num: Row number (0-indexed) ofsubplot to update.
    - col_num: Column number (0-indexed) ofsubplot to update.
    - projection: New projection type.
    - fig: Figure object (optional).
    
    Returns:
    - Updated axis object.
    """
    if fig is None:
        fig = plt.gcf()  #Get current figure if fig is not provided
        
    try: 
        #Calculateflat index ofsubplot based on row and column numbers
        flat_index = row_num * ax.shape[1] + col_num
    except IndexError:
        flat_index = row_num + col_num
        
    #Get subplot geometry information
    rows, cols, start, stop = ax.flat[flat_index].get_subplotspec().get_geometry()
    
    #Removecurrent subplot
    ax.flat[flat_index].remove()
    
    #Add a new subplot withspecified projection
    ax.flat[flat_index] = fig.add_subplot(rows, cols, flat_index + 1, projection=projection)
    
    return ax.flat[flat_index]

#Function for plotting an ellipse on a polar plot based on maj and min axis of m31
def plot_ellipse(ax, major_axis, minor_axis, PA, ax_small=False):
    
    PA_rad = np.deg2rad(PA) #convert to radians
    
    major_axis, minor_axis, PA_rad = major_axis.value, minor_axis.value, PA_rad.value
    theta = np.linspace(0, 2 * np.pi, 100)
    
    #Parametric equation of ellipse in polar coordinates
    r = (major_axis * minor_axis) / np.sqrt((minor_axis * np.cos(theta))**2 + (major_axis * np.sin(theta))**2)
    
    # Rotate ellipse
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    x_rot = x * np.cos(PA_rad) - y * np.sin(PA_rad)
    y_rot = x * np.sin(PA_rad) + y * np.cos(PA_rad)
    
    #Converting rotated coordinates back to polar
    r_rot = np.sqrt(x_rot**2 + y_rot**2)
    theta_rot = np.arctan2(y_rot, x_rot)
    
    #plot ellipse
    if ax_small: ax.plot(theta_rot, r_rot, label='ellipse')
    
    # Plotting major axis line
    if ax_small: 
        major_line_r = np.array([0, major_axis])
        major_line_radians = np.array([PA_rad, PA_rad])
        ax.plot(major_line_radians, major_line_r, color='purple', linestyle='--')#, label='Major Axis')
    else: 
        major_line_r = np.array([0, 30])
        major_line_radians = np.array([np.pi+PA_rad, np.pi+PA_rad])
        ax.plot(major_line_radians, major_line_r, color='purple', linestyle='-')
        # print(f"{np.rad2deg(major_line_radians)-180=}")
    
        major_line_radians = np.array([PA_rad, PA_rad])
        ax.plot(major_line_radians, major_line_r, color='purple', linestyle='-')#, label='Major Axis (Extended)')
        ax.plot([], [], color='purple', linestyle='--')#, label='Major Axis') #Making sure this sows in legend since i can find where to fix it

    #Plot minor axis line
    if ax_small: 
        minor_line_r = np.array([0, minor_axis])
        minor_line_radians = np.array([PA_rad, PA_rad])
        minor_line_radians = np.array([PA_rad + np.pi/2, PA_rad + np.pi/2]) #at an angle of 90 degrees
        ax.plot(minor_line_radians, minor_line_r, color='r', linestyle='--')#, label='Minor Axis')
    else: 
        minor_line_r = np.array([0, 30])
        minor_line_radians = np.array([PA_rad + np.pi/2 + np.pi, PA_rad + np.pi/2 + np.pi]) #at an angle of 90 degrees
        ax.plot(minor_line_radians, minor_line_r, color='r', linestyle='-')
        # print(f"{np.rad2deg(minor_line_radians)-180}=")
    
        minor_line_radians = np.array([PA_rad, PA_rad])

        minor_line_radians = np.array([PA_rad + np.pi/2, PA_rad + np.pi/2]) #at an angle of 90 degrees
        ax.plot(minor_line_radians, minor_line_r, color='r', linestyle='-')#, label='Minor Axis (Extended)')
        ax.plot([], [], color='r', linestyle='--')#,label='Minor Axis') #Making sure this sows in legend since i can find where to fix it

def convert_txt_to_Skycoord(file, withNames=False):
    ra_list, dec_list, name_list = [], [], []

    with open(file, 'r') as f:
        reader = csv.reader(f)
        RA, DEC, NAME = False, False, False
        
        for row in reader:
            if not row:  # Skip empty lines
                continue
            if row[0] == "RA": RA, DEC, NAME = True, False, False ; continue
            elif row[0] == "DEC": RA, DEC, NAME = False, True, False; continue
            elif row[0] == "NAME": RA, DEC, NAME = False, False, True; continue
            
            if RA: ra_list.append(row[0].strip())  # Strip any extra spaces
            elif DEC: dec_list.append(row[0].strip())
            elif NAME and withNames: name_list.append(row[0].strip())

    #Ensure lists are not empty
    if not ra_list or not dec_list: raise ValueError("RA or DEC data is missing in file.")
    
    skycoord = SkyCoord(ra=ra_list, dec=dec_list, unit=(u.hourangle, u.deg))
    
    if not withNames:  #If names are not requested, just return SkyCoord
        return skycoord
    else:  #Return SkyCoord object and list of names in tuple form
        return skycoord, name_list

def gauss2d(x, y, a, x0, y0, sx, sy):
    return a * np.exp(-((x - x0) / sx)**2 / 2. - ((y - y0) / sy)**2 / 2.)

def ra_dec_to_pixels(RA, Dec, **kw):
    
    if 'wcs' in kw: wcs = kw['wcs']
    else:    
        try: hdu = fits.open(kw["fitfile"])[0]
        except KeyError: hdu = fits.open(kw["filename"])[0]
    wcs = WCS(hdu.header)
    
    def world_to_pixel(ra_dec):
        ra, dec = ra_dec
        xpix, ypix, _, _ = wcs.world_to_pixel_values(ra, dec, 0, 0)
        return xpix, ypix

    X_pixels, Y_pixels = zip(*map(world_to_pixel, zip(RA, Dec)))

    return np.array(X_pixels), np.array(Y_pixels)

def print_smoothing_scale(delts, std_x, std_y, nsig):
    """
    Computes kernel width in pixels and smoothing scale in degrees.

    Parameters:
    - delt (float): Pixel scale in degrees per pixel (DELT1 or DELT2)
    - std_x (float): Standard deviation in pixels (X-axis).
    - std_y (float): Standard deviation in pixels (Y-axis).
    - nsig (float, optional): Smoothing extent factor. Default is 0.6.

    Returns:
    - dict: Kernel width in pixels and smoothing scale in degrees for X and Y.
    """
    #incase its reversed (a minus sign... especially for RA)
    delts = tuple(map(np.abs,(delts[0],delts[1]))) 

    #kernel width in pixels
    kernel_width_x = 2 * nsig * std_x
    kernel_width_y = 2 * nsig * std_y

    #smoothing scale in degrees
    smoothing_scale_x = kernel_width_x * delts[0]
    smoothing_scale_y = kernel_width_y * delts[1]

    print(f"Kernel width (X): {kernel_width_x:.3f} pixels")
    print(f"Kernel width (Y): {kernel_width_y:.3f} pixels")
    print(f"Smoothing scale (X): {smoothing_scale_x:.3f} degrees")
    print(f"Smoothing scale (Y): {smoothing_scale_y:.3f} degrees")

def smooth_2d_image(ra, dec, fitfile, imsize=5000, nsig=1):

    im = np.zeros((imsize, imsize), dtype=float)

    x0s, y0s = ra_dec_to_pixels(ra, dec, fitfile=fitfile)
    
    std_x, std_y = tuple(map(np.std,(x0s, y0s)))
    # print(f"{std_x=}")
    # print(f"{std_y=}")

    """IMPORTANT"""
    smooth = 1.1 #Factor for smoothing scatter plot via 2d guassian
    std_x, std_y = smooth*std_x, smooth*std_y
    
    DELTS = fits.open('LGSNLGSR.SQLGBB.FITS')[0].header['CDELT*'][:2]
    print_smoothing_scale(DELTS, std_x=std_x, std_y=std_y, nsig=nsig)
    
    # rm_m31_normalized = np.where(rm_m31 < 0, -1, np.where(rm_m31 > 0, 1, 0))
    
    sxs = [std_x]*len(x0s)
    sys = [std_y]*len(x0s)
    amps = rm_m31
        
    for x0, y0, sx, sy, amp in zip(x0s, y0s, sxs, sys, amps):
        xlo, xhi = int(x0 - nsig * sx), int(x0 + nsig * sx)
        ylo, yhi = int(y0 - nsig * sy), int(y0 + nsig * sy)

        xlo = max(xlo, 0)
        xhi = min(xhi, imsize)
        ylo = max(ylo, 0)
        yhi = min(yhi, imsize)

        # Generate grids of x and y coordinates
        imx, imy = np.meshgrid(np.arange(xlo, xhi), np.arange(ylo, yhi))

        # Ensure dimensions match
        if imx.size == 0 or imy.size == 0:
            continue

        # Calculate Gaussian distribution and add toimage
        im[ylo:yhi, xlo:xhi] += gauss2d(imx, imy, amp, x0, y0, sx, sy)

    return im


def get_width_midpoints(patchname): 
    
    #Transforming rectangle vertices to data coordinate system
    vert = patchname.get_path().vertices
    vertices = patchname.get_patch_transform().transform(vert) 
    
    # Calculate lengths of all edges
    edges = [
        (vertices[0], vertices[1]),
        (vertices[1], vertices[2]),
        (vertices[2], vertices[3]),
        (vertices[3], vertices[0])
    ]

    edge_lengths = [np.linalg.norm(v2 - v1) for v1, v2 in edges]

    #Retreiving indices of two shortest edges
    shortest_edge_indices = np.argsort(edge_lengths)[:2]

    #Calculating midpoints based on those indices that target 'edges'
    mid = [(edges[i][0] + edges[i][1])/2 for i in shortest_edge_indices]

    x = [mid[0][0], mid[1][0]]
    y = [mid[0][1], mid[1][1]]
    
    return x, y

def points_inOrOut_patch(p, patchname, mask, rm_s, rm_errs, In=True, **kw):
    """
    Handle points near poles, reflecting RA when declination approaches 90 or -90.
    """
    val_ra = 360  # RA range (0 to 360 degrees)
    tolerance = 1e-6  # Small tolerance to avoid boundary issues

    def _points_in_or_out_with_mirrored_RA(x, y, patch):
        # Transform vertices of patch
        patch_vertices = patch.get_path().vertices
        patch_vertices_transformed = patch.get_patch_transform().transform(patch_vertices)

        ra_vertices = patch_vertices_transformed[:, 0]  # RA vertices
        dec_vertices = patch_vertices_transformed[:, 1]  # Dec vertices

        # Determine where dec is greater than 90 or less than -90
        near_north_pole = dec_vertices > 90
        near_south_pole = dec_vertices < -90

        # Collect paths for standard points inside patch
        patch_path = Path(patch_vertices_transformed)
        x_wrapped = np.mod(x, val_ra)  # Wrap RA to [0, 360)

        # Fix points near RA=360 boundary
        x_wrapped = np.where(np.abs(x_wrapped - val_ra) < tolerance, 0, x_wrapped)

        # Initial mask: points within patch
        combined_mask = patch_path.contains_points(np.vstack((x_wrapped, y)).T)

        # If patch is near North Pole
        if np.any(near_north_pole):
            # Reflect RA points for dec > 90
            mirrored_ra_vertices = (val_ra - ra_vertices[near_north_pole]) % val_ra
            mirrored_dec_vertices = 180 - dec_vertices[near_north_pole]  # Reflect dec across 90

            mirrored_patch_vertices = np.column_stack((mirrored_ra_vertices, mirrored_dec_vertices))
            mirrored_patch_path = Path(mirrored_patch_vertices)

            # Apply mask for points within mirrored region
            mirrored_mask = mirrored_patch_path.contains_points(np.vstack((x_wrapped, y)).T)
            combined_mask = np.logical_or(combined_mask, mirrored_mask)

        # If patch is near South Pole
        if np.any(near_south_pole):
            # Reflect RA points for dec < -90
            mirrored_ra_vertices = (val_ra - ra_vertices[near_south_pole]) % val_ra
            mirrored_dec_vertices = -180 - dec_vertices[near_south_pole]  # Reflect dec across -90

            mirrored_patch_vertices = np.column_stack((mirrored_ra_vertices, mirrored_dec_vertices))
            mirrored_patch_path = Path(mirrored_patch_vertices)

            # Apply mask for points within mirrored region
            mirrored_mask = mirrored_patch_path.contains_points(np.vstack((x_wrapped, y)).T)
            combined_mask = np.logical_or(combined_mask, mirrored_mask)

        return combined_mask if In else np.logical_not(combined_mask)

    # Apply check for points inside patch and with mirrored RA
    inside_mask = _points_in_or_out_with_mirrored_RA(
        p.ra.deg[mask],
        p.dec.deg[mask],
        patchname
    )

    # Select points inside region
    p_inside = p[mask][inside_mask]

    # RM values and errors of those points
    rm_value_inside = rm_s[mask][inside_mask]
    rm_err_inside = rm_errs[mask][inside_mask]

    return p_inside, rm_value_inside, rm_err_inside

def confining_circle(ax, ra, dec, radius, polar_angle, positions, 
                pos_mask, neg_mask, color_pos="blue", color_neg="red", 
                transform="icrs", return_data=False, return_err=True,
                plot=(True, True), **kw):
    """
    Plot a circular region on a given axis and scatter positive and negative RM values within it.
    
    Parameters:
    ax : matplotlib.axes.Axes
       axis on which to plot.
    ra : float
        Right Ascension of center of circle in degrees.
    dec : float
        Declination of center of circle in degrees.
    radius : float
        Radius of circle in degrees.
    polar_angle : float
        Position angle of region in degrees (optional for future enhancement).
    positions : object
       ICRS coordinates of RM positions.
    pos_mask : array-like
        Mask to filter positive RM values.
    neg_mask : array-like
        Mask to filter negative RM values.
    color_pos : str, optional
        Color for positive RM values (default is "blue").
    color_neg : str, optional
        Color for negative RM values (default is "red").
    transform : str, optional
        Coordinate system for transformation (default is "icrs").
    return_data : bool, optional
        If True, return positions and RM values of points inside region (default is False).
    return_err : bool, optional
        If True, returns errors in positive and negative RM values inside region (default is True).
    plot : tuple of bool, optional
        If True, plot circle and scatter points (default is (True, True)).
    
    Returns:
    If return_data is True:
        pos_in_region, rm_value_pos, neg_in_region, rm_value_neg
    """
    rm_s, rm_errs = kw["rm_s"], kw["rm_errs"]
    
    # Create a circle (patch) for defining region
    circle = Circle((ra, dec), radius=radius, edgecolor='k', fill=False, linestyle='-', 
                    linewidth=1.1, transform=ax.get_transform('world'))
    
    if plot[0]: ax.add_patch(circle)  # To plot circle
    
    if plot[1]:  # To plot scatter points in red and blue depending on pos or neg
        alpha = 1
        # Getting positive positions and their RM values inside region
        pos_in_region, rm_value_pos, _ = points_inOrOut_patch(positions, circle, pos_mask, rm_s, rm_errs, typ='circle')
        ax.scatter(pos_in_region.ra, pos_in_region.dec,
                   transform=ax.get_transform(transform), marker='o', 
                   s=1, color=color_pos, alpha=alpha)
        
        # Getting negative positions and their RM values inside region
        neg_in_region, rm_value_neg, _ = points_inOrOut_patch(positions, circle, neg_mask, rm_s, rm_errs, typ='circle')
        ax.scatter(neg_in_region.ra, neg_in_region.dec,
                   transform=ax.get_transform(transform), marker='o', 
                   s=1, color=color_neg, alpha=alpha)
    
    # Get positive positions and their RM values inside region without plotting
    pos_in_region, rm_value_pos, rm_err_pos = points_inOrOut_patch(positions, circle, pos_mask, rm_s, rm_errs, typ='circle')
    # Get negative positions and their RM values inside region without plotting
    neg_in_region, rm_value_neg, rm_err_neg = points_inOrOut_patch(positions, circle, neg_mask, rm_s, rm_errs, typ='circle')
    
    if return_data:
        if return_err:  # Return errors (default)
            return (pos_in_region, rm_value_pos, rm_err_pos, 
                    neg_in_region, rm_value_neg, rm_err_neg)
        else:
            return (pos_in_region, rm_value_pos, 
                    neg_in_region, rm_value_neg)
        
def confining_rectangle(ax, ra, dec, width, height, angle, polar_angle, positions, 
                pos_mask, neg_mask, color_pos="blue", color_neg="red", 
                transform="icrs", return_data=False, return_err = True,
                plot=(True,True), **kw):

    rm_s, rm_errs = kw["rm_s"], kw["rm_errs"]
    
    # Createrectangle (patch) for defining region
    #Note that position of rectangel is at bottom corner
    #But box should be shifted so that point chosen is same as center of box
    bottom_left_x = ra - width / 2
    bottom_left_y = dec - height / 2
    rect = Rectangle((bottom_left_x, bottom_left_y), width=width, height=height,
                     angle=angle + polar_angle, transform=ax.get_transform('world'))

    if plot[1]: #to plot scatter in red and blue depending on pos or neg.
        # Plottingrectangle
        alpha = 1
        
        marker_size = kw["marker_size"] if "marker_size" in kw else 1
        
        # Getting positive positions and their RM values inside region
        pos_in_region, rm_value_pos, _ = points_inOrOut_patch(positions, rect, pos_mask, rm_s, rm_errs)
        
        # Getting negative positions and their RM values inside region
        neg_in_region, rm_value_neg, _ = points_inOrOut_patch(positions, rect, neg_mask, rm_s, rm_errs)
       
        vmin, vmax = -50, 50  # Maximum and minimum RM limits
        # Plotting RM positions inside region
        sctt = ax.scatter(np.concatenate([neg_in_region.ra.deg, pos_in_region.ra.deg]), 
                          np.concatenate([neg_in_region.dec.deg, pos_in_region.dec.deg]), 
                          transform=ax.get_transform(transform), marker='o',
                          cmap='brg_r', vmin=vmin, vmax=vmax,
                          c=np.concatenate([rm_value_neg, rm_value_pos]))  # Color by RM values
        
        rotation = -90  # Rotation for colorbar labels
        cbars_pad = 0.01  # Padding for colorbars
        labelpad = 20  # Padding for colorbar labels
        cbar_lab_size = 17  # Font size for cbar labels
        l_s = 14  # label size for ticklabels (numerical)

        # Colorbar for scatter plot (for all positive or negative points)
        sctt_cbar = plt.colorbar(sctt, ax=ax, fraction=0.06, pad=cbars_pad)
        sctt_cbar.set_label("RM scatter (capped)", rotation=rotation,
                            labelpad=labelpad, fontsize=cbar_lab_size)
        ticks = np.linspace(vmin, vmax, 11)
        ticklabels = [str(round(tick)) for tick in ticks]
        sctt_cbar.set_ticks(ticks)
        sctt_cbar.set_ticklabels(ticklabels)
        sctt_cbar.ax.tick_params(labelsize=l_s)
    
    # Getting positive positions and their RM values inside region without plotting
    pos_in_region, rm_value_pos, rm_err_pos = points_inOrOut_patch(positions, rect, pos_mask, rm_s, rm_errs)
    # Getting negative positions and their RM values inside region without plotting
    neg_in_region, rm_value_neg, rm_err_neg = points_inOrOut_patch(positions, rect, neg_mask, rm_s, rm_errs)

    h,w=.19,.4 #Minor shifts to encapuslate outer chosen points
    rect = Rectangle((bottom_left_x-h, bottom_left_y-w), width=width, height=height,
                     angle=angle + polar_angle, color='b',
                     fill=False, linestyle='--', linewidth=2,
                     transform=ax.get_transform('world'))
    if plot[0]: ax.add_patch(rect) #to plot rectangle

    if return_data:
        if return_err: #Also give errors (default outcome of this function)
            return (pos_in_region, rm_value_pos, rm_err_pos, 
                    neg_in_region, rm_value_neg, rm_err_neg, get_width_midpoints(rect))
        else:
            return (pos_in_region, rm_value_pos, 
                    neg_in_region, rm_value_neg, get_width_midpoints(rect))

def fit_and_plot_damped_sine_wave(X, Y, initial_guess, ax, color, prime):
    
    def func(x, A, gamma, omega, phi):
        return A * (np.exp(-gamma * x) * np.sin(omega * x + phi))

    popt, pcov = curve_fit(func, X, Y, p0=initial_guess, maxfev=int(1e5))
    A_opt, gamma_opt, omega_opt, phi_opt = popt

    x_fit = np.linspace(min(X), max(X), 500)
    y_fit = func(x_fit, A_opt, gamma_opt, omega_opt, phi_opt)

    ax.plot(x_fit, y_fit, color=color, 
            label = rf"RM = $({A_opt:.3g}) e^{{-({gamma_opt:.3g}) {prime}'}} \cdot \sin\left[{omega_opt:.3g} {prime}' + ({phi_opt:.3g})\right]$")

    # #Calculating residuals and reduced chi-square
    # residuals = Y - func(X, *popt)
    # chi_squared = np.sum(residuals ** 2)
    # dof = len(Y) - len(popt)  # Degrees of freedom
    # reduced_chi_square = chi_squared / dof
    # print(f"Reduced Chi-Square: {reduced_chi_square}")
    
def axis_transformation(points, RA_range, DEC_range):

    #Calclating gradient of line
    m = (np.diff(DEC_range)/ np.diff(RA_range))[0]
    
    #Calculating Rotational angle
    angle = np.arctan(m)
    
    #rotation matrix to be used for making line new x-axis
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    
    #Rotating points by using rotation matrix
    rotated_points = points @ rotation_matrix
    
    #Rotating original x-axis (DEC=0) for a range of RA values
    x_axis_original = np.array([[RA, 0] for RA in np.linspace(*RA_range, len(points))])
    rotated_x_axis = x_axis_original @ rotation_matrix #Rotating it to new axis frame (of actual coordinate system of new_y and new_x)
    try:
        return (rotated_x_axis[:,0], rotated_x_axis[:,1], 
            rotated_points[:,0], rotated_points[:,1], 
            angle)
    except IndexError: #If its a single point
        return (rotated_x_axis[:,0], rotated_x_axis[:,1], 
            rotated_points[0], rotated_points[1], 
            angle)

def dens_verses_rm(ra, dec, rm_values, rm_errors, above_zero=True):
    filename = 'LGSNLGSR.SQLGBB.FITS'
    hdu = fits.open(filename)[0]
    
    # Extracting HI density data
    HI_dens_data = hdu.data[0, 0, :, :]
    
    # Extract RM data positions and magnitudes
    RM_RA, RM_Dec = ra, dec
    
    # Converting RM positions to pixel coordinates
    RM_pix_x, RM_pix_y = ra_dec_to_pixels(RM_RA, RM_Dec, fitfile=filename)
    
    RM_pix_x = np.round(RM_pix_x).astype(int)
    RM_pix_y = np.round(RM_pix_y).astype(int)
    
    # Getting shape of HI density data array
    shape = HI_dens_data.shape
    
    # Finding minimum non-NaN HI density value
    min_HI_dens_value = np.nanmin(HI_dens_data)
    
    # Creating an empty array to store HI_dens_values
    HI_dens_values = np.full_like(rm_values, min_HI_dens_value, dtype=float)
    
    # Mask for valid coordinates of RM. Valid meaning within BT04 data
    valid_mask = (RM_pix_x >= 0) & (RM_pix_x < shape[1]) & (RM_pix_y >= 0) & (RM_pix_y < shape[0])

    # #This calculation shows 99% of RM points were rettained.    
    # print(sum(valid_mask)/len(valid_mask)); import sys; sys.exit() 
    
    # Getting valid coordinates that only within BT04 fit file
    valid_x = RM_pix_x[valid_mask]
    valid_y = RM_pix_y[valid_mask]
    
    # Fetching HI density values at those valid RM positions
    valid_HI_dens_values = HI_dens_data[valid_y, valid_x]
    
    # Handling RM Positions with NaN column density are overwritten to have minimum value possible based on fit file's minimum column density
    valid_HI_dens_values = np.where(np.isnan(valid_HI_dens_values), min_HI_dens_value, valid_HI_dens_values)
    
    # Assign back to HI_dens_values array with overwritten valid_HI_dens_values
    HI_dens_values[valid_mask] = valid_HI_dens_values
    
    result_tuples = list(zip(rm_values, HI_dens_values))
    
    # Extract rm and HI_dens_values from result_tuples
    rm_values, HI_dens_values = zip(*result_tuples)
    
    # #this shows that 87% of RM within BT04 data was flagged as Nan.
    # count = np.sum(np.array(HI_dens_values)< 14)
    # print(100*(1 - count/sum(valid_mask))); import sys; sys.exit()
    
    # print(min(HI_dens_values), max(HI_dens_values)); import sys; sys.exit()
    fmt = "bs" if above_zero else "rs"
    plt.errorbar(HI_dens_values, rm_values, yerr=rm_errors, fmt=fmt, 
                 markersize=3, capsize=3, ecolor='k')
    
    A = 35
    f = 0.206
    phi = 29.5
    
    n = np.linspace(min(HI_dens_values), max(HI_dens_values), int(1e5))
    
    if above_zero in [True, None]:
        plt.plot(n, -A * np.sin(2 * np.pi * f * n - phi) + 65)#,
                 # label=r"$(-A) \cdot \sin\left(2\pi f \cdot x - \phi\right) +65$")
    if above_zero in [False, None]:
        plt.plot(n, A * np.sin(2 * np.pi * f * n - phi) - 104)#,
                 # label=r"$A \cdot \sin\left(2\pi f \cdot x - \phi\right) -104$")
    
    # string = f'A = {A}, f = {f}, ' + r"$\phi=$ " + f"${phi}$"
    # plt.text(0.6, 1.03, string, transform=plt.gca().transAxes, fontsize=14,
    #          bbox=dict(facecolor='none', edgecolor='none', boxstyle='round,pad=0'))

def characterize_densVersesRM():
    plt.xlabel(r'$\log_{10}(N_H)$' + ' ' + '$[cm^{-2}]$', fontsize=14, rotation='horizontal')
    plt.ylabel(r'RM [rad/$m^2$]', fontsize=14)

    #R. Braun and D. A. Thilker: Local Group HI in WSRT wide-field survey. II FIGURE 5
    #log(N_H) limited to > 17.3
    max_y, min_y = 100, -160
    plt.axvline(17.3, ymin=min_y, ymax=max_y, color="g",
                linestyle="--")#, label = "WSRT sensitivity limit")
    plt.axvline(18.3, ymin=min_y, ymax=max_y, color = "r", 
                linestyle="--")#, label = "HIPASS sensitivity limit")
    plt.axvline(17.59, ymin=min_y, ymax=max_y, color = "b", 
                linestyle="--")#, label = "GBT sensitivity limit")
    plt.fill_betweenx([min_y, max_y], 17.3, 21, color ="g", alpha= 0.3)
    plt.fill_betweenx([min_y, max_y], 18.3, 21, color ="r", alpha= 0.5)
    plt.fill_betweenx([min_y, max_y], np.log10(3.9e17), 21, color ="blue", alpha= 0.3)
    
    plt.grid(True)
    plt.xlim(14,21)
    plt.ylim(min_y,max_y)

    plt.legend(fontsize=14, loc='upper center', bbox_to_anchor=(0.35, 1.3), 
           framealpha=0, ncol=2)

    print('successfully printed correlation')
    plt.tight_layout()

#make sure to change title if you change data being given by confining_rectangle() function.
def assessing_boxScatter_with_HI(title):
    # OpeningFITS file and extracting data
    filename = 'LGSNLGSR.SQLGBB.FITS'
    hdu = fits.open(filename)[0]
    
    # Setting upWCS projection usingheader
    wcs = WCS(hdu.header)
    
    positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0
    
    # Creating figure and axis with WCS projection
    fig, ax = plt.subplots(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y', 0, 0))
    
    #Region offilament between m31 and m33
    (filament_region_coord_pos,
      rm_val_filament_pos,
      rm_err_filament_pos,
      filament_region_coord_neg, #To be used for later[] after plot() is called
      rm_val_filament_neg,
      rm_err_filament_neg,
      rect_widthMid  #mid points of width of rectangular patch
      ) = confining_rectangle(ax, 
                            rm_s = rm_m31,
                            rm_errs = err_m31,
                            ra=m31_pos.ra.deg + 16,
                            dec=m31_pos.dec.deg - 3,
                            width=12, height=30, angle=7, 
                            polar_angle=m31_pa.value, 
                            positions=rm_pos_icrs, 
                            pos_mask=positive_mask, 
                            neg_mask = negative_mask, 
                            return_data=True)
                              
    plt.clf()
    plt.figure(figsize=(10, 6))
    plt.title(title, fontsize=28, pad=80)
    dens_verses_rm(ra=filament_region_coord_pos.ra, 
                   dec=filament_region_coord_pos.dec,
                   rm_values=rm_val_filament_pos,
                   rm_errors=rm_err_filament_pos
                   )
    dens_verses_rm(ra=filament_region_coord_neg.ra, 
                   dec=filament_region_coord_neg.dec,
                   rm_values=rm_val_filament_neg,
                   rm_errors=rm_err_filament_neg,
                   above_zero=False
                   )
    characterize_densVersesRM()
    plt.show()
    
def binned_scatter(x, y, bins):
    bin_edges = np.linspace(np.min(x), np.max(x), bins + 1)
    bin_means = []
    bin_stds = []
    for i in range(len(bin_edges) - 1):
        bin_mask = (x >= bin_edges[i]) & (x < bin_edges[i + 1])
        if np.any(bin_mask):
            bin_mean_x = np.mean(x[bin_mask])
            bin_mean_y = np.mean(y[bin_mask])
            bin_std_y = np.std(y[bin_mask])
            bin_means.append([bin_mean_x, bin_mean_y])
            bin_stds.append(bin_std_y)
    return np.array(bin_means), np.array(bin_stds), bin_edges

