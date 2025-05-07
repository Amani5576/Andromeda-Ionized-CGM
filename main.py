#Note that this script only declares variables and functions. No ouputs are produced
from astropy.table import Table
from astropy import units as u
import numpy as np
import csv
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from scipy import stats
from astropy.coordinates import SkyCoord
from matplotlib.patches import Rectangle, Circle, Path, Patch
from scipy.optimize import curve_fit
from astropy.utils.exceptions import AstropyWarning
import warnings
from scipy.interpolate import RectBivariateSpline, griddata
import os
import random

from Parsers import args #Personally made parsers

b_upper_lim = 5 * u.deg if (args.b_limit or args.b_limit_small) else 0 * u.deg

def galactic_cut(RM_lat, RM_lon, rm, rm_err, b_limit=False, extra=None, even_more=None, **kw):
    """
    Apply a Galactic latitude cut to data based on |b| > b_upper_lim.
    Includes background (but done separately if provided).
    
    Parameters:
    - RM_lat: array of Galactic latitudes
    - RM_lon: array of Galactic longitudes
    - rm: array of rotation measures
    - rm_err: array of errors on rotation measures
    - b_limit: bool, if True applies a cut of |b| > 5 degrees
    - extra: optional, tuple of separation arrays for RM and background
    
    Returns:
    - Filtered RM_lat, RM_lon, rm, rm_err arrays (and background/extra if provided)
    """

    def applying_mask(*data):
        RM_lat, *_ = data
        mask = (np.abs(RM_lat) > b_upper_lim.value)
        return [d[mask] for d in data if d is not None]
    
    #if even_more is not None:
    #    others = applying_mask(RM_lat, even_more, b_limit=b_limit)[1:]

    # print(f"{len(*(extra[:len(extra)//2] if extra else []))=}")
    #print(f"{extra[len(extra)//2]=}") if extra else None
    masked_main = applying_mask(RM_lat, RM_lon, rm, rm_err, *(extra[:len(extra)//2] if extra else []))

    if "BG_data" in kw:
        if extra:
            masked_bg = applying_mask(*kw["BG_data"], *(extra[len(extra)//2:] if extra else []))
        else:
            masked_bg = applying_mask(*kw["BG_data"])
            
        return (*masked_main, *masked_bg)

    return masked_main

#Suppress specific Astropy warnings
warnings.simplefilter('ignore', AstropyWarning)

#to be used for scaling size of scatter based on their rm value
sizes = lambda vals: 30 #np.abs(vals) * 2.3

R_vir = 300.*u.kpc #Virial Radius of Andromeda in kiloparsecs (kpc)
cutoff = 30.*u.deg #(in deg) Limit for taking to account background correction
d_m31 = 780.*u.kpc #Distance to Andromeda in kiloparsecs (kpc)

m31_maj = 73.*u.arcmin #Major axis diameter from POSS in arcmin
m31_min = 45.*u.arcmin #Minor axis diameter from POSS in arcmin
m31_pa = 37.7*u.deg #PA CCW from North from de Vaucouleurs et al. 1958.  

#Virial radius angular size
L_m31 = (np.arctan(R_vir / d_m31) * u.rad.to(u.deg)).value

#Position of M31 (in ICRS coordinates)
m31_pos = SkyCoord(ra = "00:42:44.35", dec = "+41:16:08.6", unit = (u.hourangle, u.deg), frame = 'icrs')
m31_lat = m31_pos.transform_to('galactic').b

def get_data_from_catalogue(sigma_detect_limit):
    """
    Read RM catalogue, apply detection limit and galactic cut,
    transform coordinates to ICRS and M31 frame,
    and compute separations and position angles.

    Returns:
        RM_lat, RM_lon, rm, rm_err,
        position, eq_pos, rm_m31_coord,
        m31_sep, m31_theta,
        cloud6_pos, m33_pos, m33_m31coord, m33_sep, m33_theta
    """
    t = Table.read("catalog.dat", format="ascii")
    t.rename_column("col9", "l")
    t.rename_column("col10", "b")
    t.rename_column("col17", "RM")
    t.rename_column("col18", "e_RM")
    t.keep_columns(["l", "b", "RM", "e_RM"])

    RM_unit = u.rad / u.m**2
    sig_mask = ~(sigma_detect_limit * np.abs(t["e_RM"]) > np.abs(t["RM"]))
    t["l"].unit = u.deg
    t["b"].unit = u.deg
    t["RM"].unit = RM_unit
    t["e_RM"].unit = RM_unit

    rm_raw = t["RM"]
    rm_err_raw = t["e_RM"]
    RM_lon = t["l"][sig_mask]
    RM_lat = t["b"][sig_mask]
    rm = rm_raw[sig_mask]
    rm_err = rm_err_raw[sig_mask]

    # galactic_cut to be applied by caller or in next function
    position = SkyCoord(RM_lon, RM_lat, unit="deg", frame="galactic")
    eq_pos = position.transform_to('icrs')
    new_frame_of_reference = m31_pos.skyoffset_frame()
    rm_m31_coord = eq_pos.transform_to(new_frame_of_reference)

    m31_sep = eq_pos.separation(m31_pos)
    m31_theta = eq_pos.position_angle(m31_pos)

    cloud6_pos = SkyCoord("01:08:29.6 +37:45:00", unit=(u.hourangle, u.deg), frame='icrs')

    m33_pos = SkyCoord("23.462042 30.660222", unit="deg", frame='icrs')
    m33_m31coord = m33_pos.transform_to(new_frame_of_reference)
    m33_sep = m33_m31coord.separation(m31_pos)
    m33_theta = m33_m31coord.position_angle(m31_pos)

    return (RM_lat, RM_lon, rm, rm_err,
            position, eq_pos, rm_m31_coord,
            m31_sep, m31_theta, new_frame_of_reference,
            cloud6_pos, m33_pos, m33_m31coord, m33_sep, m33_theta)


def get_CGM_and_BG_masks(rm_m31_coord, eq_pos, m31_sep, elliptic_CGM, elliptic_CGM_bg):
    """
    Compute boolean masks for CGM and background region around M31.
    Returns:
        m31_condition, bg_condition
    """
    if elliptic_CGM:
        x = rm_m31_coord.lon.arcmin
        y = rm_m31_coord.lat.arcmin
        theta = m31_pa.to(u.rad).value
        x_rot = x * np.cos(theta) + y * np.sin(theta)
        y_rot = -x * np.sin(theta) + y * np.cos(theta)
        major_deg = L_m31
        minor_deg = major_deg * (m31_min / m31_maj).value
        a = major_deg * 60
        b = minor_deg * 60
        ellipse = lambda a_, b_: (x_rot**2 / a_**2 + y_rot**2 / b_**2)
        m31_condition = ellipse(a, b) <= 1
    else:
        m31_condition = m31_sep.deg <= L_m31

    if elliptic_CGM_bg:
        bg_extent = ((cutoff.value - L_m31) * u.deg).to(u.arcmin).value
        a_outer = a + bg_extent
        b_outer = b + bg_extent
        outer_condition = ellipse(a_outer, b_outer) <= 1
        inner_condition = ellipse(a, b) >= 1
        bg_condition = inner_condition & outer_condition
    else:
        bg_condition = (m31_sep.deg > L_m31) & (m31_sep.deg < cutoff.value)

    return m31_condition, bg_condition


def apply_CGM_and_BG_masks(rm_m31_coord, eq_pos, position,
                           rm, rm_err, m31_sep,
                           m31_condition, bg_condition):
    """
    Apply masks to extract positions and RM values for CGM and background.
    Returns masked variables.
    """
    bg_pos = rm_m31_coord[bg_condition]
    bg_pos_icrs = bg_pos.transform_to('icrs')
    rm_pos = rm_m31_coord[m31_condition]
    rm_pos_icrs = rm_pos.transform_to('icrs')

    rm_pos_gal_lat = position.b.deg[m31_condition]
    rm_pos_gal_lat_bg = position.b.deg[bg_condition]

    rm_bg = rm[bg_condition]
    m31_sep_bg = (m31_sep.deg[bg_condition]) * u.deg
    err_bg = rm_err[bg_condition]

    rm_m31 = rm[m31_condition]
    m31_sep_Rvir = (m31_sep.deg[m31_condition]) * u.deg
    err_m31 = rm_err[m31_condition]

    return (bg_pos, bg_pos_icrs, rm_pos, rm_pos_icrs,
            rm_pos_gal_lat, rm_pos_gal_lat_bg,
            rm_bg, m31_sep_bg, err_bg,
            rm_m31, m31_sep_Rvir, err_m31)

# 1)Loading & transforming catalogue, get all raw coords + RMs + M31/M33 stuff
(RM_lat, RM_lon, rm, rm_err, position, eq_pos, rm_m31_coord, 
 m31_sep, m31_theta, new_frame_of_reference,
 cloud6_pos, m33_pos, m33_m31coord, m33_sep, m33_theta
) = get_data_from_catalogue(sigma_detect_limit=args.sig_limit)

# 2)Building CGM / BG masks (ellipse vs circle controlled by flags)
m31_condition, bg_condition = get_CGM_and_BG_masks(
    rm_m31_coord, eq_pos, m31_sep,
    elliptic_CGM=args.elliptic_CGM,
    elliptic_CGM_bg=args.elliptic_CGM_bg
)

# 3)Apply those masks to slice out CGM & BG subsets
(bg_pos, bg_pos_icrs, rm_pos, rm_pos_icrs, rm_pos_gal_lat, 
 rm_pos_gal_lat_bg, rm_bg, m31_sep_bg, err_bg, rm_m31, 
 m31_sep_Rvir, err_m31) = apply_CGM_and_BG_masks(
    rm_m31_coord, eq_pos, position, rm, rm_err, 
    m31_sep, m31_condition, bg_condition
)
    
#Calculating polar angles (#More variable declarations)
shift = 180*u.deg
PA_bg = m31_theta[bg_condition] + shift
PA_rm = m31_theta[m31_condition] + shift
PA_m33 = m33_theta + shift

print(f"{args.bg_corr=}")
def BG_correction(rm_coords, rm_values, bg_coords, bg_values, bg_corr=args.bg_corr, **kw):
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
    bg_corr   : bool
        Returns uncorrected RM if False (given by argument parsing --bg-corr)
    Returns:
    --------
    np.ndarray
        Background-corrected RM values.
    """
    
    if not bg_corr: #no Backgorund correction with spline fit
        return rm_values
    
    Spline_order = 1

    #Converting to degrees
    x_bg = bg_coords.ra.deg  #(N,)
    y_bg = bg_coords.dec.deg  #(N,)
    
    rm_x = rm_coords.ra.deg  #(M,)
    rm_y = rm_coords.dec.deg  #(M,)
    
    #Define a regular grid for interpolation
    grid_res = kw.get("grid_res", 50) #len(x_bg)*1 #"N" #true modifying variable... Spline_order really doesnt do much...
    x_grid = np.linspace(x_bg.min(), x_bg.max(), grid_res) #(N,)
    y_grid = np.linspace(y_bg.min(), y_bg.max(), grid_res) #(N,)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid) #each having dimensions (N,N)
    
    #Ensuring griddata inputs have correct dimensions
    bg_points = np.column_stack((x_bg, y_bg))  #(N,2) format required by griddata
    grid_points = np.column_stack((X_grid.ravel(), Y_grid.ravel()))  #(N*N,2)
    
    #Interpolating background values onto grid
    bg_grid = griddata(bg_points, bg_values, grid_points, method='cubic')  #(N*N,)
    bg_grid = bg_grid.reshape(X_grid.shape)  #Reshape back to (N,N)
    
    #Handle NaNs (replace with nearest-neighbor interpolation)
    if np.isnan(bg_grid).any():
        bg_grid = griddata(bg_points, bg_values, grid_points, method='nearest')
        bg_grid = bg_grid.reshape(X_grid.shape)

    #Fitting spline to interpolated background RM data
    fbeam = RectBivariateSpline(np.sort(y_grid), np.sort(x_grid), bg_grid,
                                kx=Spline_order, ky=Spline_order)
    
    #Interpolating background RM values at RM positions
    bg_values_interp = fbeam.ev(rm_y, rm_x)
    
    #Finally subtracting complex background (interpolated) from given rm coords.
    rm_corrected = rm_values - bg_values_interp
    
    return rm_corrected

#Conducting background subtraction.  
rm_m31 = BG_correction(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg)
rm_bg = BG_correction(bg_pos_icrs, rm_bg, bg_pos_icrs, rm_bg) #Doing same for background itself
#Validation/Sanity Check
#print(len(rm_m31), len(err_m31), len(np.power(err_m31, 2)), len(rm_bg))

#print(f"{len(m31_sep_Rvir)=}")
#Note that positions of RM must be changed based on use of function and will be changed back to ICRS coordinates
#imposing |b|>5 deg jsut after spline correction only for M31 since b_limit_small is used.
#Note, function also takes care of similar format for background region.

(rm_pos_gal_lat, rm_m31_lon, rm_m31, err_m31, m31_sep_Rvir, PA_rm,
 rm_pos_gal_lat_bg, rm_bg_lon, rm_bg, err_bg, m31_sep_bg, PA_bg
                                                                ) = galactic_cut(RM_lat=rm_pos_gal_lat,
                                                                        RM_lon=rm_pos_icrs.galactic.l.deg, 
                                                                        rm=rm_m31, rm_err=err_m31,
                                                                        BG_data= (
                                                                            rm_pos_gal_lat_bg,
                                                                            bg_pos_icrs.galactic.l.deg, 
                                                                            rm_bg, err_bg
                                                                        ),
                                                                        extra=( #Other data that needs be stripped to have same number of data points
                                                                            m31_sep_Rvir, PA_rm,
                                                                            m31_sep_bg, PA_bg
                                                                        ),
                                                                        b_limit=args.b_limit_small
                                                                        )
#print(f"{rm_pos_gal_lat[0]=}"); import sys; sys.exit()
rm_pos_gal = SkyCoord(l=rm_m31_lon, b=rm_pos_gal_lat, unit="deg", frame="galactic")
bg_pos_gal = SkyCoord(l=rm_bg_lon, b=rm_pos_gal_lat_bg, unit="deg", frame="galactic")
rm_pos_icrs = rm_pos_gal.icrs
bg_pos_icrs = bg_pos_gal.icrs
# print(f"{len(m31_sep_Rvir)=}") 

#IMPORTANT
bin_num = 10

#Calculate mean of RM values within Rvir of M31
bin_means, bin_edges_mean_m31, binnumber = stats.binned_statistic(m31_sep_Rvir, rm_m31, statistic = 'mean', bins = bin_num)

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

bin_med, bin_edges_median_m31, binnumber = stats.binned_statistic(m31_sep_Rvir, rm_m31, statistic = 'median', bins = bin_num)

bin_width = (bin_edges_mean_m31[1] - bin_edges_mean_m31[0])
bin_centers = bin_edges_mean_m31[1:] - bin_width/2

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

#SEMed = lambda sig, N: 1.253*sig/N #Using a bootstrap method of calculating standard error of median

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


#=============================================================================
#M33 Radial distance from M31 --> d_m33
#M33 Angular distance from M31 --> m33_sep 
#Background RM Radial distance from M31 --> d_bg
#Background RM Angular distance from M31 --> m31_sep_bg
#R_vir RM Radial distance from M31 --> d_rm 
#R_vir RM Angular distance from M31 --> m31_sep_Rvir
#=============================================================================


#following variables below are for M33 Aanlysis alone (Not interpolated from M31)

"""Given by: https://watermark.silverchair.com/sty1946.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA1MwggNPBgkqhkiG9w0BBwagggNAMIIDPAIBADCCAzUGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMZOR9nxmNwLQ-n5xaAgEQgIIDBkOA2_nPRhEoUdyVWWdH7HjcxkZI91FVxVd10MSR2B8JJzz4o19U3CHtrPc6xBwZV5tkx5Rn_QxT82lUN9jD7TMupO4NOSNc-QMHGiXmwMKV3QD9Q742oJ2-nsU6LP3-LDkfkuj9IlIWs1iBz3VwHJv1uUJw1ec5gomXmkmSfC338C5SuOM05Iifv7RzQON3BAPcK4kyL_L4imLPEKkIbtPvNH24gxal7OyYKvfQ6G0wnPrNFBJUAOYvpApZkq-YBpcFp35zPstuBI6pR_ngTAlU5PBY3yXsXyKZp9BYzcICEu1cUkgIELIb-2OAPk2PvytZfXPqmzq6G_fURukCPKje1zFO03UZ5_PSr1VNPv5WPNTIkXdzRGnkivyeQlg48b3skkv5Ox9yJH8qr_TOe8qmZf2fP_eDO14RVsEfu6Is_Kx7AWEbdTaRukvWNnbYsTWMZNDn8d-jcYDtSoAxDVF8y8TrYVFS6prrip_ZmnBABoAjJ2O91BV4nddoTub3XsSkWDhStt9_N-Su0FxZFLjyHTcYSZUf4nW4nykzAl5KlyHdtcQguv8b5HWuBNoJ1867irDVGNDzHPcksBmxA0ftVBvlvCtJ8TCPq5Gux0EZctAMzL_Qv7_Pwg7uDrqt9d7mBIdNKUgDbUIi4ssBarnv4vWFoN5TnbaRf0VUjxaw8pWHGKqMyXotZWWSdy9tpyd2d52lG4itwUATBhbPl9pbGfYBM0GR8rrS5eYjcOyMnk_VozW23Nc5vI3QO8NfxZ8Jd9ETQEOzXUzsa0ppsSK6AQbir3h1dGGCVMzxyyYNQM1je14NO98daWiRozsKXK6l8r8BsyOKBOqagjBDSCfC-5DGGGUZ1jx1ZICv_776aLw3b6iEv9awF6C6DUfbFxE1yntrG0cRTgZ11kYrB8ABA4LXBmUBV7FXl-f7gsNIvuhklBqIOFUXDv-R_pBeVu6K6KN1N26AzHgR-Is4-wU6TpS3ugdJRpIMezkZ9fgzVX4q5Aqlj64jLHilUeciIDnXeUr2xQ
On page 1885 (Figure 1)"""
d_m33 = 794 * u.kpc
R_vir_m33 = 160 * u.kpc

#Virial radius angular size
L_m33 = (np.arctan(R_vir_m33 / d_m33) * u.rad.to(u.deg)).value

cutoff_m33 = cutoff * (R_vir_m33 / R_vir)

def vars_before_correction_m33(eq_pos, rm, rm_err, position, R_vir_m33, d_m33, cutoff_m33):
    
    m33_pos = SkyCoord("23.462042 30.660222", unit="deg", frame='icrs')
    new_frame_of_reference = m33_pos.skyoffset_frame()
    rm_m33_coord = eq_pos.transform_to(new_frame_of_reference)

    #M33 separations
    m33_sep_new = eq_pos.separation(m33_pos)
    m33_theta = eq_pos.position_angle(m33_pos)


    #Conditions
    bg_condition_m33 = (m33_sep_new.deg > L_m33) & (m33_sep_new.deg < cutoff_m33.value)
    m33_condition = m33_sep_new.deg <= L_m33

    #Filtered positions
    bg_pos_m33 = rm_m33_coord[bg_condition_m33]
    bg_pos_icrs_m33 = bg_pos_m33.transform_to("icrs")
    rm_pos_m33 = rm_m33_coord[m33_condition]
    rm_pos_icrs_m33 = rm_pos_m33.transform_to("icrs")

    #Galactic latitudes
    rm_pos_gal_lat_m33 = position.b.deg[m33_condition]
    rm_pos_gal_lat_bg_m33 = position.b.deg[bg_condition_m33]

    #RM and error values
    rm_bg_m33 = rm[bg_condition_m33]
    m33_sep_bg = (m33_sep_new.deg[bg_condition_m33]) * u.deg
    err_bg_m33 = rm_err[bg_condition_m33]

    rm_m33 = rm[m33_condition]
    m33_sep_Rvir = (m33_sep_new.deg[m33_condition]) * u.deg
    err_m33 = rm_err[m33_condition]

    return (
        new_frame_of_reference, rm_m33_coord,
        m33_sep_new, m33_theta,
        bg_condition_m33, m33_condition,
        bg_pos_m33, bg_pos_icrs_m33, rm_pos_m33, rm_pos_icrs_m33,
        rm_pos_gal_lat_m33, rm_pos_gal_lat_bg_m33,
        rm_bg_m33, m33_sep_bg, err_bg_m33,
        rm_m33, m33_sep_Rvir, err_m33
    )

(new_frame_of_reference, rm_m33_coord,
        m33_sep_new, m33_theta,
        bg_condition_m33, m33_condition,
        bg_pos_m33, bg_pos_icrs_m33, rm_pos_m33, rm_pos_icrs_m33,
        rm_pos_gal_lat_m33, rm_pos_gal_lat_bg_m33,
        rm_bg_m33, m33_sep_bg, err_bg_m33,
        rm_m33, m33_sep_Rvir, err_m33
    ) = vars_before_correction_m33(eq_pos, rm, rm_err, position, R_vir_m33, d_m33, cutoff_m33)

rm_m33 = BG_correction(rm_pos_icrs_m33, rm_m33, bg_pos_icrs_m33, rm_bg_m33)
rm_bg_m33 = BG_correction(bg_pos_icrs_m33, rm_bg_m33, bg_pos_icrs_m33, rm_bg_m33) #Doing same for background itself

#Imposing |b|>5 deg, right after Spline fitting for m33 (similar to M31)
#print(f"{len(m33_sep_Rvir)=}")
(rm_pos_gal_lat_m33, rm_m33_lon, rm_m33, err_m33, m33_sep_Rvir,
 rm_pos_gal_lat_bg_m33, rm_bg_lon_m33, rm_bg_m33, err_bg_m33, m33_sep_bg
                                                            )= galactic_cut(RM_lat=rm_pos_gal_lat_m33,
                                                                            RM_lon=rm_pos_icrs_m33.galactic.l.deg, 
                                                                            rm=rm_m33, rm_err=err_m33,
                                                                            BG_data= (
                                                                                bg_pos_icrs_m33.galactic.b.deg,
                                                                                bg_pos_icrs_m33.galactic.l.deg, 
                                                                                rm_bg_m33, err_bg_m33
                                                                                ),
                                                                            extra=( #Other data that needs be stripped to have same number of data points
                                                                                m33_sep_Rvir,
                                                                                m33_sep_bg
                                                                                ),
                                                                                b_limit=args.b_limit_small
                                                                                )
#print(f"{len(m33_sep_Rvir)=}") ; import sys; sys.exit()
rm_pos_gal_m33 = SkyCoord(l=rm_m33_lon, b=rm_pos_gal_lat_m33, unit="deg", frame="galactic")
bg_pos_gal_m33 = SkyCoord(l=rm_bg_lon_m33, b=rm_pos_gal_lat_bg_m33, unit="deg", frame="galactic")
rm_pos_icrs_m33 = rm_pos_gal.icrs
bg_pos_icrs_m33 = bg_pos_gal.icrs

d_bg_m33 = get_projected_d(m33_sep_bg, d_m33)
d_rm_m33 = get_projected_d(m33_sep_Rvir, d_m33) #Is only for within R_vir

bin_num_m33 = bin_num #Same number of bins as for M31

#Calculate mean of RM values within Rvir of M33
bin_means_m33, bin_edges_m33, binnumber_m33 = stats.binned_statistic(m33_sep_Rvir, rm_m33, statistic = 'mean', bins = bin_num_m33)

#Using standard error of mean for error bars: 
bin_std_m33, bin_edges_m33, binnumber_m33 = stats.binned_statistic(m33_sep_Rvir, 
    rm_m33, 
    statistic = lambda rm_m33: stats.sem(rm_m33), #Standard Error ofMean
    bins = bin_num_m33)  

def plot_m31_stats(ax, d_bin_cent=d_bin_centers, **kw):

    which = kw.get("title", "both")
    color = kw.get("color", "r")
    alpha = 1
    #Extract values from Astropy Quantity objects if they exist
    d_bin_cent_values = d_bin_cent.value if hasattr(d_bin_cent, 'value') else d_bin_cent
    bin_means_values = np.absolute(bin_means.value) if hasattr(bin_means, 'value') else np.absolute(bin_means)
    bin_med_values = np.absolute(bin_med.value) if hasattr(bin_med, 'value') else np.absolute(bin_med)
    bin_std_values = bin_std.value if hasattr(bin_std, 'value') else bin_std

    #Determine markeredgecolor based on color value
    if color == "white": markeredgecolor = "k"
    elif color == "k": markeredgecolor = "white"
    else: markeredgecolor = "k" if which != "both" else None

    #Save Mean values to file
    if which == "both" or which == "Mean":
        ax.errorbar(d_bin_cent_values, bin_means_values, yerr=bin_std_values, 
                    color=color, fmt='.-', alpha=alpha, label="M31 mean", capsize=3, 
                    markeredgecolor=markeredgecolor
                   )
        #Save Mean values in a text file with appropriate units
        mean_filename = f"m31_mean_values_{bin_num}_bins.txt"
        header_mean = "Bin Center [kpc], Mean |RM| Value [rad/m^2], |RM| err [rad/m^2]"
        np.savetxt(mean_filename, np.column_stack((d_bin_cent_values, bin_means_values, bin_std_values)),
                   header=header_mean, fmt='%0.6f')

    #Save Median values to file
    if which == "both" or which == "Median":
        ax.errorbar(d_bin_cent_values, bin_med_values, yerr=bin_std_values, 
                    color=kw.get("color", 'r' if "title" in kw else 'orange'), fmt='.-', capsize=3, 
                    markeredgecolor=markeredgecolor, alpha=alpha, label="M31 median")
        #Save Median values in a text file with appropriate units
        median_filename = f"m31_median_values_{bin_num}_bins.txt"
        header_median = "Bin Center [kpc], Median |RM| Value [rad/m^2], |RM| err [rad/m^2]"
        np.savetxt(median_filename, np.column_stack((d_bin_cent_values, bin_med_values, bin_std_values)),
                   header=header_median, fmt='%0.6f')

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
    
    #Rotate ellipse
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    x_rot = x * np.cos(PA_rad) - y * np.sin(PA_rad)
    y_rot = x * np.sin(PA_rad) + y * np.cos(PA_rad)
    
    #Converting rotated coordinates back to polar
    r_rot = np.sqrt(x_rot**2 + y_rot**2)
    theta_rot = np.arctan2(y_rot, x_rot)
    
    #plot ellipse
    if ax_small: ax.plot(theta_rot, r_rot, label='ellipse')
    
    #Plotting major axis line
    if ax_small: 
        major_line_r = np.array([0, major_axis])
        major_line_radians = np.array([PA_rad, PA_rad])
        ax.plot(major_line_radians, major_line_r, color='purple', linestyle='--')#, label='Major Axis')
    else: 
        major_line_r = np.array([0, 30])
        major_line_radians = np.array([np.pi+PA_rad, np.pi+PA_rad])
        ax.plot(major_line_radians, major_line_r, color='purple', linestyle='-')
        #print(f"{np.rad2deg(major_line_radians)-180=}")
    
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
        #print(f"{np.rad2deg(minor_line_radians)-180}=")
    
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
            if not row:  #Skip empty lines
                continue
            if row[0] == "RA": RA, DEC, NAME = True, False, False ; continue
            elif row[0] == "DEC": RA, DEC, NAME = False, True, False; continue
            elif row[0] == "NAME": RA, DEC, NAME = False, False, True; continue
            
            if RA: ra_list.append(row[0].strip())  #Strip any extra spaces
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

def get_smoothing_scale(delts, std_x, std_y, nsig):
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

    return smoothing_scale_x, smoothing_scale_y #For file naming when saving plot

def smooth_2d_image(ra, dec, fitfile, rm_m31, imsize=5000, nsig=1):

    im = np.zeros((imsize, imsize), dtype=float)

    x0s, y0s = ra_dec_to_pixels(ra, dec, fitfile=fitfile)

    sft = .8 #Parametre for lightly shifting position
    x0s += sft; y0s += sft #Slight shift to be in sync with wcs axis wehn plotting RM_smoothing.py

    std_x, std_y = tuple(map(np.std,(x0s, y0s)))
    
    #Making them same from now on: 
    std = max(np.std(x0s), np.std(y0s))
    std_x = std_y = std
    # print(f"{std_x=}")
    # print(f"{std_y=}")

    """IMPORTANT"""
    if args.elliptic_CGM:
        smooth = 1.2
        #For half major axis
        # smooth = 2.5 #Factor for smoothing scatter plot via 2d guassian
    else:
        smooth = 1.1 #Factor for smoothing scatter plot via 2d guassian

    std_x, std_y = smooth*std_x, smooth*std_y
    
    DELTS = fits.open('LGSNLGSR.SQLGBB.FITS')[0].header['CDELT*'][:2]
    kernel_deg = get_smoothing_scale(DELTS, std_x=std_x, std_y=std_y, nsig=nsig)
    
    #rm_m31_normalized = np.where(rm_m31 < 0, -1, np.where(rm_m31 > 0, 1, 0))
    
    sxs = [std_x]*len(x0s)
    sys = [std_y]*len(y0s)
    amps = rm_m31
        
    for x0, y0, sx, sy, amp in zip(x0s, y0s, sxs, sys, amps):
        xlo, xhi = int(x0 - nsig * sx), int(x0 + nsig * sx)
        ylo, yhi = int(y0 - nsig * sy), int(y0 + nsig * sy)

        xlo = max(xlo, 0)
        xhi = min(xhi, imsize)
        ylo = max(ylo, 0)
        yhi = min(yhi, imsize)

        #Generate grids of x and y coordinates
        imx, imy = np.meshgrid(np.arange(xlo, xhi), np.arange(ylo, yhi))

        #Ensure dimensions match
        if imx.size == 0 or imy.size == 0:
            continue

        #Calculate Gaussian distribution and add toimage
        im[ylo:yhi, xlo:xhi] += gauss2d(imx, imy, amp, x0, y0, sx, sy)

    return im, kernel_deg

def get_width_midpoints(patchname): 
    
    #Transforming rectangle vertices to data coordinate system
    vert = patchname.get_path().vertices
    vertices = patchname.get_patch_transform().transform(vert) 
    
    #Calculate lengths of all edges
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
    val_ra = 360  #RA range (0 to 360 degrees)
    tolerance = 1e-6  #Small tolerance to avoid boundary issues

    def _points_in_or_out_with_mirrored_RA(x, y, patch):
        #Transform vertices of patch
        patch_vertices = patch.get_path().vertices
        patch_vertices_transformed = patch.get_patch_transform().transform(patch_vertices)

        ra_vertices = patch_vertices_transformed[:, 0]  #RA vertices
        dec_vertices = patch_vertices_transformed[:, 1]  #Dec vertices

        #Determine where dec is greater than 90 or less than -90
        near_north_pole = dec_vertices > 90
        near_south_pole = dec_vertices < -90

        #Collect paths for standard points inside patch
        patch_path = Path(patch_vertices_transformed)
        x_wrapped = np.mod(x, val_ra)  #Wrap RA to [0, 360)

        #Fix points near RA=360 boundary
        x_wrapped = np.where(np.abs(x_wrapped - val_ra) < tolerance, 0, x_wrapped)

        #Initial mask: points within patch
        combined_mask = patch_path.contains_points(np.vstack((x_wrapped, y)).T)

        #If patch is near North Pole
        if np.any(near_north_pole):
            #Reflect RA points for dec > 90
            mirrored_ra_vertices = (val_ra - ra_vertices[near_north_pole]) % val_ra
            mirrored_dec_vertices = 180 - dec_vertices[near_north_pole]  #Reflect dec across 90

            mirrored_patch_vertices = np.column_stack((mirrored_ra_vertices, mirrored_dec_vertices))
            mirrored_patch_path = Path(mirrored_patch_vertices)

            #Apply mask for points within mirrored region
            mirrored_mask = mirrored_patch_path.contains_points(np.vstack((x_wrapped, y)).T)
            combined_mask = np.logical_or(combined_mask, mirrored_mask)

        #If patch is near South Pole
        if np.any(near_south_pole):
            #Reflect RA points for dec < -90
            mirrored_ra_vertices = (val_ra - ra_vertices[near_south_pole]) % val_ra
            mirrored_dec_vertices = -180 - dec_vertices[near_south_pole]  #Reflect dec across -90

            mirrored_patch_vertices = np.column_stack((mirrored_ra_vertices, mirrored_dec_vertices))
            mirrored_patch_path = Path(mirrored_patch_vertices)

            #Apply mask for points within mirrored region
            mirrored_mask = mirrored_patch_path.contains_points(np.vstack((x_wrapped, y)).T)
            combined_mask = np.logical_or(combined_mask, mirrored_mask)

        return combined_mask if In else np.logical_not(combined_mask)

    #Apply check for points inside patch and with mirrored RA
    inside_mask = _points_in_or_out_with_mirrored_RA(
        p.ra.deg[mask],
        p.dec.deg[mask],
        patchname
    )

    #Select points inside region
    p_inside = p[mask][inside_mask]

    #RM values and errors of those points
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
    
    #Create a circle (patch) for defining region
    circle = Circle((ra, dec), radius=radius, edgecolor='k', fill=False, linestyle='-', 
                    linewidth=1.1, transform=ax.get_transform('world'))
    
    if plot[0]: ax.add_patch(circle)  #To plot circle
    
    if plot[1]:  #To plot scatter points in red and blue depending on pos or neg
        alpha = 1
        #Getting positive positions and their RM values inside region
        pos_in_region, rm_value_pos, _ = points_inOrOut_patch(positions, circle, pos_mask, rm_s, rm_errs, typ='circle')
        ax.scatter(pos_in_region.ra, pos_in_region.dec,
                   transform=ax.get_transform(transform), marker='o', 
                   s=1, color=color_pos, alpha=alpha)
        
        #Getting negative positions and their RM values inside region
        neg_in_region, rm_value_neg, _ = points_inOrOut_patch(positions, circle, neg_mask, rm_s, rm_errs, typ='circle')
        ax.scatter(neg_in_region.ra, neg_in_region.dec,
                   transform=ax.get_transform(transform), marker='o', 
                   s=1, color=color_neg, alpha=alpha)
    
    #Get positive positions and their RM values inside region without plotting
    pos_in_region, rm_value_pos, rm_err_pos = points_inOrOut_patch(positions, circle, pos_mask, rm_s, rm_errs, typ='circle')
    #Get negative positions and their RM values inside region without plotting
    neg_in_region, rm_value_neg, rm_err_neg = points_inOrOut_patch(positions, circle, neg_mask, rm_s, rm_errs, typ='circle')
    
    if return_data:
        if return_err:  #Return errors (default)
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
    
    #Createrectangle (patch) for defining region
    #Note that position of rectangel is at bottom corner
    #But box should be shifted so that point chosen is same as center of box
    bottom_left_x = ra - width / 2
    bottom_left_y = dec - height / 2
    rect = Rectangle((bottom_left_x, bottom_left_y), width=width, height=height,
                     angle=angle + polar_angle, transform=ax.get_transform('world'))

    if plot[1]: #to plot scatter in red and blue depending on pos or neg.
        #Plottingrectangle
        alpha = 1
        
        marker_size = kw["marker_size"] if "marker_size" in kw else 1
        
        #Getting positive positions and their RM values inside region
        pos_in_region, rm_value_pos, _ = points_inOrOut_patch(positions, rect, pos_mask, rm_s, rm_errs)
        
        #Getting negative positions and their RM values inside region
        neg_in_region, rm_value_neg, _ = points_inOrOut_patch(positions, rect, neg_mask, rm_s, rm_errs)
       
        vmin, vmax = -50, 50  #Maximum and minimum RM limits
        #Plotting RM positions inside region
        sctt = ax.scatter(np.concatenate([neg_in_region.ra.deg, pos_in_region.ra.deg]), 
                          np.concatenate([neg_in_region.dec.deg, pos_in_region.dec.deg]), 
                          transform=ax.get_transform(transform), marker='o',
                          cmap='brg_r', vmin=vmin, vmax=vmax,
                          c=np.concatenate([rm_value_neg, rm_value_pos]))  #Color by RM values
        
        rotation = -90  #Rotation for colorbar labels
        cbars_pad = 0.01  #Padding for colorbars
        labelpad = 20  #Padding for colorbar labels
        cbar_lab_size = 17  #Font size for cbar labels
        l_s = 14  #label size for ticklabels (numerical)

        #Colorbar for scatter plot (for all positive or negative points)
        sctt_cbar = plt.colorbar(sctt, ax=ax, fraction=0.06, pad=cbars_pad)
        sctt_cbar.set_label("RM scatter (capped)", rotation=rotation,
                            labelpad=labelpad, fontsize=cbar_lab_size)
        ticks = np.linspace(vmin, vmax, 11)
        ticklabels = [str(round(tick)) for tick in ticks]
        sctt_cbar.set_ticks(ticks)
        sctt_cbar.set_ticklabels(ticklabels)
        sctt_cbar.ax.tick_params(labelsize=l_s)
    
    #Getting positive positions and their RM values inside region without plotting
    pos_in_region, rm_value_pos, rm_err_pos = points_inOrOut_patch(positions, rect, pos_mask, rm_s, rm_errs)
    #Getting negative positions and their RM values inside region without plotting
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

    ##Calculating residuals and reduced chi-square
    #residuals = Y - func(X, *popt)
    #chi_squared = np.sum(residuals ** 2)
    #dof = len(Y) - len(popt)  #Degrees of freedom
    #reduced_chi_square = chi_squared / dof
    #print(f"Reduced Chi-Square: {reduced_chi_square}")
    
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
    
    #Extracting HI density data
    HI_dens_data = hdu.data[0, 0, :, :]
    
    #Extract RM data positions and magnitudes
    RM_RA, RM_Dec = ra, dec
    
    #Converting RM positions to pixel coordinates
    RM_pix_x, RM_pix_y = ra_dec_to_pixels(RM_RA, RM_Dec, fitfile=filename)
    
    RM_pix_x = np.round(RM_pix_x).astype(int)
    RM_pix_y = np.round(RM_pix_y).astype(int)
    
    #Getting shape of HI density data array
    shape = HI_dens_data.shape
    
    #Finding minimum non-NaN HI density value
    min_HI_dens_value = np.nanmin(HI_dens_data)
    
    #Creating an empty array to store HI_dens_values
    HI_dens_values = np.full_like(rm_values, min_HI_dens_value, dtype=float)
    
    #Mask for valid coordinates of RM. Valid meaning within BT04 data
    valid_mask = (RM_pix_x >= 0) & (RM_pix_x < shape[1]) & (RM_pix_y >= 0) & (RM_pix_y < shape[0])

    ##This calculation shows 99% of RM points were rettained.    
    #print(sum(valid_mask)/len(valid_mask)); import sys; sys.exit() 
    
    #Getting valid coordinates that only within BT04 fit file
    valid_x = RM_pix_x[valid_mask]
    valid_y = RM_pix_y[valid_mask]
    
    #Fetching HI density values at those valid RM positions
    valid_HI_dens_values = HI_dens_data[valid_y, valid_x]
    
    #Handling RM Positions with NaN column density are overwritten to have minimum value possible based on fit file's minimum column density
    valid_HI_dens_values = np.where(np.isnan(valid_HI_dens_values), min_HI_dens_value, valid_HI_dens_values)
    
    #Assign back to HI_dens_values array with overwritten valid_HI_dens_values
    HI_dens_values[valid_mask] = valid_HI_dens_values
    
    result_tuples = list(zip(rm_values, HI_dens_values))
    
    #Extract rm and HI_dens_values from result_tuples
    rm_values, HI_dens_values = zip(*result_tuples)
    
    ##this shows that 87% of RM within BT04 data was flagged as Nan.
    #count = np.sum(np.array(HI_dens_values)< 14)
    #print(100*(1 - count/sum(valid_mask))); import sys; sys.exit()
    
    #print(min(HI_dens_values), max(HI_dens_values)); import sys; sys.exit()
    fmt = "bs" if above_zero else "rs"
    plt.errorbar(HI_dens_values, rm_values, yerr=rm_errors, fmt=fmt, 
                 markersize=3, capsize=3, ecolor='k')
    
    A = 35
    f = 0.206
    phi = 29.5
    
    n = np.linspace(min(HI_dens_values), max(HI_dens_values), int(1e5))
    
    if above_zero in [True, None]:
        plt.plot(n, -A * np.sin(2 * np.pi * f * n - phi) + 65)#,
                 #label=r"$(-A) \cdot \sin\left(2\pi f \cdot x - \phi\right) +65$")
    if above_zero in [False, None]:
        plt.plot(n, A * np.sin(2 * np.pi * f * n - phi) - 104)#,
                 #label=r"$A \cdot \sin\left(2\pi f \cdot x - \phi\right) -104$")
    
    #string = f'A = {A}, f = {f}, ' + r"$\phi=$ " + f"${phi}$"
    #plt.text(0.6, 1.03, string, transform=plt.gca().transAxes, fontsize=14,
    #         bbox=dict(facecolor='none', edgecolor='none', boxstyle='round,pad=0'))

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
    #OpeningFITS file and extracting data
    filename = 'LGSNLGSR.SQLGBB.FITS'
    hdu = fits.open(filename)[0]
    
    #Setting upWCS projection usingheader
    wcs = WCS(hdu.header)
    
    positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0
    
    #Creating figure and axis with WCS projection
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

def curr_dir_path():
    """Returns folder path of currently running script."""
    return os.path.dirname(os.path.abspath(__file__)) + "/"

def get_distinct_colors(n, cmap_name='tab20b'):
    #Was meant to be used for RM vs PA plot
    """
    Returns `n` visually distinct colors using a qualitative matplotlib colormap.

    If n > base colormap size, it interpolates smoothly across full range of colormap.
    """
    cmap = plt.get_cmap(cmap_name)

    if hasattr(cmap, "colors"):  #For tab10, tab20, etc.
        base_colors = cmap.colors
        base_n = len(base_colors)

        if n <= base_n:
            return base_colors[:n]
        else:
            indices = np.linspace(0, base_n - 1, n)
            colors = [base_colors[int(i)] for i in indices]
            return colors
    else:
        #Continuous colormap fallback (like viridis)
        return [cmap(i / (n - 1)) for i in range(n)]

def create_annuli_binning_structure(bin_data, data, bin_num, data_errs=None, for_PA_or_b_plot=False):
    """
    Bins x data, assigns y values/errors to bins, 
    and returns structured dictionaries. 
    Very useful for Annulus analysis.
    """
    
    #bin_data = np.asarray(bin_data) #making sure its numpy

    #(bin_data=m31_distances, 
    #data=(rm_pos_gal_lat, rm_m31), 
    #bin_num=bin_num_from_main+1, 
    #for_PA_or_b_plot=True)

    assert len(bin_data) == len(data[0]) == len(data[1]), "bins are not same length in create_annuli_binning_structure()"

    #Binning edges and indices for annuli
    bin_edges = np.linspace(min(bin_data), max(bin_data), bin_num)
    bin_indices = np.digitize(bin_data, bins=bin_edges)  #Bins indexed from 1 to bin_num-1

    #Compute bin areas
    annul_area = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)  #Now a NumPy array
    if not for_PA_or_b_plot: bin_edges = bin_edges[:-1]  #Maintaining same number of elements as area

    #Initialize storage dictionaries
    Dicts = []
    for j in range(len(data)):  #In case it's for x and y axis for non-histogram/bar plots
        y_per_bin = {i: [] for i in range(0, bin_num+1)}  #Bin indices go from 1 to bin_num (instead of bin_num-1)
        if data_errs is not None:
            y_err_per_bin = {i: [] for i in range(0, bin_num+1)}
            Dicts.append([y_per_bin, y_err_per_bin])
        else:
            Dicts.append([y_per_bin,])

    #Assigning values to bins
    for i, bin_idx in enumerate(bin_indices):
        if bin_idx < bin_num:  #Ensure no out-of-range errors
            for idx, dicts in enumerate(Dicts):
                dicts[0][bin_idx].append(data[idx][i])
                if data_errs is not None:
                    dicts[1][bin_idx].append(data_errs[idx][i])

    #Convert lists to NumPy arrays (filtering out empty bins)
    for dicts in Dicts:
        dicts[0] = {k: np.array([val.value if hasattr(val, "unit") else val for val in v]) for k, v in dicts[0].items() if len(v) > 0}
        if len(dicts) > 1:
            dicts[1] = {k: np.array([val.value if hasattr(val, "unit") else val for val in dicts[1][k]]) for k in dicts[0]}

    if data_errs is not None:
        return [dicts for dicts in Dicts], annul_area, bin_edges
    return [dicts[0] for dicts in Dicts], annul_area, bin_edges

def apply_plot_attributes(push_title_up = 1.3, leg=True, xlim=(0,180),**kw):
    """
    Primarily used in plot_binned_PA() 
    but can also be used for plot_binned_gal_lat() amongst others...
    """
    plt.xlabel(kw.get("xlabel",
                      
                      #Default is meant for plotting binned PA
                      "Polar Angle " r"$(\theta)$ [$^{\circ}$] " + "(Anticlockwise from North - 37.7" + r"$^{\circ}$)"
                      ),
                      fontsize = kw["axis_lab_f_size"][0] if "axis_lab_f_size" in kw else 12)
    plt.ylabel("Rotation Measure " + r"[rad m$^{-2}$]",
                      fontsize = kw["axis_lab_f_size"][1] if "axis_lab_f_size" in kw else 12)
    if leg:
        plt.legend(fontsize = kw.get("fontsize", 9), 
                   loc = 'upper center', 
                   bbox_to_anchor = kw.get("bbox",(1.13, push_title_up)),
                    framealpha = kw.get("framealpha",0), 
                    ncols = kw.get("ncols", 6))
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.minorticks_on() 
    if xlim is not None: plt.xlim(xlim)
    if "ylim" in kw: 
        if kw["ylim"] is not None: plt.ylim(*kw["ylim"])

def get_discrete_colors(data_limits, n_bins, cmap_name):
    
    radial_lims = data_limits
    bin_edges = np.linspace(*radial_lims, n_bins)
    bin_width = bin_edges[1] - bin_edges[0]
    #bin_centers = bin_edges[:-1] + bin_width / 2  #Bin centers

    #Creating a discrete version of colormap
    base_cmap = plt.get_cmap(cmap_name)
    cmap = base_cmap(np.linspace(0, 1, n_bins))  #n_bins distinct colors
    return plt.matplotlib.colors.ListedColormap(cmap) #referencing this with variable "cmap_discrete"

