#NOTE: Increase/Decrease the value of 'sigma_detect_limit' in order to change the smoothing results
#BE CAREFUL not to oversmooth.
from personal_converter import convert_magellanic_stream_coords
from astropy.coordinates import Angle, FK4, ICRS
from astropy.time import Time
from reproject import reproject_interp

from main import (
    
# Imprting alias'
fits, WCS, plt, u, np, Table, SkyCoord,

#importing variables
m31_pos, Circle, Rectangle, rm_pos_icrs, 
rm_m31, eq_pos, m31_condition, m33_pos,
rm, bg_condition,
bg_pos_icrs,
args,


# Importing functions
smooth_2d_image, 
convert_txt_to_Skycoord, 
ra_dec_to_pixels,
curr_dir_path,
get_data_from_catalogue,
get_CGM_and_BG_masks,
apply_CGM_and_BG_masks,
get_discrete_colors
)

# # 1)Loading & transforming catalogue, get all raw coords + RMs + M31/M33 stuff
# (RM_lat, RM_lon, rm, rm_err, position, eq_pos, rm_m31_coord, 
#  m31_sep, m31_theta, new_frame_of_reference,
#  cloud6_pos, m33_pos, m33_m31coord, m33_sep, m33_theta
# ) = get_data_from_catalogue(sigma_detect_limit=args.sig_limit)

# # 2)Building CGM / BG masks (ellipse vs circle controlled by flags)
# m31_condition, bg_condition = get_CGM_and_BG_masks(
#     rm_m31_coord, eq_pos, m31_sep,
#     elliptic_CGM=args.elliptic_CGM,
#     elliptic_CGM_bg=args.elliptic_CGM_bg,
#     # cutoff=(N+10)*u.deg,
#     # L_m31=N #degrees
# )

# # 3)Apply those masks to slice out CGM & BG subsets
# (bg_pos, bg_pos_icrs, rm_pos, rm_pos_icrs, rm_pos_gal_lat, 
# rm_pos_gal_lat_bg, rm_bg, m31_sep_bg, err_bg, rm_m31, 
# m31_sep_Rvir, err_m31) = apply_CGM_and_BG_masks(
#     rm_m31_coord, eq_pos, position, rm, rm_err, 
#     m31_sep, m31_condition, bg_condition
# )

levels = 10 #For HI density contour
levels_2 = 10 #For HI density contour
cmap = 'jet'  # For HI density
rotation = -90  # Rotation for colorbar labels
cbars_pad = 0.01  # Padding for colorbars
labelpad = 20  # Padding for colorbar labels
cbar_lab_size = 13  # Font size for cbar labels
l_s = 14  # label size for ticklabels (numerical)

#Opening FITS file and extracting data
filename = 'LGSNLGSR.SQLGBB.FITS'
hdu = fits.open(filename)[0]

#Adding in new fit file
path = curr_dir_path()
filename_2 = f'{path}skv11594169016096.fits'
hdu2 = fits.open(filename_2)[0]

#Setting up WCS projection usingheader
wcs = WCS(hdu.header)

#Creating figure and axis with WCS projection
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y', 0, 0))

#2D WCS (RA-Dec only)
wcs1 = WCS(hdu.header).celestial
wcs2 = WCS(hdu2.header).celestial

target_shape = hdu.data.shape[-2:]  # Target shape from original image --> (326, 478)

#Reprojecting second image
reproj_data, _ = reproject_interp((hdu2.data, wcs2), wcs1, shape_out=target_shape)

# Overlay contours form new fit file
if args.filled_485:
    skv_cont = ax.contourf(reproj_data, levels=levels_2, cmap='viridis')
    skv_cbar = plt.colorbar(skv_cont, ax=ax, fraction=0.05, orientation='vertical')

    skv_cbar.set_label(f'Temperature [{hdu2.header["BUNIT"]}] ({hdu2.header["SURVEY"]})', 
                    labelpad=labelpad, fontsize=cbar_lab_size, rotation=rotation)

if args._485mhz:
    ax.contour(reproj_data, levels=levels_2, 
            colors='green' if args.HI else 'k', 
            linewidths=1.3)
    
print(f"{len(rm_m31)=}")
# # Gets the levels of contours data (The actual values per contour)
# # Get min and max from your data
# vmin = np.nanmin(reproj_data)
# vmax = np.nanmax(reproj_data)
# # Generate 16 contour levels
# levels_2 = np.linspace(vmin, vmax, levels_2)
# print(f"Contour temperature levels (K): {np.round(levels_2)}")

hdu_data = hdu.data[0, 0, :, :].copy()
hdu_data[hdu_data < 17] = np.nan

## Adding contours for HI data (Braun & Thilker 2004)
if args.HI:
    ax.contour(hdu_data, levels=10,
            linewidths=1, colors='black',
            linestyles='solid')

# # Gets the levels of contours data (The actual values per contour)
# # Get min and max from your hdu_data
# vmin = np.nanmin(hdu_data)  # min ignoring NaNs
# vmax = np.nanmax(hdu_data)  # max ignoring NaNs

# # Generate 16 contour levels
# levels_2 = np.linspace(vmin, vmax, level)

# Print contour temperature levels
# print(f"Contour HI levels: {np.round(levels_2,1)}"); import sys; sys.exit()

def add_resolution_to_plot(ax, pixel_rad):
    #Shifiting Parametre of encapualtion for FWHM of 49 arcsecs (or for the new 0.85 degrees)
    h = 21.5 #In deegrees

    *center_pixel, _,_ = wcs.world_to_pixel_values(m31_pos.ra.deg-h, #Move to the right of m31
                                            m31_pos.dec.deg-h*1.15, #Move to downwards from m31
                                            0 ,0) #Dealing with 4D data  
    ang_rsltn = Circle(center_pixel, 
                        radius = pixel_rad,
                        color = 'black', fill=False,
                        linestyle = '-', linewidth=2)
    ax.add_patch(ang_rsltn)#Adding circle to axis

    #Calculating lower-left corner of the square that will encapsulate ang_rsltn
    w=3 #Width paramtre to increase or decrease square that encapsulates ang_rsltn
    lower_left = (center_pixel[0] - pixel_rad*w/2, center_pixel[1] - pixel_rad*w/2)

    #Create the square patch to encapsulate the angular resoolutionn of 49 arcminutes
    square_patch = Rectangle(lower_left,
                        width=w * pixel_rad,
                        height=w * pixel_rad,
                        edgecolor='black', fill=False, 
                        linestyle='-', linewidth=2)
    ax.add_patch(square_patch)#Adding square to axis

#Creating a perfect circle to show angular resolution of BT04 data
pixel_rad_1 = (49*u.arcmin).to(u.deg).value / np.abs(hdu.header['CDELT1']) #dividing by degree per pixel scale

#For skv.....fits files header -> (COMMENT Resolution:  0.85 degrees)
# (COMMENT PixelScale:  0.3515 degrees/pixel )
pixel_rad_2 = 0.85 / np.abs(hdu.header['CDELT1'])

#Deciding which to plot
if args.HI: 
    add_resolution_to_plot(ax=ax, pixel_rad=pixel_rad_2)
else:
    add_resolution_to_plot(ax=ax, pixel_rad=pixel_rad_1)
    


ax.grid(ls='solid', alpha=.4, color="black")
f_s, tick_f_s = 15, 14
ax.coords[0].set_axislabel('RA [J2000]', fontsize=f_s)
ax.coords[0].set_ticklabel(fontsize=tick_f_s)
ax.coords[1].set_axislabel('Dec [J2000]', fontsize=f_s)
ax.coords[1].set_ticklabel(fontsize=tick_f_s)

vmin, vmax = -300, 300 #max(rm_m31) # Maximum and minimum RM limits
rm_m31_clipped = np.clip(rm_m31, vmin, vmax)

if args.scatter:
    # Marker in order to show vertical colorbar
    sctt = ax.scatter(rm_pos_icrs.ra, rm_pos_icrs.dec,
                    c=rm_m31, marker='o', s=20,
                    cmap='RdYlBu', vmin=vmin, vmax=vmax,
                    transform=ax.get_transform('world'))
    cbar = plt.colorbar(sctt, ax=ax, orientation='vertical')
    cbar.set_label(label='RM Scatter [rad/m²]', rotation=-90, labelpad=20)

# Create a scatter plot for positive and negative
positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0

ra, dec = eq_pos.ra[m31_condition], eq_pos.dec[m31_condition]
ra_bg, dec_bg = eq_pos.ra[bg_condition], eq_pos.dec[bg_condition]
# ra = np.concatenate((pos_in_region.ra, neg_in_region.ra))
# dec = np.concatenate((pos_in_region.dec, neg_in_region.dec))

# global im
# nsig is a multiplier for sigma (std) in y and x direction.
im, kernel = smooth_2d_image(ra, dec, imsize=450, kernel=2, fitfile=filename, rm_m31=rm_m31_clipped)
im_clipped = im #np.clip(im, -5000, 5000) #manually clipping based off of what i can see on defualt colorbar
# x0s_pos, y0s_pos = ra_dec_to_pixels(ra[positive_mask], dec[positive_mask], filename=filename)
# x0s_neg, y0s_neg = ra_dec_to_pixels(ra[negative_mask], dec[negative_mask], filename=filename)
if args.bg: x0s_bg, y0s_bg = ra_dec_to_pixels(ra_bg, dec_bg, filename=filename)
# #image 'im' is already trnsformed intowcs.ax's 'world'
# ax.scatter(x0s_pos, y0s_pos, marker='o', color="b", s=2)
# ax.scatter(x0s_neg, y0s_neg, marker='o', color="r", s=2)
if args.bg: ax.scatter(x0s_bg, y0s_bg, marker='x', s=.05, color='k', alpha=.6)

if args.smoothed:
    # plt.title(f"rad-limit = {N} deg")
    n_bins=10
    # Get discrete colormap and bin edges
    cmap_discrete, bin_edges = get_discrete_colors((vmin, vmax), n_bins, 'jet_r', get_edges=True)
    norm_cb = plt.matplotlib.colors.BoundaryNorm(boundaries=bin_edges, ncolors=n_bins)

    # Plot imshow
    # RM_dens = ax.imshow(im, cmap='jet_r', origin='lower', alpha=0.5, zorder=-1, vmin=vmin, vmax=vmax)
    RM_dens = ax.imshow(im, cmap=cmap_discrete, norm=norm_cb, origin='lower', alpha=0.5, zorder=-1)

    # Colorbar
    RM_dens_cbar = plt.colorbar(RM_dens, ax=ax, fraction=0.05, orientation='horizontal')
    RM_dens_cbar.set_label('Smooth RM Intensity (capped)', labelpad=labelpad, fontsize=cbar_lab_size)

    # Optional: place ticks at bin centers
    # bin_centers = bin_edges[:-1] + (bin_edges[1] - bin_edges[0]) / 2
    # RM_dens_cbar.set_ticks(bin_centers)

    # RM_dens_cbar = plt.colorbar(RM_dens, ax=ax, fraction=0.05, orientation='horizontal')
    # RM_dens_cbar.set_label('Smooth RM Intensity (capped)', labelpad=labelpad, fontsize=cbar_lab_size)

    # steps = 11
    # min_tick_1 = np.min(rm_m31_clipped) ; max_tick_1 = np.max(rm_m31_clipped) 
    # min_tick_2 = np.min(im_clipped) ;  max_tick_2 = np.max(im_clipped)
    # m = (max_tick_2 - min_tick_2) / (max_tick_1 - min_tick_1)
    # RM_dens_cbar.set_ticks(np.linspace(min_tick_2, max_tick_2, steps))
    # ticklabs = np.round(np.linspace(min_tick_2/m, max_tick_2/m,steps))
    # # ticklabs = np.round(np.linspace(min_tick_1, max_tick_1,steps))
    # RM_dens_cbar.set_ticklabels(ticklabs.astype(int))
    # RM_dens_cbar.ax.tick_params(labelsize=l_s)

#Marking M31 and M33 position
color = 'white' if args.filled_485 else "k"
sctt_m31 = ax.scatter(m31_pos.ra, m31_pos.dec,color=color, marker='*',
        transform=ax.get_transform('world'), label = "M31", s=30)
sctt_m33 = ax.scatter(m33_pos.ra, m33_pos.dec,color=color, marker='^',  
        transform=ax.get_transform('world'), label = "M33", s=20)

more=0 #Parametre to move label legends up a bit

def get_dwarf_G_info(more):
    # ICRS Coordinates of dwarf galaxies near M31
    # from paper doi:10.3847/0004-637X/816/2/81
    # with many more additons from paper: https://iopscience.iop.org/article/10.1088/0004-6256/144/1/4/pdf
    # Existing list of dwarf galaxy coordinates near M31
    global dwf_Gs
    # dwf_Gs = convert_txt_to_Skycoord("dwarf_galaxies.txt", withNames=True)
    # ax.scatter(dwf_Gs[0].ra, dwf_Gs[0].dec, 
    dwarf_galaxies = Table.read("dwarf_galaxies.dat", format='ascii')
    dwarf_galaxies.keep_columns(["l_MS", "b_MS", "Name"])
    dwarf_galaxies["l_MS"].unit, dwarf_galaxies["b_MS"].unit = u.deg, u.deg #Giving
    dwf_Gs = convert_magellanic_stream_coords(dwarf_galaxies["l_MS"], dwarf_galaxies["b_MS"])

    # Plotting Dwarf Galaxies
    ax.scatter(dwf_Gs.ra, dwf_Gs.dec, 
            #Changing color due to mentioning in http://dx.doi.org/10.3847/0004-637X/816/2/81
            #That these satellite galaxies are what could have casued the fomartion of the filament region
            color=["red" if e in ["And II", "And XV"] else "orange" for e in dwarf_galaxies["Name"]], 
            marker='^',
        transform=ax.get_transform('world'), s=100, edgecolor='k',
        label="dwarf Galaxies"
        )
    more += 0.03
if args.dwarf_G:
    get_dwarf_G_info(more)

if args.others:
    get_dwarf_G_info(more)
    # ICRS Coordinates of discrete CHVC between M33 and M31
    # from paper doi:10.3847/0004-637X/816/2/81
    cloud_coords = convert_txt_to_Skycoord("clouds.txt")

    label_clouds=True #To add label in legend for regions of clouds 
    if label_clouds:
        ax.scatter(cloud_coords.ra, cloud_coords.dec, color="yellow", 
                marker='s',
        transform=ax.get_transform('world'), s=80, edgecolor='k', 
        label="discrete clouds")
        label_clouds=False #stop repetition
        more+=0.03
    else:
        ax.scatter(cloud_coords.ra, cloud_coords.dec, color="yellow", 
                marker='s',
        transform=ax.get_transform('world'), s=80, edgecolor='k')
            
    def extract_continuum_data(file, **kw):
        #Continuum Absoroption features form paper:
        #from paper https://ui.adsabs.harvard.edu/abs/1993ApJ...405..153D/abstract
        contSources = Table.read(file, format='ascii')
        contSources.keep_columns(['NAME','RA_1950','DEC_1950','N_H','T_B_max'])
        label = kw["label_legend"] if "label_legend" in kw else ""

        contSources["Coords"] = SkyCoord(ra=contSources["RA_1950"], 
                                        dec=contSources["DEC_1950"], 
                                        unit = (u.hourangle, u.deg), 
                                        frame=FK4(equinox="B1950.0"))

        contSources["RA"], contSources["DEC"] = contSources["Coords"].ra, contSources["Coords"].dec #incase i wanna deal with one or the other.
        contSources["N_H"] = np.int64(contSources["N_H"])
        contSources["T_B_max"] = np.float64(contSources["T_B_max"])
        contSources["N_H"].unit, contSources["T_B_max"].unit = (10**(19))/(u.cm**2), u.K

        cords=contSources["Coords"]
        temp = np.array(list(map(int,([0] + list(contSources["T_B_max"])[1:]))))
        ax.scatter(cords.ra, cords.dec, 
                color='k', marker='x',transform=ax.get_transform('world'),
                s=np.int64((temp/min(temp[temp!=0]))*17), label=label
            )

    extract_continuum_data("M31_contSources.txt", label_legend= "Continuum Sources")
    extract_continuum_data("M33_contSources.txt")

    # Coordniates of Backgorund AGN 
    #from https://doi.org/10.3847/1538-4357/aa87b4
    AGN_bg = convert_txt_to_Skycoord("AGN_background.txt")
    ax.scatter(AGN_bg.ra, AGN_bg.dec, color="blue", marker='*',
        transform=ax.get_transform('world'), s=100, edgecolor='k',
        label="Background AGN")

from astropy.visualization.wcsaxes import CoordinateHelper
ra_axis: CoordinateHelper = ax.coords[0]
ra_axis.set_ticks(spacing=10 * u.deg)  #RA tick spacing
ra_axis.set_major_formatter('d.d')  #Setting format of RA axis to degrees

#ZOOMING OUT (in degrees)
ra_min, ra_max = 40, -30
dec_min, dec_max = 16, 48

# #ZOOMING IN (in degrees)
# ra_min, ra_max = 30, -4
# dec_min, dec_max = 23, 50

# Convert these to pixel coordinates using WCS
bottom_left = wcs.world_to_pixel_values(ra_min, dec_min, 0, 0)
top_right = wcs.world_to_pixel_values(ra_max, dec_max, 0, 0)

#Extracting pixel coordinates
x_min, y_min, _, _ = bottom_left
x_max, y_max, _, _ = top_right

#Setting axis limits
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

plt.legend(fontsize = 15, loc = 'upper center', bbox_to_anchor = (0.5, 1.16+more), 
            framealpha = 1 if args.filled_485 else 0, 
            ncols = 3, markerscale=1.1, facecolor="grey" if args.filled_485 else "None")
#            handles=[sctt_m31, sctt_m33, sctt_bg])

plt.tight_layout()
if args.save_plot:
    path = curr_dir_path() + "Results/"
    if args.smoothed:
        filename = f"2D_Gaussian_RM_Smoothing_kernel={tuple(np.round(kernel, 0))}_in_degrees.png"
        plt.savefig(f"{path}{filename}", dpi=600, bbox_inches="tight")
        print(f"2D_Gaussian_RM_Smoothing has been saved to {path}{filename}")
    else:
        filename = f"Temperature_Contours_within_CGM_of_M31.png"
        plt.savefig(f"{path}{filename}", dpi=600, bbox_inches="tight")
        print(f"=Temperature Contour has been saved to {path}{filename}")
    plt.close()
else:
    plt.show()
