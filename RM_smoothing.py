#NOTE: Increase/Decrease the value of 'sigma_detect_limit' in order to change the smoothing results
#BE CAREFUL not to oversmooth.
from personal_converter import convert_magellanic_stream_coords
from astropy.coordinates import Angle, FK4, ICRS
from astropy.time import Time

from main import (
    
# Imprting alias'
fits, WCS, plt, u, np, Table, SkyCoord,

#importing variables
m31_pos, Circle, Rectangle, rm_pos_icrs, 
rm_m31, eq_pos, m31_condition, m33_pos,


# Importing functions
smooth_2d_image, 
convert_txt_to_Skycoord, 
ra_dec_to_pixels
)

levels = 16 #For HI density contour
cmap = 'jet'  # For HI density
rotation = -90  # Rotation for colorbar labels
cbars_pad = 0.01  # Padding for colorbars
labelpad = 20  # Padding for colorbar labels
cbar_lab_size = 17  # Font size for cbar labels
l_s = 14  # label size for ticklabels (numerical)

#Opening FITS file and extracting data
filename = 'LGSNLGSR.SQLGBB.FITS'
hdu = fits.open(filename)[0]

#Setting up WCS projection usingheader
wcs = WCS(hdu.header)

#Creating figure and axis with WCS projection
fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y', 0, 0))

#Adding contours
ax.contour(hdu.data[0, 0, :, :], levels=levels,
           linewidths=1, colors='black',
           linestyles='solid')

#Creating a perfect circle to show angular resolution of BT04 data
pixel_rad = (49*u.arcmin).to(u.deg).value / np.abs(hdu.header['CDELT1']) #dividing by degree per pixel scale

#Shifiting Parametre of encapualtion for FWHM of 49 arcsecs
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

ax.grid(ls='solid', alpha=.4, color="black")
f_s, tick_f_s = 24, 14
ax.coords[0].set_axislabel('RA [J2000]', fontsize=f_s)
ax.coords[0].set_ticklabel(fontsize=tick_f_s)
ax.coords[1].set_axislabel('Dec [J2000]', fontsize=f_s)
ax.coords[1].set_ticklabel(fontsize=tick_f_s)

vmin, vmax = -50, 50  # Maximum and minimum RM limits

# Marker in order to show vertical colorbar
sctt = ax.scatter(rm_pos_icrs.ra, rm_pos_icrs.dec,
                  c=rm_m31, marker='', s=20,
                  cmap='brg_r', vmin=vmin, vmax=vmax,
                  transform=ax.get_transform('world'))

# Create a scatter plot for positive values
positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0

ra, dec = eq_pos.ra[m31_condition], eq_pos.dec[m31_condition]
# ra = np.concatenate((pos_in_region.ra, neg_in_region.ra))
# dec = np.concatenate((pos_in_region.dec, neg_in_region.dec))

# global im
# nsig is a multiplier for sigma (std) in y and x direction.
im = smooth_2d_image(ra, dec, imsize=450, nsig = .6, fitfile=filename)
im_clipped = np.clip(im, -5000, 5000) #manually clipping based off of what i can see on defualt colorbar

x0s, y0s = ra_dec_to_pixels(ra, dec, filename=filename)
#image 'im' is already trnsformed intowcs.ax's 'world'
# ax.scatter(x0s, y0s, marker='.', s=2)

RM_dens = ax.imshow(im_clipped, cmap='RdYlBu', origin='lower', alpha=0.5)
   
RM_dens_cbar = plt.colorbar(RM_dens, ax=ax, fraction=0.05, orientation='horizontal')
RM_dens_cbar.set_label('Smooth RM Intensity (capped)', labelpad=labelpad, fontsize=cbar_lab_size)
global min_tick_1, max_tick_1, min_tick_2, max_tick_2

steps = 11
rm_m31_clipped = np.clip(rm_m31, vmin, vmax)
min_tick_1 = np.min(rm_m31_clipped) ; max_tick_1 = np.max(rm_m31_clipped) 
min_tick_2 = np.min(im_clipped) ;  max_tick_2 = np.max(im_clipped)
RM_dens_cbar.set_ticks(np.linspace(min_tick_2, max_tick_2, steps))
ticklabs = np.round(np.linspace(min_tick_1, max_tick_1,steps))
RM_dens_cbar.set_ticklabels(ticklabs.astype(int))
RM_dens_cbar.ax.tick_params(labelsize=l_s)

#Marking M31 and M33 position
color = 'black'
sctt_m31 = ax.scatter(m31_pos.ra, m31_pos.dec,color=color, marker='*',
        transform=ax.get_transform('world'), label = "M31", s=50)
sctt_m33 = ax.scatter(m33_pos.ra, m33_pos.dec,color=color, marker='^',  
        transform=ax.get_transform('world'), label = "M33", s=50)

# ICRS Coordinates of discrete CHVC between M33 and M31
# from paper doi:10.3847/0004-637X/816/2/81
cloud_coords = convert_txt_to_Skycoord("clouds.txt")

label_clouds=True #To add label in legend for regions of clouds 
more=0 #Parametre to move label legends up a bit
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

# Plotting Dwarf Gaslaxies
ax.scatter(dwf_Gs.ra, dwf_Gs.dec, 
           #Changing color due to mentioning in http://dx.doi.org/10.3847/0004-637X/816/2/81
           #That these satellite galaxies are what could have casued the fomartion of the filament region
           color=["red" if e in ["And II", "And XV"] else "orange" for e in dwarf_galaxies["Name"]], 
           marker='^',
    transform=ax.get_transform('world'), s=100, edgecolor='k',
    label="dwrf G"
    )
more += 0.03

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

extract_continuum_data("M31_contSources.txt", label_legend= "Cont. Sources")
extract_continuum_data("M33_contSources.txt")

# # Coordniates of Backgorund AGN 
# #from https://doi.org/10.3847/1538-4357/aa87b4
# AGN_bg = convert_txt_to_Skycoord("AGN_background.txt")
# ax.scatter(AGN_bg.ra, AGN_bg.dec, color="blue", marker='*',
#     transform=ax.get_transform('world'), s=100, edgecolor='k',
#     label="Background AGN")

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

plt.legend(fontsize = 20, loc = 'upper center', bbox_to_anchor = (0.5, 1.12+more), 
            framealpha = 0, ncols = 3, markerscale=1.5)
#            handles=[sctt_m31, sctt_m33, sctt_bg])
# frame = legend.get_frame()
# frame.set_facecolor('gray')

plt.tight_layout()
print("week_2_plots part1 complete")
plt.show()





