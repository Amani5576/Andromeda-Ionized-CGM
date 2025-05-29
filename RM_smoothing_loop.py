#NOTE: Increase/Decrease the value of 'sigma_detect_limit' in order to change the smoothing results
#BE CAREFUL not to oversmooth.
from personal_converter import convert_magellanic_stream_coords
from astropy.coordinates import Angle, FK4, ICRS
from astropy.time import Time
from reproject import reproject_interp
from astropy.visualization.wcsaxes import CoordinateHelper
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.stats import sigma_clip
mpl.rcParams['font.family'] = 'Amiri'
from main import (
    
# Imprting alias'
fits, WCS, plt, u, np, Table, SkyCoord, ListedColormap,

#importing variables
m31_pos, Circle, Rectangle, args, L_m31,


# Importing functions
smooth_2d_image, 
convert_txt_to_Skycoord, 
ra_dec_to_pixels,
curr_dir_path,
get_data_from_catalogue,
get_CGM_and_BG_masks,
apply_CGM_and_BG_masks,
get_discrete_colors,
BG_correction,
galactic_cut,
simple_hist,

#Importing self-made classes
RMImageMasker
)


# 1)Loading & transforming catalogue, get all raw coords + RMs + M31/M33 stuff
(RM_lat, RM_lon, rm, rm_err, position, eq_pos, rm_m31_coord, 
m31_sep, m31_theta, new_frame_of_reference,
cloud6_pos, m33_pos, m33_m31coord, m33_sep, m33_theta
) = get_data_from_catalogue(sigma_detect_limit=args.sig_limit)

# 2)Building CGM / BG masks (ellipse vs circle controlled by flags)
m31_condition, bg_condition = get_CGM_and_BG_masks(
    rm_m31_coord, eq_pos, m31_sep,
    elliptic_CGM=args.elliptic_CGM,
    elliptic_CGM_bg=args.elliptic_CGM_bg,
    # cutoff=(N+10)*u.deg,
    # L_m31=N #degrees
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

rm_m31 = BG_correction(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg, upto_40=args.upto_40)#, All_rm=rm, mask_condition=m31_condition)
rm_bg = BG_correction(bg_pos_icrs, rm_bg, bg_pos_icrs, rm_bg, upto_40=args.upto_40)#, All_rm=rm, mask_condition=bg_condition)

rm_pos_gal_lat_old = rm_pos_gal_lat
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

rm_full_replaced = rm.copy()
# print(f"{len(rm_m31)=}")
# print(f"{len(rm_pos)=}")
if args.b_limit_small or args.b_limit:
    rm_full_replaced[m31_condition & (np.abs(position.b.value) > 5)] = rm_m31
    # print(f"{len(rm_full_replaced[m31_condition][b_lim_lat_mask])=}")
    # import sys; sys.exit()
else:
    print("Replacement succeeded 👍")
    rm_full_replaced[m31_condition] = rm_m31
    pass

# simple_hist(rm_m31, title="CGM")
# simple_hist(rm_full_replaced, title="Entire catalgoue")

levels = 10 #For HI density contour
levels_2 = 10 #For HI density contour
cmap = 'jet'  # For HI density
rotation = -90  # Rotation for colorbar labels
cbars_pad = -3  # Padding for colorbars
labelpad = 10  # Padding for colorbar labels
cbar_lab_size = 18  # Font size for cbar labels
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

cmocean_maps = ['thermal', 'haline', 'solar', 'ice', 'gray', 'oxy', 'deep', 
                'dense', 'algae', 'matter', 'turbid', 'speed', 'amp', 'tempo', 
                'rain', 'phase', 'topo', 'balance', 'delta', 'curl', 'diff', 'tarn']

for N in range(1,40):
    #Pick a radial limit of CGM just for visualization
    N = L_m31

    #Creating figure and axis with WCS projection
    fig = plt.figure(figsize=(12,12))
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
        cs = ax.contour(reproj_data, levels=levels_2, 
                colors='r' if args.HI else 'k', 
                linewidths=2.5, alpha = .4, label="485 MHz" )
        
        # #setting linestyle manually
        # for c in cs.collections:
        #     c.set_linestyle('-.')

        #Manually giving contour a label for the legend:
        _485mhz_cont = plt.Line2D([], [], color='r' if args.HI else 'k', 
                                  linewidth=2.5, label="485 MHz")
        
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
                linewidths=1, colors='k',
                linestyles='solid')
        
        #Manually giving contour a label for the legend:
        HI_cont = plt.Line2D([], [], color='k', linewidth=1, label="HI")

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
    f_s, tick_f_s = 18, 17
    ax.coords[0].set_axislabel('RA [J2000]', fontsize=f_s)
    ax.coords[0].set_ticklabel(fontsize=tick_f_s)
    ax.coords[1].set_axislabel('Dec [J2000]', fontsize=f_s)
    ax.coords[1].set_ticklabel(fontsize=tick_f_s)


    #2nd time around just to get the m31_condition for visualization
    m31_condition, bg_condition = get_CGM_and_BG_masks(
        rm_m31_coord, _, m31_sep,
        elliptic_CGM=args.elliptic_CGM,
        elliptic_CGM_bg=args.elliptic_CGM_bg,
        # cutoff=(N+10)*u.deg,
        L_m31=N #degrees
    )

    #Maximum and minimum RM limits
    vmin = np.floor(min(rm_full_replaced[m31_condition]))
    vmax = np.ceil(max(rm_full_replaced[m31_condition]))
    
    # rm_m31_clipped = np.clip(rm_m31, vmin, vmax)
    rm_m31_clipped = np.clip(rm_full_replaced, vmin, vmax)
    gauss_clipped_rm = sigma_clip(rm_full_replaced[m31_condition], 
                                  sigma=4, 
                                  maxiters=1e2)

    vmin = np.floor(min(gauss_clipped_rm))
    vmax = np.ceil(max(gauss_clipped_rm))

    ra, dec = eq_pos.ra[m31_condition], eq_pos.dec[m31_condition]

    # print(f"{sum(rm_m31_clipped)}")
    if args.scatter:

        vmin, vmax = np.round((min(rm_m31), max(rm_m31))) #max(rm_m31) # Maximum and minimum RM limits wihtin CGM of M31
        rm_m31_clipped = np.clip(rm_full_replaced, vmin, vmax)
        rm_m31_clipped = rm_m31_clipped[m31_condition]

        n_bins_sctt = 6 if args.percent_cbar else 15

        custom_colors = [
            "#cdba40",  # yellow
            "#995833",   # red
            "#7109C687",  # dark blue
            "#595959",  # black
            # "#dcbeff",  # lavender
            "#f58231",  # orange
            "#3cb44b",  # green
            "#4363d8",  # blue
        ]
        cmap = ListedColormap(custom_colors)

        # Get discrete colormap and bin edges
        cmap_discrete_sctt, bin_edges_sctt = get_discrete_colors((vmin, vmax), n_bins_sctt, cmap if args.percent_cbar else 'tab20', 
                                                       get_edges=True, percentiles=args.percent_cbar, data=rm_m31_clipped)

        if args.percent_cbar:
            norm_cb = plt.matplotlib.colors.BoundaryNorm(bin_edges_sctt, ncolors=cmap_discrete_sctt.N)
        else:
            norm_cb = plt.matplotlib.colors.BoundaryNorm(boundaries=bin_edges_sctt, ncolors=n_bins_sctt)
        
        # Scatter plot with discrete colormap
        sctt = ax.scatter(
            # eq_pos.ra, eq_pos.dec,
            # c=rm_full_replaced, 
            ra,dec,
            c=rm_m31_clipped, 
            marker='o', s=20,
            cmap=cmap_discrete_sctt, norm=norm_cb,
            transform=ax.get_transform('world'))

        # Colorbar for scatter
        if args.percent_cbar:
            cbar = plt.colorbar(sctt, ax=ax, spacing='proportional', boundaries=bin_edges_sctt, orientation='horizontal', fraction=0.05)
            cbar.ax.tick_params(labelrotation=90)
        else:
            cbar = plt.colorbar(sctt, ax=ax, orientation='horizontal', fraction=0.05)
        
        cbar.set_label(label=r'RM Scatter [$rad \text{ } m^{-2}$]', labelpad=labelpad, fontsize=cbar_lab_size)
        cbar.ax.tick_params(labelsize=l_s)

    # Create a scatter plot for positive and negative
    positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0

    ra_bg, dec_bg = eq_pos.ra[bg_condition], eq_pos.dec[bg_condition]
    # ra = np.concatenate((pos_in_region.ra, neg_in_region.ra))
    # dec = np.concatenate((pos_in_region.dec, neg_in_region.dec))

    # nsig is a multiplier for sigma (std) in y and x direction.
    im, kernel, imsize = smooth_2d_image(ra, dec, kernel=2, fitfile=filename, rm_m31=rm_m31_clipped)
    im_clipped = im #np.clip(im, -5000, 5000) #manually clipping based off of what i can see on defualt colorbar
    # x0s_pos, y0s_pos = ra_dec_to_pixels(ra[positive_mask], dec[positive_mask], filename=filename)
    # x0s_neg, y0s_neg = ra_dec_to_pixels(ra[negative_mask], dec[negative_mask], filename=filename)
    if args.bg: x0s_bg, y0s_bg = ra_dec_to_pixels(ra_bg, dec_bg, filename=filename)
    # #image 'im' is already trnsformed intowcs.ax's 'world'
    # ax.scatter(x0s_pos, y0s_pos, marker='o', color="b", s=2)
    # ax.scatter(x0s_neg, y0s_neg, marker='o', color="r", s=2)
    if args.bg: ax.scatter(x0s_bg, y0s_bg, marker='x', s=.05, color='k', alpha=.6)

    if args.smoothed:
        # plt.title(f"rad-limit = {N} deg", fontsize=20)
        n_bins=15
        # Get discrete colormap and bin edges

        #Creating masker and apply it for shaded grey region on outskirts of R_vir of M31
        masker = RMImageMasker(ra, dec, L_m31_deg=L_m31, fitfile=filename, m31_pos=m31_pos, outskirts_RM_value=vmin - 1000)
        im_shaded_outskirts = masker.apply_mask(im, imsize=imsize)

        cmap_discrete, bin_edges = get_discrete_colors((vmin, vmax), n_bins, 'delta', get_edges=True, use_cmocean=True)
        norm_cb = plt.matplotlib.colors.BoundaryNorm(boundaries=bin_edges, ncolors=n_bins)

        # Plot imshow
        # RM_dens = ax.imshow(im, cmap='jet_r', origin='lower', alpha=0.5, zorder=-1, vmin=vmin, vmax=vmax)
        # RM_dens = ax.imshow(im, cmap=cmap_discrete, norm=norm_cb, origin='lower', alpha=0.5, zorder=-1)
        RM_dens = ax.imshow(im_shaded_outskirts, cmap=cmap_discrete, norm=norm_cb, origin='lower', alpha=0.5, zorder=-1)

        # Colorbar
        RM_dens_cbar = plt.colorbar(RM_dens, ax=ax, fraction=0.05, orientation='horizontal')
        RM_dens_cbar.set_label(r'Smooth RM Intensity (4$\sigma$ clipped) [$rad \text{ } m^{-2}$]', labelpad=labelpad, fontsize=cbar_lab_size)
        RM_dens_cbar.ax.tick_params(labelsize=l_s)
        
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

        #Adding in filled region of grey jsut for the legend
        outskirts = plt.fill_between([],[], label=r"> R$_{vir}$", color="#888f9f")

        if args.add_gal_grid:

            #Adding an extra coordinate grid (e.g., Galactic)
            overlay = ax.get_coords_overlay('galactic')
            overlay.grid(color='k', ls='--', alpha=.4)
            overlay[0].set_axislabel('Galactic Longitude (l)', fontsize=f_s, fontname='Cantarell')
            overlay[1].set_axislabel('Galactic Latitude (b)', fontsize=f_s, fontname='Cantarell')

            l_axis: CoordinateHelper = overlay[0]
            b_axis: CoordinateHelper = overlay[1]
            l_axis.set_ticks(spacing=10 * u.deg)  #l tick spacing
            b_axis.set_ticks(spacing=10 * u.deg)  #b tick spacing
            l_axis.set_ticklabel(fontsize=tick_f_s)
            b_axis.set_ticklabel(fontsize=tick_f_s)

    #Marking M31 and M33 position
    color = color = 'white' if args.filled_485 else ("#0000ff" if args.scatter else ('#ff0000' if args.smoothed else None))
    sctt_m31 = ax.scatter(m31_pos.ra, m31_pos.dec, color=color, marker='s',
            transform=ax.get_transform('world'), label = "M31", s=100, alpha=1, zorder=10, edgecolor='k')
    sctt_m33 = ax.scatter(m33_pos.ra, m33_pos.dec, color=color, marker='^',
            transform=ax.get_transform('world'), label = "M33", s=100, alpha=1, zorder=10, edgecolor='k')

    more=.09 if args.smoothed else 0 #Parametre to move label legends up a bit

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

    ra_axis: CoordinateHelper = ax.coords[0]
    ra_axis.set_ticks(spacing=10 * u.deg)  #RA tick spacing
    ra_axis.set_major_formatter('d.d')  #Setting format of RA axis to degrees

    #ZOOMING OUT (in degrees)
    ra_min, ra_max = 45 if args.scatter else 40, -30
    dec_min, dec_max = 12 if args.scatter else 16, 48
    
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


    # List of possible handles (could include None)
    potential_handles = [
        sctt_m31,
        sctt_m33,
        _485mhz_cont if args._485mhz else None,
        HI_cont if args.HI else None,
        outskirts if args.smoothed else None,
    ]

    # Filter out any that are None
    handles = [h for h in potential_handles if h is not None]

    plt.legend(fontsize = 20, loc = 'upper center', bbox_to_anchor = (0.5, 1.1+more), 
                framealpha = 1 if args.filled_485 else 0, 
                ncols = len(handles), markerscale=1.1, facecolor="grey" if args.filled_485 else "None",
               handles=handles)

    #Inner Plot for scatter subplot
    if args.scatter:
        inset_ax = inset_axes(ax, width="22%", height="22%", bbox_to_anchor=(-0.72, -0.71, 1, 1), bbox_transform=ax.transAxes)
        simple_hist(rm_m31_clipped, #gauss_clipped_rm,
                    N=N, 
                    inset           =   inset_ax, 
                    n_bins          =   n_bins_sctt,
                    cmap_discrete   =   cmap_discrete_sctt, 
                    bin_edges       =   bin_edges_sctt)

    plt.tight_layout()
    if args.save_plot:
        path = curr_dir_path() + "Results/"
        # path = curr_dir_path() + "Results/Changing_CGM_4/"
        if args.smoothed:
            fname = f"2D_Gaussian_RM_Smoothing_kernel={tuple(np.round(kernel, 0))}_in_degrees.png"
            # fname = f"2D_Gaussian_RM_Smoothing_kernel={tuple(np.round(kernel, 0))}_in_degrees_{N}.png"
            plt.savefig(f"{path}{fname}", dpi=600, bbox_inches="tight")
            print(f"2D_Gaussian_RM_Smoothing has been saved to {path}{fname}")
        else:
            fname = f"Temperature_Contours_within_CGM_of_M31.png"
            plt.savefig(f"{path}{fname}", dpi=600)
            print(f"=Temperature Contour has been saved to {path}{fname}")
        plt.close()
    else:
        plt.show()

    # simple_hist(rm_m31_clipped[m31_condition], N=N)
    # simple_hist(gauss_clipped_rm, N=N) #NOTE: Gaussiuan clipping of points is only for smoothing visuals. Not for scatter plot

    break #Just for CGM of M31 (Nothing smnaller or larger)

