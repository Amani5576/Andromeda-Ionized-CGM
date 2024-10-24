from main import (
    
# Importing functions
confining_rectangle,

#Importing package alias'
fits, WCS, plt, np, u, Table, SkyCoord,

#Importing variables
rm_m31, m31_pos, m31_pa, rm_pos_icrs, err_m31
)

def plot_patch_on_m31():

    levels = 16 #For HI density contour

    #Opening FITS file and extracting data
    filename = 'LGSNLGSR.SQLGBB.FITS'
    hdu = fits.open(filename)[0]

    #Setting up WCS projection usingheader
    wcs = WCS(hdu.header)

    #Creating figure and axis with WCS projection
    fig = plt.figure(figsize=(9.5, 10))
    ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y', 0, 0))

    #Adding contours
    ax.contour(hdu.data[0, 0, :, :], levels=levels,
               linewidths=1, colors='black',
               linestyles='solid')

    ax.grid(ls='solid', alpha=.4, color="black")
    f_s, tick_f_s = 16, 10
    ax.coords[0].set_axislabel('Right Ascension (J2000)', fontsize=f_s)
    ax.coords[0].set_ticklabel(fontsize=tick_f_s)
    ax.coords[1].set_axislabel('Declination (J2000)', fontsize=f_s)
    ax.coords[1].set_ticklabel(fontsize=tick_f_s)
    ax.tick_params(axis='x', labelsize=f_s/1.2)
    ax.tick_params(axis='y', labelsize=f_s)

    # Create a scatter plot for positive values
    positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0

    # This is done in assessing_boxScatter_with_HI() but plotting to show box and line is done here.
    #Region away form filament (focusing on regions closer to m31)
    #Creating scatter plot that is transformed uses these variables.
    (filament_region_coord_pos, #To be used for later after plot() is called
     rm_val_filament_pos, #To be used for later after plot() is called
     rm_err_filament_pos, #To be used for later after plot() is called
     filament_region_coord_neg, #To be used for later[] after plot() is called
     rm_val_filament_neg,  #To be used for later after plot() is called
     rm_err_filament_neg,  #To be used for later after plot() is called
     rect_widthMid  #mid points of width of rectangular patch
     ) = confining_rectangle(ax, marker_size=40,
                            ra=m31_pos.ra.deg+6.4,
                            dec=m31_pos.dec.deg-1,
                            width=10, height=17, angle=7, 
                            polar_angle=m31_pa.value, 
                            positions=rm_pos_icrs, 
                            pos_mask=positive_mask, 
                            neg_mask = negative_mask,
                            rm_s=rm_m31, 
                            rm_errs=err_m31,
                            return_data=True,
                            plot=(True, #plotting box 
                                 True) #plotting scatter within box
                           )
    
    #Plot line passing through the Rectangular_filament region
    # This should be commented out if the plotting of the box is false
    ax.plot(np.array(rect_widthMid[0])-.2, rect_widthMid[1], color='k', linestyle='-',
            transform=ax.get_transform("world"), linewidth=2)# label="Custom axis")
    
    # Adjust the origin for the arrows
    origin_x = m31_pos.ra.deg-10.8
    origin_y = m31_pos.dec.deg+2.2
    
    arrow_length = 2
    
    # Adding the orthogonal arrows onto the WCS axis
    hw = .5
    ax.arrow(origin_x, origin_y, arrow_length, -1.6, color='red', head_width=hw, transform=ax.get_transform("world"), zorder=2)
    ax.arrow(origin_x, origin_y, 1.9, arrow_length, color='red', head_width=hw, transform=ax.get_transform("world"), zorder=2)
    
    # Adding labels for x' and y'
    ax.text(origin_x +1, origin_y-1.5, "$x'$", fontsize=20, color='red', transform=ax.get_transform("world"))
    ax.text(origin_x+1, origin_y+1.6 , "$y'$", fontsize=20, color='red', transform=ax.get_transform("world"))
    
    
    from astropy.visualization.wcsaxes import CoordinateHelper
    ra_axis: CoordinateHelper = ax.coords[0]
    ra_axis.set_ticks(spacing=2 * u.deg)  #RA tick spacing
    # ra_axis.set_major_formatter('d.d')  #Setting format of RA axis to degrees

    #ZOOMING (in degrees)
    ra_min, ra_max = 19, -3
    dec_min, dec_max = 31, 50
 
    # Convert these to pixel coordinates using WCS
    bottom_left = wcs.world_to_pixel_values(ra_min, dec_min, 0, 0)
    top_right = wcs.world_to_pixel_values(ra_max, dec_max, 0, 0)
   
    #Extracting pixel coordinates
    x_min, y_min, _, _ = bottom_left
    x_max, y_max, _, _ = top_right

    #Setting axis limits
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    # plt.legend(fontsize = 20, loc = 'upper center', bbox_to_anchor = (0.5, 1.12), 
    #             framealpha = 0, ncols = 3, markerscale=1.5)
    #            handles=[sctt_m31, sctt_m33, sctt_bg])
    # frame = legend.get_frame()
    # frame.set_facecolor('gray')
    
    plt.tight_layout()
    print("main complete")
    plt.show()
        
plot_patch_on_m31()

