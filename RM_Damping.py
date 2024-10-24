from main import (
    
# Importing functions
confining_rectangle, 
axis_transformation, 
binned_scatter, 
fit_and_plot_damped_sine_wave,

#Importing package alias'
fits, WCS, plt, np,

#Importing variables
rm_m31, m31_pos, m31_pa, rm_pos_icrs, err_m31

)

filename = 'LGSNLGSR.SQLGBB.FITS'
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y', 0, 0))
positive_mask = rm_m31 > 0 ; negative_mask = rm_m31 < 0

(filament_region_coord_pos,
 rm_val_filament_pos,
 rm_err_filament_pos,
 filament_region_coord_neg,
 rm_val_filament_neg,
 rm_err_filament_neg,
 rect_widthMid  #mid points of width of rectangular patch
 ) = confining_rectangle(ax, 
                        ra=m31_pos.ra.deg+6,
                        dec=m31_pos.dec.deg-1.5,
                        width=10, height=17, angle=7, 
                        polar_angle=m31_pa.value, 
                        positions=rm_pos_icrs, 
                        pos_mask=positive_mask, 
                        neg_mask = negative_mask,
                        rm_s=rm_m31, 
                        rm_errs=err_m31,
                        return_data=True,
                        plot=(False, #plotting box 
                             False) #plotting scatter within box
                       )
plt.close() #removing figure from showing since only axis was needed here.
def plot_binned(ax, coord, rm_values, rm_errs, bins, axis_label, combine, plot_type='scatter', markersize=80, x_prime=False, **kw):
    
    if combine: #Setting up combined verisons of the binned data
        coord = np.append(coord[0],coord[1])
        rm_values = np.append(rm_values[0], rm_values[1])
        rm_errs = np.append(rm_errs[0], rm_errs[1])

        binned_xy, binned_err, _ = binned_scatter(coord, rm_values, bins)
        binned_prime_axis, binned_RM = binned_xy[:,0], binned_xy[:,1]
    else:
        binned_pos, binned_pos_err, _ = binned_scatter(coord[0], rm_values[0], bins)
        binned_neg, binned_neg_err, _ = binned_scatter(coord[1], rm_values[1], bins)
        
    if plot_type == 'scatter':
        if combine:
            # print(binned_prime_axis[5]); import sys; sys.exit()
            ax.scatter(binned_prime_axis, binned_RM, marker='s', color="white", edgecolors="b", s=markersize)
            ax.errorbar(binned_prime_axis, binned_RM, yerr=binned_err, fmt="", markersize=3, capsize=3, ecolor='red')
        else:
            ax.scatter(binned_pos[:,0], binned_pos[:,1], marker='s', label='Positive RM', color="blue", s=markersize)
            ax.scatter(binned_neg[:,0], binned_neg[:,1], marker='s', label='Negative RM', color='red', s=markersize+2)
            ax.errorbar(binned_pos[:,0], binned_pos[:,1], yerr=binned_pos_err, fmt="", markersize=3, capsize=3, ecolor='blue')
            ax.errorbar(binned_neg[:,0], binned_neg[:,1], yerr=binned_neg_err, fmt="", markersize=3, capsize=3, ecolor='red')
            
        #Drawing x-axis line at y=0
        ax.axhline(y=0, color='black', linewidth=1)
        
        #Drawing x-axis line at y=minimum y value to enclose the bottom region
        wf = kw["width_factor"] if "width_factor" in kw else 0
        y_minium = -60 if wf ==.14 else -80
        x_lab_pad = 220 if wf ==.14 else 260
        ax.axhline(y=y_minium, color='black', linewidth=1)
        
        #Hiding right and top spines
        ax.spines['bottom'].set_position('zero')
        ax.set_xlabel(axis_label, fontsize=25, labelpad=x_lab_pad)
        zoomout_x = 1
        
    elif plot_type == 'bar':
        zoomout_x = 2
        if combine:
            wf = kw["width_factor"] if "width_factor" in kw else .8
            width = (binned_prime_axis[1] - binned_prime_axis[0]) * wf  # Bar width
            X = binned_prime_axis - width/2
            ax.bar(X, binned_RM, width=width, color='gray',  edgecolor='k')
            xmin, xmax = min(binned_prime_axis), max(binned_prime_axis)
            if wf == .14:  #for x prime
                initial_guesses = [
                    [60, 0.01, 0.6, 0.2]
                ]
                
                for initial_guess in initial_guesses:
                    fit_and_plot_damped_sine_wave(X, binned_RM, initial_guess, ax, color="r", prime='x')
                # fit_and_plot_damped_sine_wave(binned_prime_axis, binned_RM, [50, 0.1, 0.5, 0], ax)
            else:
                
                initial_guesses = [
                    [5.92, 0.0123, 0.5, 0],
                    [25, 0.02, 0.8, 1],
                ]
                col = ["r", "b"]
                
                for idx, initial_guess in enumerate(initial_guesses):
                    fit_and_plot_damped_sine_wave(X, binned_RM, initial_guess, ax, color=col[idx], prime='y')

                # fit_and_plot_damped_sine_wave(binned_prime_axis, binned_RM, [45, -0.05, 1, 0], ax)
            
        else:
            width = (binned_pos[:,0][1] - binned_pos[:,0][0]) * wf  # Bar width
            ax.bar(v1:=binned_pos[:,0] - width/2, binned_pos[:,1], width=width, label='Positive RM', color='blue', edgecolor='blue')
            ax.bar(v2:=binned_neg[:,0] - width/2, binned_neg[:,1], width=width, label='Negative RM', color='red', edgecolor='red', alpha=0.7)
            xmin=np.min(np.append(v1,v2))-zoomout_x
            xmax=np.max(np.append(v1,v2))+zoomout_x
        
        #Adding horizontal line to simply show the where RM=0
        ax.hlines(y=0, xmin=xmin, xmax=xmax, color='k', linewidth=.6)
        ax.set_xlabel(axis_label, fontsize=25)
        
    if combine:
        ax.set_xlim(min(binned_prime_axis)-zoomout_x, 
                    max(binned_prime_axis)+zoomout_x)
    else:
        ax.set_xlim(min(np.min(binned_pos[:,0]), np.min(binned_neg[:,0]))-zoomout_x, 
                    max(np.max(binned_pos[:,0]), np.max(binned_neg[:,0]))+zoomout_x)
    ax.minorticks_on()
    ax.set_ylabel("RM "+"[rad/$m^2$]", fontsize=20)
    ax.tick_params(labelsize=14)
    
    if x_prime:
        ax.legend(fontsize=16, loc='upper right', framealpha=1, ncol=3)
    else:
        ax.legend(fontsize=16, loc='upper left', framealpha=1, ncol=1)
        
    plt.grid(True)

def create_plots(plot_type, combine=False): 
    #If combine == True then then all points (regardless of positive/negative) 
    # will be treated as one list of elements rather than searate negative and positive points

    plt.figure(figsize=(15, 15))
    bins = 20

    # plt.suptitle(f'boxed region focused around m31(max {bins=})', fontsize=25)
    # plt.suptitle(f'Within Filament Region ({bins=})', fontsize=35)

    #Preparing the retreived coordinates to be 
    #used properly in axis_transformation() function
    pos_coords = np.vstack([filament_region_coord_pos.ra.deg,
                            filament_region_coord_pos.dec.deg]).T

    neg_coords = np.vstack([filament_region_coord_neg.ra.deg,
                            filament_region_coord_neg.dec.deg]).T

    #Getting coordinates and axis after transformtion
    #DEC = -RA * tan(angle) #for x'
    #DEC = -RA * cot(angle) #for y'
    (new_x_pos, new_y_pos,
     new_RA_points_pos, new_DEC_points_pos,
     angle
     ) = axis_transformation(points=pos_coords,
                            RA_range=rect_widthMid[0],
                            DEC_range=rect_widthMid[1]
                            )
    (new_x_neg, new_y_neg,
     new_RA_points_neg, new_DEC_points_neg,
     angle
     ) = axis_transformation(points=neg_coords,
                            RA_range=rect_widthMid[0],
                            DEC_range=rect_widthMid[1]
                            )

    #Get M31 position in transformed axis of RA (x-prime)
    #Getting coordinates and axis after transformtion
    (_,_,
     new_RA_point_m31, new_DEC_point_m31, _
     ) = axis_transformation(
         points=[m31_pos.ra.deg, m31_pos.dec.deg],
                RA_range=rect_widthMid[0],
                DEC_range=rect_widthMid[1]
                )

    #Plotting RM values against their RA (in new rotated axis)
    ax1 = plt.subplot(211)

    # print(new_RA_point_m31, new_DEC_point_m31); import sys; sys.exit()
    #plotting line to indicate position of M31 along x-prime
    ax1.axvline(x=new_RA_point_m31, 
                ymin=np.min(np.append(rm_val_filament_pos, rm_val_filament_neg)),
                ymax=np.max(np.append(rm_val_filament_pos, rm_val_filament_neg)), 
                linestyle="--", color='k')#, label="m31")
    
    plot_binned(ax1, (new_RA_points_pos, new_RA_points_neg), 
                (rm_val_filament_pos, rm_val_filament_neg), 
                (rm_err_filament_pos, rm_err_filament_neg), 
                bins=bins, axis_label="x\'", #"$ or $DEC = -RA \cdot \\tan( {np.rad2deg(angle):.1f} )$", 
                plot_type=plot_type, x_prime=True, 
                combine=combine, width_factor=.14)

    #Plotting RM values against their Dec
    ax2 = plt.subplot(212)

    #plotting line to indicate position of M31 along y-prime
    ax2.axvline(x=new_DEC_point_m31, 
                ymin=np.min(np.append(rm_val_filament_pos, rm_val_filament_neg)),
                ymax=np.max(np.append(rm_val_filament_pos, rm_val_filament_neg)), 
                linestyle="--", color='k')

    plot_binned(ax2, (new_DEC_points_pos, new_DEC_points_neg), 
                (rm_val_filament_pos, rm_val_filament_neg), 
                (rm_err_filament_pos, rm_err_filament_neg), 
                bins=bins, axis_label="y\'", #"$ or $DEC = RA \cdot \\cot( {np.rad2deg(angle):.1f} )$", 
                plot_type=plot_type, combine=combine, width_factor=.8)

    plt.tight_layout(pad=5)
    plt.show()

create_plots(plot_type='scatter', combine=True)
create_plots(plot_type='bar', combine=True)