from main import (
    
#importing variables
L_m31, cutoff, rm_m31, rm_bg, rm_m31_coord, PA_bg, 
m31_sep_bg, PA_m33, m31_condition, m33_m31coord, PA_rm,
m31_sep_Rvir, m33_sep, m31_maj, m31_min, m31_pa,
  
#Importing function
update_projection, 
plot_ellipse
)

from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

# Fill a range of radial distances for background region
r_min, r_max = L_m31, cutoff.value  #Defining a radial range to fill
theta_fill = np.linspace(0, 2 * np.pi, 100)  #Coveringentire polar angle

                        #away #Neutral #Towards
cmap = ListedColormap(['red','none','blue'])

#Plotting all projections available
projections = ['polar']#[None, 'polar', 'lambert'] # Getting a list of all available 2D and 3D projection names

len_proj = len(projections)
figsize, add_leg, add_title = [(15,7.5), True, False] if len_proj == 1 else [(10,14), False,False] 

#Creating a figure with subplots for each projection
fig, ax = plt.subplots(nrows=len_proj, 
                        ncols=2, 
                        figsize=figsize)  

rm_m31_normalized = np.where(rm_m31 < 0, -1, np.where(rm_m31 > 0, 1, 0))
rm_bg_normalized = np.where(rm_bg < 0, -1, np.where(rm_bg > 0, 1, 0))

D = "N" #Starting direction for polar plot
#Iterating over each projection
for i, pro_i in enumerate(projections):
    #Updating first column subplot in current row
    
    #Updating second column subplot in current row
    ax1 = update_projection(ax, row_num = i, col_num = 1, 
                            projection = pro_i, 
                            fig=fig)
    if pro_i != "polar":       
        ax1.scatter(rm_m31_coord.lon.deg[m31_condition], 
                rm_m31_coord.lat.deg[m31_condition],
                c= rm_m31, cmap='Spectral', 
                s=2, marker='.',
                vmin=-100, vmax=100)
        ax1.plot(m33_m31coord.lon.deg, m33_m31coord.lat.deg, 'ks', label="m33")
        
        if pro_i == None:
            ax1.set_xlabel("galactic Longitude")
            ax1.set_ylabel("galactic Latitude")
            
    else:
        vmin, vmax = -1.5, 1.5
        # Fill a range of radial distances for background region
        ax1.fill_between(theta_fill, r_min, r_max, color='grey', 
                          alpha=.3)#, label='Background region')

        #plot rm within R_vir
        m31_rvir_sctt_norm = ax1.scatter(PA_rm, m31_sep_Rvir,
                    c= rm_m31_normalized, cmap=cmap, 
                    s=2, marker='.')
        
        m31_rvir_cbar = plt.colorbar(m31_rvir_sctt_norm, ax=ax1, cmap=cmap, 
                      orientation = "horizontal", 
                      fraction=.05)

        #Colorbars for none-normalized background and Rvir RM's are
        #stillsame due to vmin and vmax parametres

        ticks = [-1, -.66, 0, .66, 1]
        ticklabels = ['','-1', '0', '1','']
        m31_rvir_cbar.set_label(r"RM$_{Norm}$ [rad/$m^2$]", labelpad=1)
        m31_rvir_cbar.set_ticks(ticks)
        m31_rvir_cbar.set_ticklabels(ticklabels)
        
        #plot background scatter on grey region
        ax1.scatter(PA_bg, m31_sep_bg,
                    c= rm_bg_normalized, cmap=cmap, 
                    s=2, marker='.')
        
        #Mark M33
        ax1.plot(PA_m33, m33_sep.value, 'ks')#, label="m33")
        
        #Plotting major and minor axis
        mag_maj_axis = m31_maj.to(u.deg) #Magnitude of major axis of m31
        mag_min_axis = m31_min.to(u.deg) #Magnitude of minor axis of m31
        
        ax1.set_rlabel_position(115)
        ax1.set_theta_zero_location(D)  # 0 degrees attop (North)
        ax1.set_theta_direction(1)  # Degrees increase counterclockwise
    
    ax1.scatter([0], [0], color='y', marker='s', edgecolor='k', 
                   zorder=10)#,label="m31")
    
    if add_title: ax1.set_title(f'{pro_i} - Graph 2')

    ax2 = update_projection(ax, row_num = i, 
                        col_num = 0, 
                        projection = pro_i, fig=fig)
        
    if pro_i != "polar":
        ax2.scatter(rm_m31_coord.lon.deg[m31_condition], 
            rm_m31_coord.lat.deg[m31_condition], 
            c=rm_m31_normalized, marker='.', 
            s=2, cmap=cmap)
        ax2.plot(m33_m31coord.lon.deg, m33_m31coord.lat.deg, 'ks', label="m33")
        
        if pro_i == None:
            ax2.set_xlabel("galactic Longitude")
            ax2.set_ylabel("galactic Latitude")
    else:
        vmin, vmax = -100, 100
        cmap = 'Spectral'
        # Fill a range of radial distances for background region
        ax2.fill_between(theta_fill, r_min, r_max, color='grey', 
                          alpha=.3)#, label='Background region')

        #plot rm within R_vir
        m31_rvir_sctt = ax2.scatter(PA_rm, m31_sep_Rvir,
                    c=rm_m31, marker='.', 
                    s=2, cmap=cmap, vmin=vmin, vmax=vmax)
        
        #Colorbars for nomrlaized background and Rvir RM's are same
        #Hence no need to make a separate one.
        m31_rvir_cbar = plt.colorbar(m31_rvir_sctt, ax=ax2, cmap=cmap, 
                      orientation = "horizontal", 
                      fraction=.05)
        m31_rvir_cbar.set_label( r"RM [rad/$m^2]$", labelpad=1)
        
        #plot background region
        ax2.scatter(PA_bg, m31_sep_bg, marker='.',
                    c= rm_bg_normalized, 
                    s=2, cmap=cmap)
        
        #Mark M33
        ax2.plot(PA_m33, m33_sep.value, 'ks', label="m33")
        
        # Plotellipse representingmajor and minor axes based on its PA
        plot_ellipse(ax=ax2, major_axis = mag_maj_axis,
                      minor_axis = mag_min_axis , PA = m31_pa)
        
        ax2.set_rlabel_position(115)
        ax2.set_theta_zero_location(D)  # 0 degrees attop (North)
        ax2.set_theta_direction(1)  # Degrees increase counterclockwise
    
    ax2.scatter([0], [0], color='y', marker='s', edgecolor='k', 
                   zorder=10, 
                   label="m31")
       
    if add_title: ax2.set_title(f'{pro_i} - Graph 1')
    
#Adjusting subplot layout to fit intofigure
if add_leg: 
    fig.legend(fontsize=16, loc="upper left", 
               framealpha=0, ncols=1, borderpad=.4)
plt.tight_layout(pad=0)
plt.show()

#Plotting all projections available
projections = ['polar']#[None, 'polar', 'lambert'] # Getting a list of all available 2D and 3D projection names

fig, ax = plt.subplots(nrows=len_proj, 
                        ncols=2, 
                        figsize=(15,7.5))  

#Iterating over each projection
for i, pro_i in enumerate(projections):
    #Updating first column subplot in current row
    
    ax1 = update_projection(ax, row_num = i, 
                        col_num = 0, 
                        projection = pro_i, fig=fig)

    #Plotting major and minor axis
    mag_maj_axis = m31_maj.to(u.deg) #Magnitude of major axis of m31
    mag_min_axis = m31_min.to(u.deg) #Magnitude of minor axis of m31
    
    # Plotellipse representingmajor and minor axes based on its PA
    plot_ellipse(ax=ax1, major_axis = mag_maj_axis,
                  minor_axis = mag_min_axis , PA = m31_pa,
                  ax_small=True)
    
    ax1.set_theta_zero_location(D)  # 0 degrees attop (North)
    ax1.set_theta_direction(1)  # Degrees increase counterclockwise
    
    #MArk m31 with yellow swuare encapsulated with black oultine.
    ax1.scatter([0], [0], color='y', marker='s', edgecolor='k', 
                zorder=10, 
                label="m31")
    
    # ax1.set_title("Zoomed in without RM scatter")
    ax1.set_rlabel_position(293)
    
    #Updating second column subplot in current row
    ax2 = update_projection(ax, row_num = i, col_num = 1, 
                            projection = pro_i, 
                            fig=fig)

    # Fill a range of radial distances for background region
    ax2.fill_between(theta_fill, r_min, r_max, color='grey', 
                      alpha=.3)#, label='Background region')

    #plot rm within R_vir
    m31_rvir_sctt = ax2.scatter(PA_rm, m31_sep_Rvir,
                c= rm_m31, cmap='Spectral', 
                s=2, marker='.',
                vmin=-100, vmax=100)
    
    #plot background region
    ax2.scatter(PA_bg, m31_sep_bg,
                c= rm_bg, cmap='Spectral', 
                s=2, marker='.',
                vmin=-100, vmax=100)
    
    #Colorbars for none-normalized background and Rvir RM's are
    #stillsame due to vmin and vmax parametres
    m31_rvir_cbar = plt.colorbar(m31_rvir_sctt, ax=ax2, cmap=cmap, 
                  orientation = "horizontal", 
                  fraction=.05)
    m31_rvir_cbar.set_label( "RM values", labelpad=1)
    
    
    #Plotting major and minor axis
    mag_maj_axis = m31_maj.to(u.deg) #Magnitude of major axis of m31
    mag_min_axis = m31_min.to(u.deg) #Magnitude of minor axis of m31
    
    plot_ellipse(ax=ax2, major_axis = mag_maj_axis,
                  minor_axis = mag_min_axis , PA = m31_pa)
    
    ax2.set_rlabel_position(115)
    ax2.set_theta_zero_location(D)  # 0 degrees attop (North)
    ax2.set_theta_direction(1)  # Degrees increase counterclockwise
    
    #Mark M33 with star
    ax2.plot(PA_m33, m33_sep.value, 'ks')#, label="m33")

    # ax2.scatter([0], [0], color='y', marker='s', edgecolor='k', 
    #             zorder=10, 
    #             label="m31")
    
    ax2.set_title("Zoomed out with RM scatter")

if add_leg: 
    fig.legend(fontsize=12, bbox_to_anchor=(.47, 1.01), 
               framealpha=0, ncols=7, borderpad=.4)
plt.tight_layout(pad=3.5)
plt.show()