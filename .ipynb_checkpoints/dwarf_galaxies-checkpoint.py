from astropy.coordinates import Galactocentric
from astropy.utils.exceptions import AstropyWarning
import warnings

# Suppress specific Astropy warnings
warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from personal_converter_copy import generate_marker_color_pairs, convert_string_where_applicable
from main import (

#importing functions    
get_angle_sep_from_distance, 
get_projected_d, 

#importing alias'
SkyCoord, u, np, stats, plt, Table,

#Importing varibales
rm, rm_err, eq_pos, L_m31, m31_pos
)

import pickle

dwarf_galaxies = Table.read("dwarf_galaxies.dat", format='csv')

convert_string_where_applicable(dwarf_galaxies)
Dwf_Gs_Loc = SkyCoord(ra = dwarf_galaxies["RA"], dec=dwarf_galaxies["DEC"], 
                  unit=(u.deg, u.deg), frame="icrs")

# Set units for specific columns
dwarf_galaxies["l_MS"].unit = u.deg
dwarf_galaxies["b_MS"].unit = u.deg
dwarf_galaxies["RA"].unit = u.hourangle
dwarf_galaxies["DEC"].unit = u.deg
dwarf_galaxies["D_dwarf"].unit = u.kpc
dwarf_galaxies["R"].unit = u.kpc
dwarf_galaxies["X"].unit = u.kpc
dwarf_galaxies["Y"].unit = u.kpc
dwarf_galaxies["M_*"].unit = 10**5 * u.M_sun
dwarf_galaxies["M_200"].unit = 10**8 * u.M_sun
dwarf_galaxies["R_200"].unit = u.kpc
dwarf_galaxies["v_200"].unit = u.km / u.s
dwarf_galaxies["v_LSR"].unit = u.km / u.s
dwarf_galaxies["v_m31"].unit = u.km / u.s

# In[]:
# =============================================================================
# #This shows that the Radius to M31 is not always larger than R_200
# dwarf_galaxies["R"] >  dwarf_galaxies["R_200"]
# array([False, False,  True,  True, False,  True,  True,  True,  True,
#         True,  True,  True,  True,  True,  True,  True,  True,  True,
#         True,  True,  True,  True,  True,  True,  True,  True])
# =============================================================================


def get_dwarf_galaxy_stats(eq_pos, rm, dwarf_galaxies):
    
    """This function assesses dwarf galaxy stats within 'R_200'.
    
    return_projection_data if True returns lists of each dwarf containing both 
    projected separation and RM for plotting (similar to M31).
    """
    
    #Creating duplicates 
    eq_pos_s = [eq_pos for i in range(len(dwarf_galaxies))] # In case it's mutable
    
    #Getting dwarf galaxy coordinates in terms of their RA and Dec
    dwarf_frame = Dwf_Gs_Loc.skyoffset_frame()  # A bunch of offset frames from each individual dwarf galaxy
    
    # Calculating the angular separations of each RM coordinate relative to dwarf galaxy
    dwarf_Angular_Separations_s = [eq_p.separation(Dwf_Gs_Loc[i]) for i, eq_p in enumerate(eq_pos_s)]
    
    # Get radial limit (R_200) in degrees of separation
    R_200_limits = get_angle_sep_from_distance(R_proj=dwarf_galaxies["R_200"],
                                               d=dwarf_galaxies["D_dwarf"])  # In degrees away from center
    
    # Filter dwarf angular separations based on R_200 limits in degrees
    dwarf_ang_sep_within_R_200_limit_s = [
        d_ang_sep[d_ang_sep.deg < R_200_limits[i].value]
        for i, d_ang_sep in enumerate(dwarf_Angular_Separations_s)
    ]
    
    #Background (Amount of background is similar in ratio to that of M31)
    cutoff_dwf_G = 30 * (R_200_limits / L_m31)
    
    #Getting background RM's for each dwarf galaxy
    dwarf_ang_sep_bg_s = [
        d_ang_sep[(d_ang_sep.deg > R_200_limits[i].value) & (d_ang_sep.deg < cutoff_dwf_G[i].value)]
        for i, d_ang_sep in enumerate(dwarf_Angular_Separations_s)
    ]

    dwarf_rm_coords = [eq_p.transform_to(dwarf_frame[i]) for i, eq_p in enumerate(eq_pos_s)]
    
    # Filter RM coordinates near dwarf within relative R_200 limit
    dwarf_rm_coords_within_R_200 = [
        dw_rm_coords[dw_ang_sep.deg < R_200_limits[i].value]
        for i, (dw_ang_sep, dw_rm_coords) in enumerate(list(zip(dwarf_Angular_Separations_s, dwarf_rm_coords)))
    ]
    
    dwrf_G_condition_s = [d_ang_sep.deg < R_200_limits[i].value
                          for i, d_ang_sep in enumerate(dwarf_Angular_Separations_s)
                          ]
    dwrf_G_bg_condition_s = [(d_ang_sep.deg > R_200_limits[i].value) & (d_ang_sep.deg < cutoff_dwf_G[i].value)
                             for i, d_ang_sep in enumerate(dwarf_Angular_Separations_s)]

    # Applying conditions to filter RM positions and making sure they are in ICRS format
    dwrf_G_bg_pos_s = [eq_pos[dwrf_G_bg_condition] for dwrf_G_bg_condition in dwrf_G_bg_condition_s] 
    dwrf_G_pos_s = [eq_pos[dwrf_G_condition] for dwrf_G_condition in dwrf_G_condition_s]
    dwrf_G_bg_pos_icrs = [dwrf_G_bg_pos.transform_to("icrs") for dwrf_G_bg_pos in dwrf_G_bg_pos_s] 
    dwrf_G_pos_icrs = [dwrf_G_pos.transform_to("icrs") for dwrf_G_pos in dwrf_G_pos_s]
    
    # Applying the same conditions to get the RM values per dwarf
    dwrf_G_bg_rm_val_s = [rm[dwrf_G_bg_condition] for dwrf_G_bg_condition in dwrf_G_bg_condition_s] 
    dwrf_G_rm_val_s = [rm[dwrf_G_condition] for dwrf_G_condition in dwrf_G_condition_s]
    
    # Applying the same conditions to get the ERRORS of the RM values per dwarf
    dwrf_G_bg_rm_err_s = [rm_err[dwrf_G_bg_condition] for dwrf_G_bg_condition in dwrf_G_bg_condition_s] 
    dwrf_G_rm_err_s = [rm_err[dwrf_G_condition] for dwrf_G_condition in dwrf_G_condition_s]

    # Making sure all sub-elements are of type list and flattened to a single list of sublists
    dwrf_G_bg_rm_val_s_flat = [list(sublist) for sublist in dwrf_G_bg_rm_val_s]
    dwrf_G_rm_val_s_flat = [list(sublist) for sublist in dwrf_G_rm_val_s]
    
    dwrf_G_rm_bg_mean_s = [np.mean(sublist) for sublist in dwrf_G_bg_rm_val_s]  # Getting mean RM of each background RM of dwarf 
    dwrf_G_rm_bg_mean_s = np.nan_to_num(dwrf_G_rm_bg_mean_s)  # Converting any NaN to 0
    
    # Correction by background subtraction for each dwarf galaxy
    dwrf_G_rm_val_s = [dwrf_G_rm_val - dwrf_G_rm_bg_mean_s[i] for i, dwrf_G_rm_val in enumerate(dwrf_G_rm_val_s)]
    dwrf_G_rm_mean_s = [np.mean(sublist) for sublist in dwrf_G_rm_val_s]  # Getting mean RM of each dwarf within R_200_limit 
    dwrf_G_rm_mean_s = np.nan_to_num(dwrf_G_rm_mean_s)  # Converting any NaN to 0
    dwrf_G_rm_std_s = [np.std(sublist) for sublist in dwrf_G_rm_val_s]  # Getting standard deviation of RM within R_200_limit
    dwrf_G_rm_std_s = np.nan_to_num(dwrf_G_rm_std_s)  # Converting any NaN to 0
    
    # Angular separation within R_200_limit
    dwrf_G_sep_R_200_limit = [(Ang_Sep.deg[dwrf_G_condition]) * u.deg 
                              for Ang_Sep, dwrf_G_condition in
                              list(zip(dwarf_Angular_Separations_s, dwrf_G_condition_s))]
    
    # Angular separation for points between R_200_limit and background limit
    dwrf_G_sep_bg_R_200_limit = [(Ang_Sep.deg[dwrf_G_bg_condition]) * u.deg 
                              for Ang_Sep, dwrf_G_bg_condition in
                              list(zip(dwarf_Angular_Separations_s, dwrf_G_bg_condition_s))]

    return (dwrf_G_rm_mean_s, dwrf_G_rm_std_s,
        dwrf_G_rm_val_s, dwrf_G_bg_rm_val_s, 
        dwrf_G_sep_R_200_limit, dwrf_G_sep_bg_R_200_limit)


# =============================================================================
# #Checking whether tolerance value gives accurate RM collection
# print([len(dwarf_rm_coords_within_R_200[i]) for i in range(len(dwarf_rm_coords_within_R_200))])
# [191, 188, 13, 22, 57, 17, 38, 16, 32, 31, 117, 113, 13, 12, 22, 19, 10, 14, 9, 24, 10, 33, 53, 49, 21, 7]
# print([len(rm_within_R_200[i]) for i in range(len(rm_within_R_200))])
# [191, 188, 13, 22, 57, 17, 38, 16, 32, 31, 117, 113, 13, 12, 22, 19, 10, 14, 9, 24, 10, 33, 53, 49, 21, 7]
# =============================================================================

# In[]:

#These can take some time and that is why its been pickled after analysis below
#Its been commented out cause the data can just be accesed with picking instead of rerunning the code.

(rm_within_R_200_means, rm_within_R_200_stds,
dwrf_G_R_200_rm_val_s, dwrf_G_R_200_bg_rm_val_s,
dwrf_G_sep_R_200, dwrf_G_sep_bg_R_200) = get_dwarf_galaxy_stats(eq_pos, rm, dwarf_galaxies)

#Creating file name to store the pickled data of dwarf statistical data
data = (rm_within_R_200_means, rm_within_R_200_stds) #Virial Radius of dwarf galaxy

#dumping the data (Pickling the data)
with open("rm_within_R_200_data.pkl", 'wb') as f: #opening in binary mode
    pickle.dump(data, f)
    print(f"Data successfully pickled to {f}")

#Loading the Pickle data
with open("rm_within_R_200_data.pkl", 'rb') as f: #opening in binary mode
    data = pickle.load(f)
    print(f"Pickled data successfully loaded from {f}")
rm_within_R_200_means, rm_within_R_200_stds = data

# In[]:
    
def plot_dwarf_galaxy_analysis():

    def get_plotting_specs_for_dwarf(ax, x_data, xlabel, show_legend=False, is_first_col=False):
        
        if "M" in xlabel: #If dealing with Mass
            if "200" in xlabel:
                x_data *= 1e8 #Get it into just units of solar mass alone
            else:#For stellar mass
                x_data *= 1e5 #Get it into just units of solar mass alone
            x_data=np.log10(x_data)
                
        #Set the y-label only for the first column
        if is_first_col:
            ax.set_ylabel(r"Mean RM [rad/m$^{-2}$]")
        # else:
            # ax.set_yticklabels([])
        
        ax.grid(True)
        
        for i, name in enumerate(dwarf_galaxies["Name"]):
            color = "red" if name in ["LGS 3", "IC 10", "Lac I"] else "k"
            marker = "^" if name == "LGS 3" else "s"
            if name == "Lac I": marker = "x"
            ax.scatter(x_data[i], rm_within_R_200_means[i], color=color, marker=marker)
    
        # Plotting data within R_200
        ax.errorbar(x_data, rm_within_R_200_means, yerr=rm_within_R_200_stds, 
        fmt="none", ecolor=(1, 0, 0, .3), capsize=2.5)#, label="Dwarf G " + "$(R_{200})$")
        
        # Adding a horizontal line at y=0
        ax.axhline(0, color='k', linestyle='-', linewidth=1)
        
        # Adding the legend at the top of the figure
        if show_legend:
            ax.figure.legend(loc='upper center', fontsize='large', 
                             ncol=2, bbox_to_anchor=(0.5, 1.05), framealpha=0)
         
        # if "M" == xlabel[0]: ax.set_xscale('log')
        ax.set_xlabel(xlabel)
    
    fig, axs = plt.subplots(2, 2, figsize=(15, 7))
    get_plotting_specs_for_dwarf(axs[0, 0], dwarf_galaxies["M_200"], r"log($M_{200}$) [M$_{\odot}$]", show_legend=True, is_first_col=True)
    get_plotting_specs_for_dwarf(axs[0, 1], dwarf_galaxies["v_200"], r"v$_{200}$ [km/s]", show_legend=False)
    get_plotting_specs_for_dwarf(axs[1, 0], dwarf_galaxies["v_m31"], r"v$_{m31}$ [km/s]", show_legend=False, is_first_col=True)
    get_plotting_specs_for_dwarf(axs[1, 1], dwarf_galaxies["M_*"], r"log(M$_*$) [M$_{\odot}$]", show_legend=False)
    
    plt.tight_layout()
    plt.show()

def best_fit_line_plot(x, y, y_err, ax, get_grad=False, just_plot=False, show_residuals=False, **kw):
    if isinstance(x, u.Quantity):
        x = x.value

    coeff_s, cov_matrix = np.polyfit(x, y, 1, cov=True)  # 1 for linear fit
    best_fit_line = np.polyval(coeff_s, x)
    m, c = coeff_s

    perr = np.sqrt(np.diag(cov_matrix))
    m_err, c_err = perr

    if get_grad:
        return m, m_err

    if just_plot:
        plt.plot(x, best_fit_line, color='r', linestyle=":",
                 label=r"$RM$ = " + f"({m:.2g}"+r"$\pm$" +f"{m_err:.0g})" + r"$r_{Dwarf}$" + f" + ({c:.0g}" + r"$\pm$"+ f"{c_err:.0g})",
                 zorder=2, alpha=0.6, **kw)
        return

    ax = kw.get("ax", plt.gca())
    ax.plot(x, best_fit_line, color='r', linestyle=":",
            label=r"$RM$ = " + f"({m:.2g}"+r"$\pm$" +f"{m_err:.0g})" + r"$r_{Dwarf}$" + f" + ({c:.0g}" + r"$\pm$"+ f"{c_err:.0g})", 
            zorder=2, alpha=0.6, **kw)

    if show_residuals:
        residuals = y - best_fit_line
        ax.errorbar(x, residuals, yerr=y_err, fmt='o', color='b', alpha=0.5,
                    label='Residuals')
        ax.axhline(0, color='k', linestyle='--', linewidth=0.7)  # Zero line

    ax.legend()

# In[The 4 subplots]:
plot_dwarf_galaxy_analysis() #Plotting the M_200, V_m31, M_*, etc against RM

# In[]:
    
def get_bin_mean_std(X, Y, bin_num=22, get_edges=True, norm=False, **kw):
    # Using standard error of the mean for error bars within the R_200 limit
    bin_stds, _, _ = stats.binned_statistic(X, Y,
                                            statistic=lambda x: stats.sem(x), 
                                            bins=bin_num)  # Standard Error of Mean

    bin_means, bin_edges, _ = stats.binned_statistic(X, Y, 
                                                    statistic='mean', 
                                                    bins=bin_num)
    if norm:
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    else:
        bin_centers = get_projected_d(0.5 * (bin_edges[1:] + bin_edges[:-1]) * u.deg, kw["D"])
    
    L = bin_means, bin_stds, bin_edges, bin_centers
    if get_edges: return L[:-1]#just the bin_edges for unique bin_centres calculation
    return L 
    
Dwarf_Plots=[]
# Loop through each sublist for storing
for i in range(len(dwrf_G_sep_R_200)):
    
    # Check if arrays are empty and skip iteration if they are
    if len(dwrf_G_sep_R_200[i]) == 0 or len(dwrf_G_R_200_rm_val_s[i]) == 0:
        print(f"Skipping index {i} due to empty R_200 data.")
        continue
    
    R_200 = dwarf_galaxies["R_200"][i]
    D = dwarf_galaxies["D_dwarf"][i] * u.kpc
    
    # Using standard error of the mean for error bars within the R_200 limit
    (bin_means_R_200, 
     bin_std_R_200, 
     bin_edges_R_200) = get_bin_mean_std(X=dwrf_G_sep_R_200[i], 
                                         Y=dwrf_G_R_200_rm_val_s[i], 
                                         bin_num=22, D=D)


    # Calculating the projected distances for R_200
    d_rm_R_200 = get_projected_d(dwrf_G_sep_R_200[i], D)
    d_bg_R_200 = get_projected_d(dwrf_G_sep_bg_R_200[i], D)
    d_bin_centers_R_200 = get_projected_d(0.5 * (bin_edges_R_200[1:] + bin_edges_R_200[:-1]) * u.deg, D)

    
    Name = dwarf_galaxies["Name"][i]
    dist_to_m31 = dwarf_galaxies["R"][i]
    
    # Create subplots for R_200
    fig, ax = plt.subplots()
    
    # Plotting linear trend of corrected RM for Dwarf with limit = R_200
    grad = best_fit_line_plot(x=d_rm_R_200, ax=ax, 
                              y=dwrf_G_R_200_rm_val_s[i], 
                              y_err = bin_std_R_200,
                              get_grad=True)
    
    info = (d_rm_R_200, 
            dwrf_G_R_200_rm_val_s[i], 
            d_bg_R_200, 
            dwrf_G_R_200_bg_rm_val_s[i], 
            bin_means_R_200, 
            bin_std_R_200, 
            bin_edges_R_200, 
            Name, 
            R_200,
            d_bin_centers_R_200)
    
    #Storing data and gradient for recreating plots later in particualar order
    Dwarf_Plots.append((info, grad, dist_to_m31))
    
    plt.tight_layout()
    plt.close(fig)

    
#Sorting data based on smallest to largest gradient before starting to plot
Dwarf_Plots_sorted = sorted(Dwarf_Plots, key=lambda x: x[-1])

# In[Recreating the plot (Individual subplots)]:

num_plots = len(Dwarf_Plots_sorted)
num_cols = 5
num_rows = (num_plots + num_cols - 1) // num_cols  #Ceiling division to get enough rows

fig, axs = plt.subplots(num_rows, num_cols, figsize=(5 * num_rows, 34))
axs = axs.flatten() #Flattening 2D array of axes for easier indexing

handles, labels = [], []
for i, ((d_rm_R_200, 
          rm_vals, 
          d_bg_R_200, 
          bg_val, 
          bin_means, 
          bin_std, 
          bin_edges, 
          Name, 
          R_200,
          d_bin_centers), grad, D_to_m31) in enumerate(Dwarf_Plots_sorted):
    
    ax = axs[i]  # Select the appropriate axis for the current plot
   
    if i == 1:
        lab_1 = "RM" + r" [rad m$^{-2}$]" ; lab_2 = 'Background' + r" [rad m$^{-2}$]" ; lab_3 = r"$R_{200}$"+ " [kpc]"; lab_4 = "Mean" + r" [rad m$^{-2}$]"
    else:
        lab_1, lab_2, lab_3, lab_4= None, None, None, None
    
    # Plotting for R_200
    size=20
    ax.scatter(d_rm_R_200, rm_vals, marker='o', label=lab_1, s=size)
    ax.scatter(d_bg_R_200, bg_val, marker='o', label=lab_2, s=size)
    
    ax.plot([R_200, R_200], [150., -150.], linestyle='-.', color="purple", label=lab_3)
    
    ax.errorbar(d_bin_centers, bin_means, yerr=bin_std, fmt='k.-', label=lab_4)
    
    best_fit_line_plot(x=d_rm_R_200, y=rm_vals, y_err=bin_std, ax=ax)
    
    ax.text(0.5, 1.05, f'{Name} (m={grad[0]:.3g})', rotation='horizontal', fontsize=20, ha='center', va='bottom', transform=ax.transAxes)

    if i == 1:
        ax.legend(fontsize = 40, loc = 'upper center', bbox_to_anchor = (.5, 1.9), 
                        framealpha = 0, ncol = 4)
    
fig.text(-0.05, .5, 'RM'+ r" [rad m$^{-2}$]", va='center', rotation='vertical', fontsize=50) #Giving Common y_label
fig.text(.24, -0.03, 'Projected Separation from Dwarf Galaxy [kpc]', va='center', rotation='horizontal', fontsize=50) #Giving Common x_label

#Single legend outside of the plots
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.01, 1))

# Removing any empty subplots if num_plots < num_cols * num_rows
for j in range(num_plots, len(axs)):
    fig.delaxes(axs[j])

plt.tight_layout(h_pad=8, w_pad=-70)
plt.show()

# In[ Plotting all dwarf RM data all together before analysis]:
    
plt.figure(figsize=(10,5))
for i, ((d_rm_R_200,  #0
          rm_vals,    #1
          d_bg_R_200, #2
          bg_val,     #3
          bin_means,  #4
          bin_std,    #5 or -5
          bin_edges,  #6 or -4
          Name,       #7 or -3
          R_200,      #8 or -2
          d_bin_centers) #9 or -1
          , grad, D_to_m31) in enumerate(Dwarf_Plots_sorted):
        
        if i == 0:
            lab_1 = "RM" ; lab_2 = 'Background' ; lab_3 = r"$R_{200}$"
        else:
            lab_1, lab_2, lab_3 = None, None, None
    
        size=10
        
        d_rm_R_200 = d_rm_R_200/R_200
        
        plt.scatter(d_rm_R_200, rm_vals, marker='x', label=lab_1, s=size)

# ax.text(0.5, 1.05, f'm={grad:.3g})', rotation='horizontal', fontsize=35)plt.legend(fontsize = 40, loc = 'upper center', bbox_to_anchor = (.5, 1.9), 
#                 framealpha = 0, ncol = 3)
plt.xlabel(r'Normalized r$_{Dwarf G}$ [kpc]', 
           va='center', rotation='horizontal', fontsize=10,
           labelpad=15)
plt.ylabel(r'RM [rad m$^2$]', va='center', rotation='vertical', fontsize=15)

Dwarf_Plots_sorted_old = Dwarf_Plots_sorted
Dwarf_Plots_sorted = list(zip(*Dwarf_Plots))[0]
Dwarf_Plots_sorted = list(zip(*Dwarf_Plots_sorted))

#NORMALIZING the distance by R_200 (Same done for the bin_centres)
Dwarf_Plots_sorted[0] = [ d_proj/r_200 for d_proj, r_200 in list(zip(Dwarf_Plots_sorted[0], Dwarf_Plots_sorted[-2]))] 
Dwarf_Plots_sorted[-1] = [ d_bin_cen/r_200 for d_bin_cen, r_200 in list(zip(Dwarf_Plots_sorted[-1], Dwarf_Plots_sorted[-2]))] 

# Assuming you're trying to flatten the data for plotting
d_rm_R_200_flat = [x.value for sublist in list(Dwarf_Plots_sorted[0]) for x in sublist]
rm_vals_flat = [x for sublist in list(Dwarf_Plots_sorted[1]) for x in sublist]
d_bin_centers_flat = [x.value for sublist in list(Dwarf_Plots_sorted[-1]) for x in sublist]
bin_means_flat = [x for sublist in list(Dwarf_Plots_sorted[4]) for x in sublist]
bin_std_flat = [x for sublist in list(Dwarf_Plots_sorted[5]) for x in sublist]

#Sortign out data based on d_bin_centers: (getting chronological order)
combined = list(zip(d_bin_centers_flat, d_rm_R_200_flat, rm_vals_flat, bin_means_flat, bin_std_flat))
combined_sorted = sorted(combined, key=lambda x: x[0])
d_bin_centers_flat, d_rm_R_200_flat, rm_vals_flat, bin_means_flat, bin_std_flat = zip(*combined_sorted)

# Call your best fit line plot function
best_fit_line_plot(x=d_rm_R_200_flat, y=rm_vals_flat, y_err=bin_std_flat, ax=fig, just_plot=True)

# Using standard error of the mean for error bars within the R_200 limit
(bin_means_flat, 
 bin_std_flat, 
 _,#dont need the edges
d_bin_centers_flat) = get_bin_mean_std(X=d_rm_R_200_flat, 
                                       Y=rm_vals_flat, 
                                       bin_num=22, 
                                       get_edges=False, 
                                       norm=True)
plt.errorbar(d_bin_centers_flat, 
             bin_means_flat, 
             yerr=bin_std_flat, fmt='k.-', label="mean")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


def get_relative_size_for_scatter(max_size, length):
    # Larger size of scatter plot means closer to M31
    ratio_dist = dwarf_galaxies['D_dwarf']/max(dwarf_galaxies['D_dwarf'])
    ratio_distance = np.round((np.argsort(list(ratio_dist)) + 1)/length * max_size, 2)
    return ratio_distance

plt.figure(figsize=(15,8))
mca = generate_marker_color_pairs(len(Dwarf_Plots_sorted_old))
for i, ((d_rm_R_200,  #0
          rm_vals,    #1
          d_bg_R_200, #2
          bg_val,     #3
          bin_means,  #4
          bin_std,    #5 or -5
          bin_edges,  #6 or -4
          Name,       #7 or -3
          R_200,      #8 or -2
          d_bin_centers) #9 or -1
          , grad, D_to_m31) in enumerate(Dwarf_Plots_sorted_old):
        
        # if i == 0:
        #     lab_1 = "RM" ; lab_2 = 'Background' ; lab_3 = r"$R_{200}$"
        # else:
        #     lab_1, lab_2, lab_3 = None, None, None
    
        
        size= get_relative_size_for_scatter(max_size=160, length=len(Dwarf_Plots_sorted_old))[i]
        
        d_rm_R_200 = d_rm_R_200/R_200
        
        plt.scatter(d_rm_R_200, rm_vals, marker=mca[i][0],
                    color=mca[i][1], alpha=mca[i][2], label=Name, s=size)
        

# ax.text(0.5, 1.05, f'm={grad:.3g})', rotation='horizontal', fontsize=35)plt.legend(fontsize = 40, loc = 'upper center', bbox_to_anchor = (.5, 1.9), 
#                 framealpha = 0, ncol = 3)
plt.axvline(x=1, ymax=max(rm_vals), ymin=min(rm_vals), color="b", linestyle="-.", label=r"$R_{200}$ Normalized")
plt.xlabel(r'Normalized r$_{Dwarf G}$ [kpc]',
           va='center', rotation='horizontal', fontsize=15,
           labelpad=15)
plt.ylabel('RM [rad m$^2$]', va='center', rotation='vertical', fontsize=15)

# Call your best fit line plot function
best_fit_line_plot(x=d_rm_R_200_flat, y=rm_vals_flat, y_err=bin_std_flat , ax=fig, just_plot=True)

# Using standard error of the mean for error bars within the R_200 limit
(bin_means_flat,
 bin_std_flat,
 _,#dont need the edges
d_bin_centers_flat) = get_bin_mean_std(X=d_rm_R_200_flat,
                                       Y=rm_vals_flat,
                                       bin_num=22,
                                       get_edges=False, 
                                       norm=True)
plt.errorbar(d_bin_centers_flat, 
             bin_means_flat, 
             yerr=bin_std_flat, fmt='k.-', label="mean", capsize=2,
             zorder=1)

plt.grid(True)
plt.legend(loc='upper center', fontsize=11.3,
                 ncols=8, bbox_to_anchor=(0.5, 1.25), framealpha=0)

plt.tight_layout()
plt.show()
