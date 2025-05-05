#Importing package alias'
from main import np, stats
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--grid-res', type=int, help='Increasing Grid Resolution for BG data - for spline fitting')
    parser.add_argument('--save-plot', action='store_true', help='Saves the plot to Results folder rather than just showing the plot')

    args = parser.parse_args()


#Importing vairables
from main import (
m31_sep_Rvir, rm_m31, bin_num, d_bg, 
rm_bg, d_rm, R_vir, d_bin_centers, bin_means,

rm_pos_icrs, 
rm, m31_condition,
bg_condition,
bg_pos_icrs,
rm_pos_gal_lat, rm_pos_gal_lat_bg,
bin_num as bin_num_from_main,
m31_lat, #Galactic Latitude of M31 in degrees

#Importing functions
BG_correction, curr_dir_path, create_annuli_binning_structure, apply_plot_attributes
)

#Makiing sure its untainted rm_m31 data
rm_m31 = rm[m31_condition] 
rm_bg = rm[bg_condition]

grid_res = 50 if args.grid_res is None else args.grid_res
rm_m31 = BG_correction(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg, grid_res = grid_res)
rm_bg = BG_correction(bg_pos_icrs, rm_bg, bg_pos_icrs, rm_bg, grid_res = grid_res)


#Using standard error of mean for error bars: 
bin_std, bin_edges, binnumber = stats.binned_statistic(m31_sep_Rvir, 
    rm_m31, 
    statistic = lambda rm_m31: stats.sem(rm_m31), #Standard Error ofMean
    bins = bin_num)  

plt.figure(figsize=(8,8))
plt.plot(d_rm, rm_m31, 'x', markersize=3)
plt.plot(d_bg, rm_bg, 'x', markersize=3)
plt.plot([R_vir.value, R_vir.value], [1000., -1000.], 'r--', linewidth=2)
plt.xlim(0., np.max(d_bg.value))
plt.ylim(-R_vir.value, R_vir.value) #within Virial radius
plt.errorbar(d_bin_centers, bin_means, yerr = bin_std, fmt = 'k.-')
plt.xlabel(r'$r_{M31}$ [kpc]', fontsize=10)
plt.ylabel(r'RM [rad/$m^2$]', fontsize=10)
plt.grid(True, which="both")

if args.save_plot: #Save individual plots
    path = curr_dir_path() + "Results/"
    plt.savefig(f"{path}RM_vs_M31_radius_grid_res={grid_res}.png", dpi=600, bbox_inches="tight")
    print(f"RM_vs_M31_radius has been saved to {path}")
else:
    plt.show()

#________________________________________________________________________
#New version below gives color to RM's pertainting to a particular bin
#________________________________________________________________________

[distance_per_b_bin, rm_per_b_bin], _, bin_edges = create_annuli_binning_structure(bin_data=rm_pos_gal_lat, 
                                                                                        data=(d_rm, rm_m31), 
                                                                                        bin_num=bin_num_from_main+1, 
                                                                                        for_PA_or_b_plot=True)
def RM_vs_Dist_gal_lat(distances, rms, gal_latitudes, save_plot=False):
    """
    Scatter plot of RM vs distance from M31, colored by Galactic Latitude (b).
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Scatter background RMs (if available)
    ax.scatter(d_bg, rm_bg, marker="x", alpha=1, s=14, c="k", label="Background")

    # Normalize and colormap
    norm = plt.Normalize(vmin=-51, vmax=10)  # Galactic latitude spans [-90, 90]
    cmap = plt.cm.tab20

    # Main scatter: RM vs distance, colored by galactic latitude
    sc = ax.scatter(distances, rms, c=gal_latitudes, cmap=cmap, norm=norm,
                    s=24, alpha=0.9, edgecolors='none')
    
    plt.errorbar(d_bin_centers, bin_means, yerr = bin_std, fmt = 'k.-')

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Galactic Latitude [deg]", rotation=-90, labelpad=20)

    # Optional extras
    ax.axvline(x=R_vir.value, linestyle="--", color="r", label="M31 $R_{vir}$")

    ax.set_xlabel(r'$r_{M31}$ [kpc]')
    ax.set_ylabel(r'RM [rad/$m^2$]')
    # ax.grid(True, which="both")

    apply_plot_attributes(leg=True, bbox = (0.5, 1.1), xlim=(0,450), ylim=(-300,300))

    if save_plot:
        path = curr_dir_path() + "Results/"
        fname = f"{path}RM_vs_M31_distance_colored_by_gal_lat.png"
        plt.savefig(fname, dpi=600, bbox_inches="tight")
        print(f"Saved to {fname}")
        plt.close()
    else:
        plt.show()


RM_vs_Dist_gal_lat(d_rm, rm_m31, rm_pos_gal_lat)#save_plot=True)


