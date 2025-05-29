#Importing package alias'
from main import np, stats
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse

"""NOTE: 
rotation measures landing in CGM are dependent on whether galactic binning
is commented out or not in main.py
"""
#Importing vairables
from main import (
bin_num, d_bg, R_vir, 
d_bin_centers,
bin_means, d_m31, L_m31,
args,

#Importing functions
BG_correction, curr_dir_path, 
create_annuli_binning_structure, 
apply_plot_attributes,
get_projected_d,
get_data_from_catalogue,
get_CGM_and_BG_masks,
apply_CGM_and_BG_masks
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
    elliptic_CGM_bg=args.elliptic_CGM_bg
)

# 3)Apply those masks to slice out CGM & BG subsets
(bg_pos, bg_pos_icrs, rm_pos, rm_pos_icrs, rm_pos_gal_lat, 
 rm_pos_gal_lat_bg, rm_bg, m31_sep_bg, err_bg, rm_m31, 
 m31_sep_Rvir, err_m31) = apply_CGM_and_BG_masks(
    rm_m31_coord, eq_pos, position, rm, rm_err, 
    m31_sep, m31_condition, bg_condition
)

grid_res = 50 if args.grid_res is None else args.grid_res
rm_m31 = BG_correction(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg, grid_res = grid_res)
rm_bg = BG_correction(bg_pos_icrs, rm_bg, bg_pos_icrs, rm_bg, grid_res = grid_res)
d_rm = get_projected_d(m31_sep_Rvir, d_m31) #Is only for within R_vir
d_bg = get_projected_d(m31_sep_bg, d_m31) #Is only for within R_vir

# print(f"{len(rm_m31)=}"); import sys; sys.exit()

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


#New version below gives color to RM's pertainting to a particular bin

d_s = np.concatenate([d_rm,d_bg])
rm_s = np.concatenate([rm_m31,rm_bg])
rm_pos_gal_lats = np.concatenate([rm_pos_gal_lat, rm_pos_gal_lat_bg]) #Inclusive of its background (if used...)

# print(f"{len(d_s)=}")
# print(f"{len(rm_s)=}")
# print(f"{len(rm_pos_gal_lats)=}") ; import sys; sys.exit()
    
[distance_per_b_bin, rm_per_b_bin], _, bin_edges = create_annuli_binning_structure(bin_data=rm_pos_gal_lats, 
                                                                                        data=(d_s, rm_s), 
                                                                                        bin_num=10, 
                                                                                        for_PA_or_b_plot=True #Work similar to how its been used before
                                                                                        )

def plot_binned_gal_lat(D, RM):
    """
    Raw RM against distance colored by galactic latitude (binned),
    with discrete colormap (10 galactic-b-bins) and mean errorbars.
    """
    gal_lat_lims = (-50,10)
    num_bins = 12 +1 # + 1 cause of using linspace

    # Distance bins
    bin_edges = np.linspace(*gal_lat_lims, num_bins)
    bin_width = bin_edges[1] - bin_edges[0]
    bin_centers = bin_edges - bin_width / 2

    # Discrete colors using ListedColormap
    cmap = plt.cm.get_cmap("tab20", num_bins)  # or any qualitative map
    listed_cmap = mcolors.ListedColormap(cmap.colors[:num_bins])
    norm = mcolors.BoundaryNorm(boundaries=np.append(bin_edges, bin_edges[-1]+bin_width), ncolors=num_bins)

    fig, ax = plt.subplots(1, 1, figsize=(14, 10))

    # # Background points
    # ax.scatter(d_bg, rm_bg, marker="x", alpha=1, s=20, c="k", label="Background")

    # Plot each galactic latitude bin with discrete color
    for bin_idx in range(1, len(RM) + 1):
        if bin_idx in RM:
            color = listed_cmap(bin_idx - 1)  # offset index for colormap
            ax.scatter(D[bin_idx], RM[bin_idx], marker="o", alpha=0.7, s=24,
                       edgecolors='none', color=color)

    # Colorbar for Galactic latitude bins
    sm = plt.cm.ScalarMappable(cmap=listed_cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_ticks(bin_edges)
    cbar.set_ticklabels([f"{lat:.1f}"+ r"$^{\circ}$" for lat in bin_edges])
    cbar.set_label("Galactic Latitude [deg]", rotation=-90, labelpad=20)

    # Vertical line at M31's latitude
    ax.axvline(x=R_vir.value, linestyle="--", color="k", label=r"$b_{M31}$")

    # Plot error bars
    ax.errorbar(d_bin_centers, bin_means, yerr=bin_std, fmt='k.-', label='Binned data with error bars')

    # Apply styling and layout
    apply_plot_attributes(xlim=(0, 450), ylim=(-300, 300),
                          xlabel=r'$r_{M31}$ [kpc]',
                          push_title_up=1.1, ncols=1)

    # Save or show
    if args.save_plot:
        path = curr_dir_path() + "Results/"
        fname = f"{path}RM_vs_M31_radius_colored_by_galactic_b.png"
        plt.savefig(fname, dpi=600, bbox_inches="tight")
        print(f"Galactic Latitude vs RM plot saved in {path}")
        plt.close()
    else:
        plt.show()

plot_binned_gal_lat(distance_per_b_bin, rm_per_b_bin)

