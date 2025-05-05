# RM_vs_M33_radius.py

# Importing package alias'
from main import np, stats, u
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse

# Importing variables fromk main
from main import (
m33_sep_Rvir, rm_m33, bin_num, d_bg_m33, 
rm_bg_m33, d_rm_m33, R_vir_m33,  bg_pos_m33,
rm_pos_icrs_m33, bin_means_m33,

rm_pos_icrs_m33, 
rm, m33_condition,
bg_condition_m33,
bg_pos_icrs_m33,
rm_pos_gal_lat_m33, rm_pos_gal_lat_bg_m33,
args,

#importing functions from main
BG_correction, curr_dir_path, 
stats, 
create_annuli_binning_structure, 
apply_plot_attributes,
)

# Making sure it's untainted rm_m33 data
rm_m33 = rm[m33_condition]
rm_bg_m33 = rm[bg_condition_m33]

grid_res = 50 if args.grid_res is None else args.grid_res
rm_m33 = BG_correction(rm_pos_icrs_m33, rm_m33, bg_pos_icrs_m33, rm_bg_m33, grid_res=grid_res)
rm_bg_m33 = BG_correction(bg_pos_icrs_m33, rm_bg_m33, bg_pos_icrs_m33, rm_bg_m33, grid_res=grid_res)

# Bin RM within M33 CGM
bin_num = 10
bin_means, bin_edges, binnumber = stats.binned_statistic(
    d_rm_m33, rm_m33, statistic='mean', bins=bin_num)

bin_std, _, _ = stats.binned_statistic(
    d_rm_m33, rm_m33, statistic=lambda x: stats.sem(x), bins=bin_num)

bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

# Plot base RM vs r_M33
plt.figure(figsize=(8, 8))
plt.plot(d_rm_m33, rm_m33, 'x', markersize=3, label="RM within " + fr"R$_{{vir_{{M33}}}}$")
plt.plot(d_bg_m33, rm_bg_m33, 'x', markersize=3, label="Background")
plt.plot([R_vir_m33.value, R_vir_m33.value], [1000., -1000.], 'r--', linewidth=2, label="Rvir M33")

plt.errorbar(bin_centers, bin_means, yerr=bin_std, fmt='k.-', label='Binned Mean')# $\pm$ SEM')
plt.xlim(0., np.max(d_bg_m33.value))
plt.ylim(-R_vir_m33.value, R_vir_m33.value)
plt.xlabel(r'$r_{M33}$ [kpc]', fontsize=10)
plt.ylabel(r'RM [rad/$m^2$]', fontsize=10)
plt.grid(True, which="both")
plt.legend()

if args.save_plot:
    path = curr_dir_path() + "Results/"
    plt.savefig(f"{path}RM_vs_M33_radius_grid_res={grid_res}.png", dpi=600, bbox_inches="tight")
    print(f"RM_vs_M33_radius has been saved to {path}")
else:
    plt.show()

rm_pos_lat_all = np.concatenate([rm_pos_gal_lat_m33, rm_pos_gal_lat_bg_m33])
d_all = np.concatenate([d_rm_m33, d_bg_m33])
rm_all = np.concatenate([rm_m33, rm_bg_m33])

[distance_per_b_bin, rm_per_b_bin], _, bin_edges = create_annuli_binning_structure(
                                                                    bin_data=rm_pos_lat_all,
                                                                    data=(d_all, rm_all),
                                                                    bin_num=10,
                                                                    for_PA_or_b_plot=True
                                                                )

def plot_binned_gal_lat_m33(D, RM):
    gal_lat_lims = (-50, 10)
    num_bins = 12 + 1
    bin_edges = np.linspace(*gal_lat_lims, num_bins)
    bin_width = bin_edges[1] - bin_edges[0]

    cmap = plt.cm.get_cmap("tab20", num_bins)
    listed_cmap = mcolors.ListedColormap(cmap.colors[:num_bins])
    norm = mcolors.BoundaryNorm(boundaries=np.append(bin_edges, bin_edges[-1] + bin_width), ncolors=num_bins)

    fig, ax = plt.subplots(1, 1, figsize=(14, 10))

    for bin_idx in range(1, len(RM) + 1):
        if bin_idx in RM:
            color = listed_cmap(bin_idx - 1)
            ax.scatter(D[bin_idx], RM[bin_idx], marker="o", alpha=0.7, s=24,
                       edgecolors='none', color=color)

    sm = plt.cm.ScalarMappable(cmap=listed_cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_ticks(bin_edges)
    cbar.set_ticklabels([f"{lat:.1f}" + r"$^\circ$" for lat in bin_edges])
    cbar.set_label("Galactic Latitude [deg]", rotation=-90, labelpad=20, fontsize=15)

    ax.axvline(x=R_vir_m33.value, linestyle="--", color="k", label=fr"$R_{{vir, M33}} \approx$"+f"{R_vir_m33} -  Putman et al. (2018)")
    ax.errorbar(bin_centers, bin_means, yerr=bin_std, fmt='k.-', label='Binned data with error bars')

    apply_plot_attributes(xlim=(0,max(d_bg_m33.value)),
                        xlabel=fr'$r_{{M33}}$ [kpc]',
                        axis_lab_f_size=(13,13), #For axis label font sizes
                        bbox=(.5, 1.1), ncols=1, fontsize=14)

    if args.save_plot:
        path = curr_dir_path() + "Results/"
        fname = f"{path}RM_vs_M33_radius_colored_by_galactic_b.png"
        plt.savefig(fname, dpi=600, bbox_inches="tight")
        print(f"Galactic Latitude vs RM (M33) plot saved in {path}")
        plt.close()
    else:
        plt.show()

plot_binned_gal_lat_m33(distance_per_b_bin, rm_per_b_bin)