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
rm, 
m31_condition,
bg_condition,
bg_pos_icrs,


#Importing functions
BG_correction, curr_dir_path
)

grid_res = 50 if args.grid_res is None else args.grid_res
rm_m31 = rm[m31_condition]
rm_bg = rm[bg_condition]
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
