from scipy.stats import gaussian_kde
from matplotlib.colors import LogNorm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--save-plot', action='store_true', help='Save the density plot of RM from random sources vs M31')
parser.add_argument('--cumulative', action='store_true', help='Convert the density plot to a cumulative probability map')
args = parser.parse_args()

from M31_signal_vs_entire_sky import (
    np, plt,
    all_d_bin_centers as bin_centers,
    all_means as mn,
    all_medians as md,
    all_bin_stds as std,
    number_of_patches,
    plot_m31_stats, curr_dir_path
)

print("Extracted pickled data for M31 background density plot.")

# def fill_all_arrays_with_initial_0(arrays):
#     arr_1 = arrays[0]

#     Arr = []
#     for idx, arr in enumerate(arrays):
#         #Adding 0 (or 0-array) to ensure kde-plot starts from point (0,0)
#         Arr.append(
#         np.insert(np.array(arr), 0, np.zeros_like(arr_1), axis=0) 
#         if len(arr) > 0 else np.array(arr) #checking if array is empty
#         )
#     return tuple(Arr)

def filter_nans(*arrays):
    """
    Filters out indices where any of the input arrays contain NaN values.
    Parameters:
        *arrays: Multiple NumPy arrays (must have the same length).
    Returns:
        A tuple of filtered NumPy arrays with no NaN-containing rows.
    """
    arrays = [np.array(arr) for arr in arrays]  # Ensure all inputs are NumPy arrays
    
    # Stack them to check for NaNs across all arrays
    stacked = np.column_stack(arrays)  # Shape: (N, num_arrays)
    
    # Create a mask where no NaNs are present in any row
    mask = ~np.isnan(stacked).any(axis=1)
    
    # Apply mask to each array to filter out bad row
    # Also making sure the x-axis starts from 0 by insering 0 as an initial for all arrays.
    return tuple(arr[mask] for arr in arrays)

global Z_mean, Z_med
def give_cbar_properties(cbar, Z, typ):

    ticks = np.linspace(0, 1, 11) if args.cumulative else np.linspace(Z.min(), Z.max(), 11)  # Use actual data range
    cbar.set_ticks(ticks)
    if args.cumulative:
        cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])
        label = f"{'Cumulative Probability' if args.cumulative else 'Density'} ({typ})"
    else:
        cbar.set_ticklabels([f"{tick*10**5:.0f}" for tick in ticks])
        label = f"{typ} Density " + r"[$10^{-5}\text{ kpc}^{-1} \cdot (\text{rad/m}^2)^{-1}$]"
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(label + "\n(10K Random " + r"R$_{vir}$)", rotation=-90, labelpad=24 if not args.cumulative else 27, fontsize=10)

# for i in [bin_centers, mn, md, std]: print(np.shape(i))
# bin_centers, mn, md, std = fill_all_arrays_with_initial_0([bin_centers, mn, md, std])
# for i in [bin_centers, mn, md, std]: print(np.shape(i))

convert = lambda L: np.array(L).ravel()
x, std = convert(bin_centers), convert(std)
y_mean, y_med = convert(mn), convert(md)

x, y_mean, y_med, std = filter_nans(x, y_mean, y_med, std)

#Generating a Kernel Density Estimation of the points.
xy_mean, xy_med = np.vstack([x, y_mean]), np.vstack([x, y_med])

#Creating grid for evaluation
grid_num = int(1e2)
x_grid = np.linspace(0, 300, grid_num)
y_mean_grid, y_med_grid = np.linspace(0, max(y_mean), grid_num), np.linspace(0, max(y_med), grid_num)
X, Y_mean = np.meshgrid(x_grid, y_mean_grid)
X, Y_med = np.meshgrid(x_grid, y_med_grid)

bw = 0.07 if not args.cumulative else 0.085
kde_mean  = gaussian_kde(xy_mean, bw_method=bw, weights=std)
kde_med = gaussian_kde(xy_med, bw_method=bw, weights=std)

Z_mean = kde_mean(np.vstack([X.ravel(), Y_mean.ravel()])).reshape(X.shape)
Z_med = kde_med(np.vstack([X.ravel(), Y_med.ravel()])).reshape(X.shape)

if args.cumulative:
    """
    Note: 
    Z_mean[::-1, :] --> reverses the order of the rows (flips the array vertically)"
    np.cumsum(..., axis=0) --> Cumulative sum along rows
    [::-1, :] --> Last part just flips back to original order of rows
    """
    Z_mean = np.cumsum(Z_mean[::-1, :], axis=0)[::-1, :]
    Z_med = np.cumsum(Z_med[::-1, :], axis=0)[::-1, :]
    Z_mean /= Z_mean.max()
    Z_med /= Z_med.max()
    
levels = 100

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

def plot_density(ax, X, Y, Z, title, xlabel, ylabel=None):
    """Plots a contour plot with colorbar and labels. (As well as M31 data)"""
    im = ax.contourf(X, Y, Z, levels=levels, cmap="viridis")
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_ylim(0, 200)
    ax.set_xlim(0, 300)

    plot_m31_stats(ax) #Line plots of Mean and Median for |RM| vs proj distance within CGM of M31
    ax.legend()

    cbar = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.01)
    give_cbar_properties(cbar, Z, title)

# Use the function for Mean and Median plots
plot_density(axes[0], X, Y_mean, Z_mean, "Mean", r'R$_{projected}$ [kpc]', "|RM| [rad " + r"m$^2$]")
plot_density(axes[1], X, Y_med, Z_med, "Median", r'R$_{projected}$ [kpc]')

plt.tight_layout(w_pad=2)
if args.save_plot:
    path = curr_dir_path() + "Results/"
    filename = f"{'Cumulative_' if args.cumulative else ''}Density_{bw=}.png"
    plt.savefig(f"{path}{filename}", dpi=600, bbox_inches="tight")
    print(f"{'Cumulative ' if args.cumulative else ''}Density Plot has been saved to {path}{filename}")
else:
    plt.show()

