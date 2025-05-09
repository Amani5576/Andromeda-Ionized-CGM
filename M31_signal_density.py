from scipy.interpolate import RegularGridInterpolator
from scipy.stats import gaussian_kde
from matplotlib.colors import LogNorm

from M31_signal_vs_entire_sky import (
    np, plt,
    all_d_bin_centers as bin_centers, #For all Random Rvir
    all_means as mn,
    all_medians as md,
    all_bin_stds as std,
    bin_num_from_main,
    number_of_patches,
    d_bin_centers_m31, #Just for M31 alone
    rm_m31,
    args,
    plot_m31_stats, curr_dir_path
)

from main import (
    bin_means as bin_means_m31, 
    bin_med as bin_med_m31,
    get_discrete_colors
)
print("Extracted pickled data for M31 background density plot.")

def annotate(ax, x_vals, y_vals, vals, color):
    """
    Annotates each point with its fraction above M31 value or 
    jsut M31 value if not cumulative.
    """
    
    for x, y, f in zip(x_vals, y_vals, vals):
        x = x.value # Can only apply 'add' function to dimensionless quantities
        ax.text(x-5, y + 13, f"{f:.3f}", fontsize=8, color=color, alpha=1,
                bbox=dict(
                facecolor="white",   # Background color of the text box
                edgecolor="none",    # No border
                boxstyle="round,pad=0.2",  # Rounded box with slight padding
                alpha=0.7             # Transparency for blending
            ))


def cumulative_fraction_above_m31(X, Y, Z_cum, m31_x, m31_y):
    """
    For each radius (m31_x), returns  fraction of KDE above  M31 |RM| value (m31_y),
    using  cumulative KDE Z_cum.

    Assumes Z_cum has been computed with cumulative sum from high to low y (RM) values.

    Parameters:
    -----------
    X, Y : 2D meshgrid of x and y
    Z_cum : 2D cumulative KDE array (same shape as X and Y)
    m31_x, m31_y : arrays of bin centers and corresponding M31 |RM| values (mean or median)

    Returns:
    --------
    fractions : np.array of cumulative fractions above M31 |RM| value at each radius
    """

    x_grid = X[0, :]
    y_grid = Y[:, 0]

    # Ensure M31 arrays are within bounds
    m31_x = np.clip(m31_x.value, x_grid.min(), x_grid.max())
    m31_y = np.clip(m31_y, y_grid.min(), y_grid.max())

    interp_func = RegularGridInterpolator((y_grid, x_grid), Z_cum)

    points = np.vstack([m31_y, m31_x]).T
    fractions = interp_func(points)

    return fractions

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
    Filters out indices where any of  input arrays contain NaN values.
    Parameters:
        *arrays: Multiple NumPy arrays (must have  same length).
    Returns:
        A tuple of filtered NumPy arrays with no NaN-containing rows.
    """
    arrays = [np.array(arr) for arr in arrays]  # Ensure all inputs are NumPy arrays
    
    # Stack m to check for NaNs across all arrays
    stacked = np.column_stack(arrays)  # Shape: (N, num_arrays)
    
    # Create a mask where no NaNs are present in any row
    mask = ~np.isnan(stacked).any(axis=1)
    
    # Apply mask to each array to filter out bad row
    # Also making sure  x-axis starts from 0 by insering 0 as an initial for all arrays.
    return tuple(arr[mask] for arr in arrays)

global Z_mean, Z_med
def give_cbar_properties(cbar, Z, typ=False):

    ticks = np.linspace(0, 1, 11) if args.cumulative else np.linspace(Z.min(), Z.max(), 11)  # Use actual data range
    cbar.set_ticks(ticks)
    if args.cumulative:
        typ = f"({typ})" if typ else ""
        cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])
        label = f"{'Cumulative Probability' if args.cumulative else 'Density'} {typ}"
    else:
        cbar.set_ticklabels([f"{tick*10**5:.1f}" for tick in ticks])
        label = f"{typ if typ else ''} Density " + r"[$10^{-5}\text{ kpc}^{-1} \cdot (\text{rad/m}^2)^{-1}$]"
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(label + f"\n({number_of_patches} Random " + r"R$_{vir}$)", rotation=-90, labelpad=24 if not args.cumulative else 27, fontsize=10)

# for i in [bin_centers, mn, md, std]: print(np.shape(i))
# bin_centers, mn, md, std = fill_all_arrays_with_initial_0([bin_centers, mn, md, std])
# for i in [bin_centers, mn, md, std]: print(np.shape(i))

convert = lambda L: np.array(L).ravel()
x, std = convert(bin_centers), convert(std)
y_mean, y_med = convert(mn), convert(md)

x, y_mean, y_med, std = filter_nans(x, y_mean, y_med, std)

#Generating a Kernel Density Estimation of  points.
xy_mean, xy_med = np.vstack([x, y_mean]), np.vstack([x, y_med])

#Creating grid for evaluation
grid_num = int(1e2)
x_grid = np.linspace(0, 300, grid_num)
y_upper = 220 if args.cumulative else 300 #truly decides the y-limit of the grid and data points
y_mean_grid, y_med_grid = (np.linspace(0, 
                                       y_upper, #max(y_mean), 
                                       grid_num), 
                           np.linspace(0, 
                                       y_upper, #max(y_med), 
                                       grid_num))
X, Y_mean = np.meshgrid(x_grid, y_mean_grid)
X, Y_med = np.meshgrid(x_grid, y_med_grid)

#For condition |b| > 5 --> if args.bw else (0.15 if args.cumulative else 0.1)
if bin_num_from_main == 10:# <1687: #|b| > 5 limit imposed
    bw = args.bw if args.bw else 0.15 #this is for 10 radial bins
else:
    bw = args.bw if args.bw else 0.07 #this is for 22 radial bins

kde_mean  = gaussian_kde(xy_mean, bw_method=bw, weights=std)
kde_med = gaussian_kde(xy_med, bw_method=bw, weights=std)

Z_mean = kde_mean(np.vstack([X.ravel(), Y_mean.ravel()])).reshape(X.shape)
Z_med = kde_med(np.vstack([X.ravel(), Y_med.ravel()])).reshape(X.shape)

if args.cumulative:
    """
    Note: 
    Z_mean[::-1, :] --> reverses  order of  rows (flips  array vertically)"
    np.cumsum(..., axis=0) --> Cumulative sum along rows
    [::-1, :] --> Last part just flips back to original order of rows
    """
    Z_mean = np.cumsum(Z_mean[::-1, :], axis=0)[::-1, :]
    Z_med = np.cumsum(Z_med[::-1, :], axis=0)[::-1, :]
    Z_mean /= np.max(Z_mean, axis=0)
    Z_med /= np.max(Z_med, axis=0)
    
    #Getting fraction of KDE above M31's value
    frac_above_mean = cumulative_fraction_above_m31(X, Y_mean, Z_mean, d_bin_centers_m31, np.abs(bin_means_m31))
    frac_above_median = cumulative_fraction_above_m31(X, Y_med, Z_med, d_bin_centers_m31, np.abs(bin_med_m31))

    # print(f"{frac_above_mean=}")
    # print(f"{frac_above_median=}")

levels = 100
levels_2 = 20 #For additional distinct lined contour
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

def plot_density(ax, X, Y, Z, title, xlabel, ylabel=None, show_legend=False):
    """Plots a contour plot with colorbar and labels. (As well as M31 data)"""
    
    if args.cumulative:
        cmap_discrete = get_discrete_colors(cmap_name="gray_r", 
                                           n_bins=levels_2, 
                                           data_limits=(0,1))
    else: cmap_discrete = "viridis"
    im = ax.contourf(X, Y, Z, levels=levels, cmap=cmap_discrete)

    #Additional distinct contour:
    ax.contour(X, Y, Z, levels=levels_2, cmap="gray" if args.cumulative else "gray_r", linewidths=.5)

    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.set_ylim(0, y_upper)
    ax.set_xlim(0, 300)

    plot_m31_stats(ax, title=title, color="r" if args.cumulative else "white")  # Line plots of Mean and Median for |RM| vs proj distance within CGM of M31
    
    # if show_legend:
    #     fig.legend(framealpha=0, bbox_to_anchor=(.6, 1.04), ncols =2)
    ax.legend(facecolor="white", framealpha=.4)
    
    if title == "Median":
        cbar = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.01)
        give_cbar_properties(cbar, Z)

    if args.cumulative:
        #Adding in text for fraction of KDE above the M31
        if title == "Mean":
            annotate(ax, d_bin_centers_m31, np.abs(bin_means_m31), frac_above_mean, "k")
        elif title == "Median":
            annotate(ax, d_bin_centers_m31, np.abs(bin_med_m31), frac_above_median, "k")
    else: #Just the Density
        if title == "Mean":
            annotate(ax, d_bin_centers_m31, np.abs(bin_means_m31), np.abs(bin_means_m31), "k")
        elif title == "Median":
            annotate(ax, d_bin_centers_m31, np.abs(bin_med_m31), np.abs(bin_med_m31), "k")


if args.mean and not args.median:
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_density(ax, X, Y_mean, Z_mean, "Mean", r'R$_{projected}$ [kpc]', "|RM| [rad " + r"m$^{-2}$]", show_legend=True)

elif args.median and not args.mean:
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_density(ax, X, Y_med, Z_med, "Median", r'R$_{projected}$ [kpc]', "|RM| [rad " + r"m$^{-2}$]", show_legend=True)

else:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    plot_density(axes[0], X, Y_mean, Z_mean, "Mean", r'R$_{projected}$ [kpc]', "|RM| [rad " + r"m$^{-2}$]", show_legend=True)
    plot_density(axes[1], X, Y_med, Z_med, "Median", r'R$_{projected}$ [kpc]')

plt.tight_layout(w_pad=2)
if args.save_plot:
    path = curr_dir_path() + "Results/"
    filename = f"""{'Cumulative_' if args.cumulative else ''}Density_{bw=}_{args.patch_num}_patches_{
        'mean' if args.mean else ('median' if args.median else "")
        }.png"""
    plt.savefig(f"{path}{filename}", dpi=600, bbox_inches="tight")
    print(f"{'Cumulative ' if args.cumulative else ''}Density Plot has been saved to {path}{filename}")
else:
    plt.show()
