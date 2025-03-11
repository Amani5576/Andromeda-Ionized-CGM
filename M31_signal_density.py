from scipy.stats import gaussian_kde

from M31_signal_vs_entire_sky import(
#importing alias'
np, plt,

#Importing Variables
all_d_bin_centers as bin_centers,
all_means as mn,
all_medians as md,
all_bin_stds as std,

#importing functions
plot_m31_stats, indiv_bg_corr
);

def fill_all_arrays_with_initial_0(arrays):
    arr_1 = arrays[0]

    Arr = []
    for idx, arr in enumerate(arrays):
        #Adding 0 (or 0-array) to ensure kde-plot starts from point (0,0)
        Arr.append(
        np.insert(np.array(arr), 0, np.zeros_like(arr_1), axis=0) 
        if len(arr) > 0 else np.array(arr) #checking if array is empty
        )
    return tuple(Arr)

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
def give_cbar_properties(cbar):
    # vmin, vmax = 0, max(Z_mean.max(), Z_med.max()) #colorbar range base don larges RM value in all datasets
    ticks = np.linspace(0, 1, 10)  #Generating 10 tick positions
    ticklabels = np.linspace(0, 1, len(ticks))  #Normalizing tick labels from 0 to 1
    
    cbar.set_ticks(ticks)  #Tick positions
    cbar.set_ticklabels([f"{tick:.1f}" for tick in ticklabels])  #Applying tick labels
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Density", rotation=-90, labelpad=20, fontsize=17)

# for i in [bin_centers, mn, md, std]: print(np.shape(i))
bin_centers, mn, md, std = fill_all_arrays_with_initial_0([bin_centers, mn, md, std])
# for i in [bin_centers, mn, md, std]: print(np.shape(i))

convert = lambda L: np.array(L).ravel()
x, std = convert(bin_centers), convert(std)
y_mean, y_med = convert(mn), convert(md)

x, y_mean, y_med, std = filter_nans(x, y_mean, y_med, std)

#executing background correction
y_mean = indiv_bg_corr(y_mean, bin_cent=x)
y_med = indiv_bg_corr(y_med, bin_cent=x)

#Generating a Kernel Density Estimation of the points.
xy_mean, xy_med = np.vstack([x, y_mean]), np.vstack([x, y_med])

#Creating grid for evaluation
grid_num = int(1e2)
x_grid = np.linspace(min(x), max(x), grid_num)
y_mean_grid, y_med_grid = np.linspace(min(y_mean), max(y_mean), grid_num), np.linspace(min(y_med), max(y_med), grid_num)
X, Y_mean = np.meshgrid(x_grid, y_mean_grid)
X, Y_med = np.meshgrid(x_grid, y_med_grid)

kde_mean  = gaussian_kde(xy_mean, bw_method=bw, weights=std)
kde_med = gaussian_kde(xy_med, bw_method=bw, weights=std)

Z_mean = kde_mean(np.vstack([X.ravel(), Y_mean.ravel()])).reshape(X.shape)
Z_med = kde_med(np.vstack([X.ravel(), Y_med.ravel()])).reshape(X.shape)
Z_max = max(Z_mean.max(), Z_med.max()) #To be used for normalisaztion
levels = 50

# Create subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

#Plots the data from intial starterkit (So nothing new here)
plot_m31_stats(axes[0]) 
plot_m31_stats(axes[1]) 


# PLOTTING THE MEAN (using contourf)
im1 = axes[0].contourf(X, Y_mean, Z_mean/Z_mean.max(), levels=levels, cmap="viridis")
# axes[0].scatter(x, y_mean, s=5, alpha=0.3, c='k')  # Original points
axes[0].set_title("Mean")
axes[0].set_xlabel("X-axis")
axes[0].set_ylabel("Y-axis")
axes[0].set_ylim(0, 85)
axes[0].set_xlim(0,300)

# cont_cbar_mean = fig.colorbar(im1, ax=axes[0], fraction=0.06, pad=0.01)
# give_cbar_properties(cont_cbar_mean)

# PLOTTING THE MEDIAN (using contourf)
im2 = axes[1].contourf(X, Y_med, Z_med/Z_med.max(), levels=levels, cmap="viridis")
# axes[1].scatter(x, y_med, s=5, alpha=0.3, c='k')  # Original points
axes[1].set_title("Median")
axes[1].set_xlabel("X-axis")
axes[1].set_ylim(0, 85)
axes[1].set_xlim(0,300)

#undoing repeat of y-axis ticks and ticklabels
axes[1].set_yticks([])
axes[1].set_yticklabels([])

cont_cbar_med = fig.colorbar(im2, ax=axes[1], fraction=0.06, pad=0.01)
give_cbar_properties(cont_cbar_med)

plt.tight_layout(w_pad=1)

plt.show()

