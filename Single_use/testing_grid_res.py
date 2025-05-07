# def BG_correction_checking(rm_coords, rm_values, bg_coords, bg_values, **kw):
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from astropy.io import fits
#     from astropy.wcs import WCS
#     from scipy.interpolate import griddata
#     from matplotlib.colors import Normalize

#     curr_dir_path = lambda: '/home/amani/Documents/MASTERS_UCT/GitHub_stuff/Andromeda-Ionized-CGM/'

#     def get_wcs(filename):
#         hdu = fits.open(filename)[0]
#         return WCS(hdu.header)

#     def plot_BG_density(rm_coords, bg_coords, bg_values, fits_file, **kw):
#         """
#         Plots:
#         - Background RM density using interpolation.
#         - Background RM positions as scatter points.
#         - Uses WCS projection from FITS file.

#         Parameters:
#         -----------
#         rm_coords : SkyCoord
#             Sky coordinates of RM measurements.
#         bg_coords : SkyCoord
#             Sky coordinates of background RM measurements.
#         bg_values : np.ndarray
#             RM values at background positions.
#         fits_file : str
#             Path to the FITS file for WCS projection.

#         Returns:
#         --------
#         None
#         """
#         # Convert coordinates to degrees
#         x_bg = bg_coords.ra.deg
#         y_bg = bg_coords.dec.deg

#         # Invert RA greater than 180 to negative RA + 180
#         x_bg = np.where(x_bg > 180, x_bg - 360, x_bg)

#         # Define grid
#         grid_res = int(kw.get("grid_res", 50))
#         # print(f"{grid_res=}")
#         x_grid = np.linspace(x_bg.min(), x_bg.max(), grid_res)
#         y_grid = np.linspace(y_bg.min(), y_bg.max(), grid_res)
#         X_grid, Y_grid = np.meshgrid(x_grid, y_grid)

#         # Interpolate using cubic method
#         bg_points = np.column_stack((x_bg, y_bg))
#         grid_points = np.column_stack((X_grid.ravel(), Y_grid.ravel()))
#         bg_grid = griddata(bg_points, bg_values, grid_points, method='cubic')
#         bg_grid = bg_grid.reshape(X_grid.shape)

#         # Handle NaNs with nearest-neighbor interpolation
#         if np.isnan(bg_grid).any():
#             bg_grid = griddata(bg_points, bg_values, grid_points, method='nearest')
#             bg_grid = bg_grid.reshape(X_grid.shape)

#         # Plot
#         fig, ax = plt.subplots(figsize=(10, 10))

#         # Density plot (background)
#         norm = Normalize(vmin=np.nanmin(bg_grid), vmax=np.nanmax(bg_grid))
#         im = ax.imshow(bg_grid, extent=[x_bg.min(), x_bg.max(), y_bg.min(), y_bg.max()], 
#                     origin='lower', cmap='viridis', norm=norm)

#         # Scatter plot of background RM points (foreground)
#         ax.scatter(x_bg, y_bg, color='red', s=0.4, label="RM Background Points", alpha=0.6)

#         # Set RA limits from -180 to 180
#         ax.set_xlim(-40, 60)

#         # Set the plot limits to the range of background coordinates
#         ax.set_ylim(y_bg.min(), y_bg.max())

#         # Axis labels & title
#         ax.set_xlabel("RA (deg)")
#         ax.set_ylabel("Dec (deg)")
#         ax.set_title(f"Grid resolution: {grid_res}")

#         # Add legend
#         ax.legend(bbox_to_anchor=(0.35, 1.3))

#         # Add colorbar
#         plt.colorbar(im, ax=ax, label="Background RM Density", shrink=.5)

#         # Save the figure
#         path = curr_dir_path() + "Results/BG_Res_Test/"
#         plt.savefig(f"{path}BG_grid_res={grid_res}.png", dpi=600, bbox_inches="tight")


#     plot_BG_density(rm_coords, bg_coords, bg_values, "LGSNLGSR.SQLGBB.FITS", grid_res=kw["grid_res"])

# BG_correction_checking(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg, grid_res = 50)

# for N in np.linspace(1,50, 51):
#     BG_correction_checking(rm_pos_icrs, rm_m31, bg_pos_icrs, rm_bg, grid_res = N)

    # for N in np.linspace(10,1000,80):
#     BG_correction_checking(bg_pos_icrs, rm_bg, grid_res = N)

# def BG_correction_checking(bg_coords, bg_values, **kw):
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from astropy.io import fits
#     from astropy.wcs import WCS
#     from scipy.interpolate import griddata
#     from matplotlib.colors import Normalize

#     curr_dir_path = lambda: '/home/amani/Documents/MASTERS_UCT/GitHub_stuff/Andromeda-Ionized-CGM/'

#     def get_wcs(filename):
#         hdu = fits.open(filename)[0]
#         return WCS(hdu.header)

#     def plot_BG_density(bg_coords, bg_values, fits_file, **kw):
#         # Convert coordinates to degrees
#         x_bg = bg_coords.ra.deg
#         y_bg = bg_coords.dec.deg

#         # Define grid
#         grid_res = int(kw.get("grid_res", 50))
#         print(f"{grid_res=}")
#         x_grid = np.linspace(x_bg.min(), x_bg.max(), grid_res)
#         y_grid = np.linspace(y_bg.min(), y_bg.max(), grid_res)
#         X_grid, Y_grid = np.meshgrid(x_grid, y_grid)

#         # Interpolate using cubic method
#         bg_points = np.column_stack((x_bg, y_bg))
#         grid_points = np.column_stack((X_grid.ravel(), Y_grid.ravel()))
#         bg_grid = griddata(bg_points, bg_values, grid_points, method='cubic')
#         bg_grid = bg_grid.reshape(X_grid.shape)

#         # Handle NaNs with nearest-neighbor interpolation
#         if np.isnan(bg_grid).any():
#             bg_grid = griddata(bg_points, bg_values, grid_points, method='nearest')
#             bg_grid = bg_grid.reshape(X_grid.shape)

#         # WCS projection from FITS file
#         wcs = get_wcs(fits_file)

#         # Plot
#         fig = plt.figure(figsize=(16, 10))
#         ax = fig.add_subplot(111, projection=get_wcs("LGSNLGSR.SQLGBB.FITS"), slices=('x', 'y', 0, 0))

#         # Density plot (background)
#         im = ax.imshow(bg_grid, origin='lower', cmap='viridis', transform=ax.get_transform('world'))

#         # Scatter plot of background RM points (foreground)
#         ax.scatter(x_bg, y_bg, color='red', s=0.4, label="RM Background Points", alpha=0.6, transform=ax.get_transform('world'))

#         # Axis labels & titles
#         ax.set_xlabel("RA (deg)")
#         ax.set_ylabel("Dec (deg)")
#         ax.set_title(f"{grid_res=}")
#         ax.legend()

#         # Add colorbar
#         plt.colorbar(im, ax=ax, label="Background RM Density")

#         # Save the figure
#         path = curr_dir_path() + "Results/BG_Res_Test/"
#         plt.savefig(f"{path}wcs_attempt.png", dpi=600, bbox_inches="tight")


#     plot_BG_density(bg_coords, bg_values, "LGSNLGSR.SQLGBB.FITS", grid_res=kw["grid_res"])
# BG_correction_checking(bg_pos_icrs, rm_bg, grid_res = 10)

# # for N in np.linspace(10,1000,80):
# #     BG_correction_checking(bg_pos_icrs, rm_bg, grid_res = N)
