import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def plot_grid(files, row_titles, col_titles, sup_title,
              folder='./Results/', output_path='./Results/combined_subplots.png', 
              dpi=200, figsize=(8, 14)):
    """
    Plot a grid of images (8 rows x 2 cols) with shared row and column titles,
    hide axes, and save to file.

    Parameters:
    - files: list of 16 filenames in row-major order
    - row_titles: list of 8 strings for each row label
    - col_titles: list of 2 strings for each column label
    - sup_title: string for the overall figure title
    - folder: path to image folder
    - output_path: where to save the combined figure
    - dpi: resolution for saving
    - figsize: figure size tuple
    """
    fig, axes = plt.subplots(4, 2, figsize=figsize)

    # Column titles
    for ax, ct in zip(axes[0], col_titles):
        ax.set_title(ct, pad=12, fontsize=12)

    # Fill grid
    for i, (row_axes, rt) in enumerate(zip(axes, row_titles)):
        # Row title
        row_axes[0].text(-0.1, 0.5, rt,
                         va='center', ha='right', rotation='vertical',
                         fontsize=12, transform=row_axes[0].transAxes)
        # Two images per row
        for ax, fname in zip(row_axes, files[2*i:2*i+2]):
            img = mpimg.imread(f"{folder}{fname}")
            ax.imshow(img)
            ax.axis('off')

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    fig.suptitle(sup_title, fontsize=14)
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved combined figure to {output_path}")


# Define file lists and titles for standard and elliptical CGM
standard_files = [
    'RM_vs_M31.png',
    'RM_vs_M31_colored_by_galactic_b.png',
    'RM_vs_M31_upto_40.png',
    'RM_vs_M31_colored_by_galactic_b_upto_40.png',
    'RM_vs_M31_upto_40_lat_cuts.png',
    'RM_vs_M31_colored_by_galactic_b_upto_40_lat_cuts.png',
    'RM_vs_M31_lat_cuts.png',
    'RM_vs_M31_colored_by_galactic_b_lat_cuts.png',
]
elliptical_files = [
    'RM_vs_M31_elliptical.png',
    'RM_vs_M31_colored_by_galactic_b_elliptical.png',
    'RM_vs_M31_upto_40_elliptical.png',
    'RM_vs_M31_colored_by_galactic_b_upto_40_elliptical.png',
    'RM_vs_M31_upto_40_lat_cuts_elliptical.png',
    'RM_vs_M31_colored_by_galactic_b_upto_40_lat_cuts_elliptical.png',
    'RM_vs_M31_lat_cuts_elliptical.png',
    'RM_vs_M31_colored_by_galactic_b_lat_cuts_elliptical.png',
]
row_titles = [
    "BG cutoff for spline grid ≤ 30°",
    "BG cut-off for spline grid ≤ 40°",
    "Lat Cuts ≤ 5° \n BG cut-off for spline grid ≤ 40°",
    "Lat Cuts ≤ 5° \n BG cut-off for spline grid ≤ 30°"
]
col_titles = ["No Galactic b coloring", "Colored by Galactic b-bins"]
sup_titles = ["Circular CGM", "Elliptic CGM" ]
#Plot and save both grids seperately
plot_grid(standard_files, row_titles, col_titles, sup_title = sup_titles[0],
            output_path='./Results/Upto_40_deg_cut_off_testing.png')
plt.close()
plot_grid(elliptical_files, row_titles, col_titles, sup_title = sup_titles[1],
            output_path='./Results/Upto_40_deg_cut_off_testing_elliptical.png')