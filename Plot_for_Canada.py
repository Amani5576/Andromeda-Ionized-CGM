import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import argparse
import os

# ========== Argument Parsing ==========
parser = argparse.ArgumentParser(description="Display and optionally save side-by-side comparison of two images.")
parser.add_argument('--save-plot', action='store_true', help="Save the final image instead of just showing it.")
args = parser.parse_args()

# ========== Image Paths ==========
img1_path = 'Results/scatter_HI_contour_no_BG_corr.png'
img2_path = 'Results/2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png'

# ========== Load Images ==========
img1 = mpimg.imread(img1_path)
img2 = mpimg.imread(img2_path)

# ========== Create Subplots ==========
fig, axes = plt.subplots(1, 2, figsize=(12, 8))  # 1 row, 2 columns

# First subplot
axes[0].imshow(img1)
axes[0].axis('off')

# Second subplot
axes[1].imshow(img2)
axes[1].axis('off')

# ========== Layout and Save Logic ==========
plt.tight_layout()

if args.save_plot:
    results_dir = os.path.join(os.getcwd(), "Results")
    os.makedirs(results_dir, exist_ok=True)
    
    fname = "side_by_side_image_comparison.png"
    save_path = os.path.join(results_dir, fname)
    plt.savefig(save_path, dpi=600, bbox_inches="tight")
    print(f"✅ Saved comparison to: {save_path}")
    
    plt.close()
else:
    plt.show()
