from PIL import Image
import glob
import re
import os

# Base directory
base_dir = "./Results/Changing_CGM_4/"

# Get all histogram files and extract numbers
simple_hist_files = glob.glob(os.path.join(base_dir, "Simple_hist_*.png"))
hist_numbers = []

for fname in simple_hist_files:
    match = re.search(r"Simple_hist_(\d+)\.png", fname)
    if match:
        hist_numbers.append(int(match.group(1)))

if not hist_numbers:
    print("No histogram files found!")
    exit()

max_num = max(hist_numbers)

# Loop through and combine images efficiently
for i in range(1, max_num + 1):
    gauss_path = os.path.join(base_dir, f"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees_{i}.png")
    hist_path = os.path.join(base_dir, f"Simple_hist_{i}.png")
    out_path = os.path.join(base_dir, f"combined_{i}.png")

    if not os.path.exists(gauss_path) or not os.path.exists(hist_path):
        print(f"Skipping {i}: missing image.")
        continue

    try:
        # Open and convert to RGB just in case
        img1 = Image.open(gauss_path).convert("RGB")
        img2 = Image.open(hist_path).convert("RGB")

        # Resize to same height
        #if img1.size[1] != img2.size[1]:
        #    common_height = min(img1.size[1], img2.size[1])
        #    img1 = img1.resize((int(img1.size[0] * common_height / img1.size[1]), common_height))
        #    img2 = img2.resize((int(img2.size[0] * common_height / img2.size[1]), common_height))

        # Create new image and paste side by side
        total_width = img1.size[0] + img2.size[0]
        combined = Image.new("RGB", (total_width, img1.size[1]))
        combined.paste(img1, (0, 0))
        combined.paste(img2, (img1.size[0], 0))

        # Save and cleanup
        combined.save(out_path)
        img1.close()
        img2.close()
        combined.close()
        del img1, img2, combined

    except Exception as e:
        print(f"Error on pair {i}: {e}")
