import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import argparse
import os

def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoothed', action='store_true')
    parser.add_argument('--CGM-dispersion', action='store_true')  # Optional, you can drop this if unused
    parser.add_argument('--temp', action='store_true')   # Optional, you can drop this if unused
    return parser

def get_variant_suffix(args):
    if args.smoothed:
        return '_not_smootheD', 'debugging_assesment_not_smoothed.png'
    elif args.temp:
        return '_smoothed', 'debugging_assesment_smoothed.png'
    elif args.CGM_dispersion:
        return '_disp', 'debugging_assesment_dispersion.png'
    else:
        return '', 'debugging_assesment.png' #Cumulative Frequency

def main():
    parser = build_parser()
    args = parser.parse_args()
    suffix, outname = get_variant_suffix(args)

    base_path = './Results/'

    # 3 rows × 2 columns: (Circular, Elliptical)
    variants = [
        ('nothing_no_bg',            'elliptical-CGM_no_bg',           'No lat cuts, No BG Corr'),
        ('nothing',                 'elliptical-CGM',                  'No lat cuts, With BG Corr'),
        ('b-limit-small',          'b-limit-small-elliptical-CGM',    'Lat cuts, With BG Corr'),
    ]

    fig, axs = plt.subplots(3, 2, figsize=(12, 15))

    for i, (circ_prefix, ellip_prefix, row_title) in enumerate(variants):
        for j, prefix in enumerate([circ_prefix, ellip_prefix]):
            filename = os.path.join(base_path, prefix + suffix + '.png')
            ax = axs[i][j]
            if os.path.exists(filename):
                img = mpimg.imread(filename)
                ax.imshow(img)
                ax.set_title(f"{'Circular' if j == 0 else 'Elliptical'} CGM\n{row_title}")
            else:
                print(f"Warning: Missing file {filename}")
                ax.text(0.5, 0.5, 'Image not found', ha='center', va='center')
            ax.axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(base_path, outname), dpi=600)
    plt.close()

if __name__ == '__main__':
    main()
