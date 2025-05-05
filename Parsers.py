import argparse
def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sig-limit', type=int, default=0, 
                        help='Sigma detection limitation on all RM errors. e.g. 5 sigma or 3 sigma threshholds (Default is 0 sigma)')
    
    parser.add_argument('--bg-corr', action='store_true', help='Background Correction ujsing Spline fit (Default=Fasle)')
    parser.add_argument('--b-limit', action='store_true', help='|b|> 5 degrees before Correction and Analysis')
    parser.add_argument('--b-limit-small', action='store_true', help='|b|> 5 degrees are only taking to account after background spline-fit corrections')
    parser.add_argument('--elliptic-CGM', action='store_true', help='Apply an Elliptic CGM to Andromeda. (not inclusive to Random R_vir)')
    parser.add_argument('--patch-num', type=int, default=int(1e4),
                        help='Declare number of patches (default is 10,000)')
    parser.add_argument('--original-plot', action='store_true', help='Plots the original plot from Honours (RM against projected distance of M31)')
    parser.add_argument('--pickle', action='store_true', help='Overwrites the pickled data where new random patches are produced and processed. (Save manually from jupyter ilifu to local)')
    parser.add_argument('--test-patches', action='store_true', help='Testing by showing patches on sphere as they get smaller')
    parser.add_argument('--CGM-dispersion', action='store_true', help='Also give dispersion plot of Rotation Measure within Halo and the average Random R_vir')
    parser.add_argument('--annuli-anal', action='store_true', help='Conducting annulus analysis with histograms of RANDOM patches in the sky')
    parser.add_argument('--annuli-video', action='store_true', help='Creating video of change in Rm per annulus for mean and median')
    parser.add_argument('--m31-annuli-anal', action='store_true', help='Conducting annulus analysis with histograms for M31 halo')
    parser.add_argument('--overplot', action='store_true', help='Enable overplotting; All radial annuli histograms on one plot. (only works if --annuli-anal is set)')
    parser.add_argument('--seaborn-hist', action='store_true', help='Enable overplotting; All radial annuli histograms on one plot. (only works if --annuli-anal is set)')
    parser.add_argument('--rm-vs-PA', action='store_true', help='plot RM as a function of Polar Angle in the anticlockwise direction with M31 as the reference frame')
    parser.add_argument('--m31-ks-test', action='store_true', help='Perform KS-Test Between random regions in the sky and that of M31')
    parser.add_argument('--rm-vs-gal-lat', action='store_true', help='Plotting RM against galactic latitude for M31 (inclusive of its background)')
    parser.add_argument('--rm-per-patch-hist', action='store_true', help='Histogram of how many RM points landed in each Random Virial Radius')
    parser.add_argument('--plot-random-cosdec', action='store_true', help='Plot histogram of number of random points per 0.1 bin in cos(|Dec|), split by positive and negative declinations')
    parser.add_argument('--cdf-anal', action='store_true', help='Making a Cumulative Density Plot for Random RM sources in the sky and m31')
    parser.add_argument('--print-ks-test', action='store_true', help="Prints the Multiple Tests conducted on the Random Sample from the Sky vs M31's CGM.")
    
    parser.add_argument('--cumulative', action='store_true', help='Convert  density plot to a cumulative probability map')
    parser.add_argument('--bw', type=float, help='Bandwidth for Cumulative/Probablity Density plots')
    
    parser.add_argument('--grid-res', type=int, help='Increasing Grid Resolution for BG data - for spline fitting')
    
    parser.add_argument('--save-plot', action='store_true', help='Saves the plot to Results folder rather than just showing the plot (To be used with --original-plot)')
    
    parser.add_argument('--mean', action='store_true', help='Choosing to only display results fro Mean')
    parser.add_argument('--median', action='store_true', help='Choosing to only display results fro Median')

    parser.add_argument('--smoothed', action='store_true', help='Resulting smoothed RM points within the Halo of M31')
    parser.add_argument('--others', action='store_true', help='Showing extra elements; dwarf galaxies, CHVCs etc')
    parser.add_argument('--dwarf-G', action='store_true', help='Showing only dwarf galaxies')
    parser.add_argument('--filled-485', action='store_true', help='Filling the contours for the 485 MHz')
    parser.add_argument('--_485mhz', action='store_true', help='Non-filled contours for the 485 MHz')
    parser.add_argument('--HI', action='store_true', help='Adding in the log_10(HI) contours')
    parser.add_argument('--scatter', action='store_true', help='Adding in scatter of RM and its BG on WCS axis plot')


    return parser


parser = build_parser()
args = parser.parse_args()

#Ensuring some arguments are only used when --annuli-anal is enabled
if not args.annuli_anal and not args.m31_annuli_anal:
    if args.overplot:
        parser.error("--overplot requires --annuli-anal or --m31-annuli-anal to be set.")
    if args.annuli_video:
        parser.error("--annuli-video requires --annuli-anal to be set.") #At leas tfor now it does... could expand to m31_annuli_anal....
elif args.overplot and args.annuli_video: #Only make a video if not overplotting (or superimposing plots)
    parser.error("--overplot cannot be done with --annuli-video")
if args.seaborn_hist and not args.overplot:
    parser.error("--seaborn-hist can only be used with --m31-annuli-anal")

if args.annuli_anal and args.m31_annuli_anal:
    parser.error("To lessen confusion please either use --annuli_anal or --m31-annuli-anal. Not Both")

if args.b_limit and args.b_limit_small:
    parser.error("Cant use both galactic latitude limits.")