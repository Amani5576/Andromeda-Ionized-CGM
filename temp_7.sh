#Scatter plot with HI contours on M31. No BG correctiona nd no Latitude cuts
python3 RM_smoothing_loop.py --save-plot --median --HI --scatter #--percent-cbar
mv ./Results/Temperature_Contours_within_CGM_of_M31.png ./Results/scatter_HI_contour_no_BG_corr.png

#Latitude cuts and BG spline fit correction on CGM of M31. 485 Mhz contours with Gaussian smoothing for visualisation and glactic coordinate grid superimposed on ICRS grid
python3 RM_smoothing_loop.py --bg-corr --b-limit-small --save-plot --median --smoothed --_485mhz --add-gal-grid  

python3 Plot_for_Canada.py --save-plot