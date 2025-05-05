#No Latitude cuts or BG correction on Circular CGM
python3 RM_smoothing.py --save-plot --_485mhz --HI --smoothed #--scatter
mv Results/"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png" Results/nothing_no_bg_smoothed.png

#No Latitude cuts but with BG correction on Circular CGM
python3 RM_smoothing.py --bg-corr --save-plot --_485mhz --HI --smoothed #--scatter
mv Results/"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png" Results/nothing_smoothed.png

#Latitude cuts and BG correction on Circular CGM
python3 RM_smoothing.py --b-limit-small --bg-corr --save-plot --_485mhz --HI --smoothed #--scatter
mv Results/"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png" Results/b-limit-small_smoothed.png

#No Latitude cuts or BG correction on Elliptical CGM
python3 RM_smoothing.py --elliptic-CGM --bg-corr --save-plot --_485mhz --HI --smoothed #--scatter
mv Results/"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png" Results/elliptical-CGM_no_bg_smoothed.png

#No Latitude cuts but BG correction on Elliptical CGM
python3 RM_smoothing.py --elliptic-CGM --bg-corr --save-plot --_485mhz --HI --smoothed #--scatter
mv Results/"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png" Results/elliptical-CGM_smoothed.png

#Latitude cuts and BG correction on Elliptical CGM
python3 RM_smoothing.py --b-limit-small --elliptic-CGM --bg-corr --save-plot --_485mhz --HI --smoothed #--scatter
mv Results/"2D_Gaussian_RM_Smoothing_kernel=(2.0, 2.0)_in_degrees.png" Results/b-limit-small-elliptical-CGM_smoothed.png

python3 debugging_cumulative_density.py --smoothed
