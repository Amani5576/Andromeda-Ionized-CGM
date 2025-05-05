#No latitude cuts and no BG correction on circular CGM
python3 RM_smoothing.py --save-plot --_485mhz --HI --scatter
mv Results/Temperature_Contours_within_CGM_of_M31.png Results/nothing_no_bg_not_smootheD.png

#No latitude cuts but with BG correction on circular CGM
python3 RM_smoothing.py --bg-corr --save-plot --_485mhz --HI --scatter
mv Results/Temperature_Contours_within_CGM_of_M31.png Results/nothing_not_smootheD.png

#Latitude cuts and BG correction on circular CGM
python3 RM_smoothing.py --b-limit-small --bg-corr --save-plot --_485mhz --HI --scatter
mv Results/Temperature_Contours_within_CGM_of_M31.png Results/b-limit-small_not_smootheD.png

#No latitude cuts and no BG correction on elliptical CGM
python3 RM_smoothing.py --elliptic-CGM --save-plot --_485mhz --HI --scatter
mv Results/Temperature_Contours_within_CGM_of_M31.png Results/elliptical-CGM_no_bg_not_smootheD.png

#BG correction without latitude cuts on elliptical CGM
python3 RM_smoothing.py --elliptic-CGM --bg-corr --save-plot --_485mhz --HI --scatter
mv Results/Temperature_Contours_within_CGM_of_M31.png Results/elliptical-CGM_not_smootheD.png

#Latitude cuts and BG correction on ellitpical CGM
python3 RM_smoothing.py --b-limit-small --elliptic-CGM --bg-corr --save-plot --_485mhz --HI --scatter
mv Results/Temperature_Contours_within_CGM_of_M31.png Results/b-limit-small-elliptical-CGM_not_smootheD.png

python3 debugging_cumulative_density.py --temp