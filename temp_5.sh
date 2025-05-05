#No latitude cuts and no BG correction on circular CGM
python3 M31_signal_density.py --save-plot --median --cumulative --patch-num 10000
mv Results/Cumulative_Density_bw=0.15_10000_patches_median.png Results/nothing_no_bg.png

#BG correction without latitude cuts on circular CGM
python3 M31_signal_density.py --bg-corr --save-plot --median --cumulative --patch-num 10000
mv Results/Cumulative_Density_bw=0.15_10000_patches_median.png Results/nothing.png

#Latitude cuts with BG correcttion on circular CGM
python3 M31_signal_density.py --b-limit-small --bg-corr --save-plot --median --cumulative --patch-num 10000
mv Results/Cumulative_Density_bw=0.15_10000_patches_median.png Results/b-limit-small.png

#No latitude cuts and no BG correction on elliptical CGM
python3 M31_signal_density.py --elliptic-CGM --save-plot --median --cumulative --patch-num 10000
mv Results/Cumulative_Density_bw=0.15_10000_patches_median.png Results/elliptical-CGM_no_bg.png

#No Latitude cuts with backgorund corrrection
python3 M31_signal_density.py --elliptic-CGM --bg-corr --save-plot --median --cumulative --patch-num 10000
mv Results/Cumulative_Density_bw=0.15_10000_patches_median.png Results/elliptical-CGM.png

#Latitude cuts with background corrrection
python3 M31_signal_density.py --b-limit-small --elliptic-CGM --bg-corr --save-plot --median --cumulative --patch-num 10000
mv Results/Cumulative_Density_bw=0.15_10000_patches_median.png Results/b-limit-small-elliptical-CGM.png

python3 debugging_cumulative_density.py
