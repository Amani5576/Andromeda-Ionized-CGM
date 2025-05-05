#No latitude cuts and no BG correction on circular CGM
python3 M31_signal_vs_entire_sky.py --save-plot --median --CGM-dispersion --patch-num 10000
mv Results/Dispersion_Median.png Results/nothing_no_bg_disp.png

#BG correction without latitude cuts on circular CGM
python3 M31_signal_vs_entire_sky.py --bg-corr --save-plot --median --CGM-dispersion --patch-num 10000
mv Results/Dispersion_Median.png Results/nothing_disp.png

#BG correction and Latitude cuts on Circualr CGM
python3 M31_signal_vs_entire_sky.py --b-limit-small --bg-corr --save-plot --median --CGM-dispersion --patch-num 10000
mv Results/Dispersion_Median.png Results/b-limit-small_disp.png

#No latitude cuts and no BG correction on elliptical CGM
python3 M31_signal_vs_entire_sky.py --elliptic-CGM --save-plot --median --CGM-dispersion --patch-num 10000
mv Results/Dispersion_Median.png Results/elliptical-CGM_no_bg_disp.png

#BG correction without latitude cuts on elliptical CGM
python3 M31_signal_vs_entire_sky.py --elliptic-CGM --bg-corr --save-plot --median --CGM-dispersion --patch-num 10000
mv Results/Dispersion_Median.png Results/elliptical-CGM_disp.png

#BG correction with latitude cuts on elliptical CGM
python3 M31_signal_vs_entire_sky.py --b-limit-small --elliptic-CGM --bg-corr --save-plot --median --CGM-dispersion --patch-num 10000
mv Results/Dispersion_Median.png Results/b-limit-small-elliptical-CGM_disp.png

python3 debugging_cumulative_density.py --CGM-dispersion

