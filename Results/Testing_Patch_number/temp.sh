python3 M31_signal_density.py --pickle --cumulative --save-plot
python3 M31_signal_density.py --save-plot
echo "for 10K:" > kstest-log.txt
python3 M31_signal_vs_entire_sky.py --m31-ks-test | grep P-value >> kstest-log.txt
python3 M31_signal_density.py --save-plot
rm -rf ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10000_patches)'
mkdir ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10000_patches)'
mv ../*.pkl ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10000_patches)'
echo "_______________________________________________________________"


python3 M31_signal_density.py --pickle --cumulative --patch-num 1000 --save-plot
echo "for 1000:" >>kstest-log.txt
python3 M31_signal_vs_entire_sky.py --m31-ks-test | grep P-value >> kstest-log.txt
python3 M31_signal_density.py --patch-num 1000 --save-plot
rm -rf ../'removed_|b|_<_5_(10_radial_bins,_1st_order,1000_patches)'
mkdir ../'removed_|b|_<_5_(10_radial_bins,_1st_order,1000_patches)'
mv ../*.pkl ../'removed_|b|_<_5_(10_radial_bins,_1st_order,1000_patches)'
echo "_______________________________________________________________"


python3 M31_signal_density.py --pickle --cumulative --patch-num 100 --save-plot
echo "for 100:" >>kstest-log.txt
python3 M31_signal_vs_entire_sky.py --m31-ks-test | grep P-value >> kstest-log.txt
python3 M31_signal_density.py --patch-num 100 --save-plot
rm -rf ../'removed_|b|_<_5_(10_radial_bins,_1st_order,100_patches)'
mkdir ../'removed_|b|_<_5_(10_radial_bins,_1st_order,100_patches)'
mv ../*.pkl ../'removed_|b|_<_5_(10_radial_bins,_1st_order,100_patches)'
echo "_______________________________________________________________"


python3 M31_signal_density.py --pickle --cumulative --patch-num 10 --save-plot
echo "for 10:" >> kstest-log.txt
python3 M31_signal_vs_entire_sky.py --m31-ks-test | grep P-value >> kstest-log.txt
python3 M31_signal_density.py --patch-num 10 --save-plot
rm -rf ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10_patches)'
mkdir ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10_patches)'
mv ../*.pkl ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10_patches)'
echo "_______________________________________________________________"