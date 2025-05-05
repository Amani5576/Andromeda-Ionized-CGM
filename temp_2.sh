rm ../*.pkl
cp -r ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10000_patches)'/* ..
python3 M31_signal_vs_entire_sky.py --m31-ks-test --patch-num 10000
echo "_______________________________________________________________"


rm ../*.pkl
cp ../'removed_|b|_<_5_(10_radial_bins,_1st_order,1000_patches)'/* ..
python3 M31_signal_vs_entire_sky.py --m31-ks-test  --patch-num 1000
echo "_______________________________________________________________"

rm ../*.pkl
cp ../'removed_|b|_<_5_(10_radial_bins,_1st_order,100_patches)'/* ..
python3 M31_signal_vs_entire_sky.py --m31-ks-test  --patch-num 100
echo "_______________________________________________________________"

rm ../*.pkl
cp ../'removed_|b|_<_5_(10_radial_bins,_1st_order,10_patches)'/* ..
python3 M31_signal_vs_entire_sky.py --m31-ks-test  --patch-num 10
echo "_______________________________________________________________"
