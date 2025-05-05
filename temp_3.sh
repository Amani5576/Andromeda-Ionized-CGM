python3 M31_signal_vs_entire_sky.py --pickle 
mv ../*.pkl ../'NON-removed_|b|_<_5_(10_radial_bins,_1st_order,10000_patches)'

python3 M31_signal_vs_entire_sky.py --pickle --patch-num 1000 
mv ../*.pkl ../'NON-removed_|b|_<_5_(10_radial_bins,_1st_order,1000_patches)'

python3 M31_signal_vs_entire_sky.py --pickle  --patch-num 100 
mv ../*.pkl ../'NON-removed_|b|_<_5_(10_radial_bins,_1st_order,100_patches)'

python3 M31_signal_vs_entire_sky.py --pickle  --patch-num 10 
mv ../*.pkl ../'NON-removed_|b|_<_5_(10_radial_bins,_1st_order,10_patches)'