echo "for 10:" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_vs_entire_sky.py --pickle --b-limit-small --patch-num 10 --m31-ks-test --rm-per-patch-hist --save-plot | grep "KS Test" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_density.py --b-limit-small --patch-num 10 --save-plot
python3 M31_signal_density.py --b-limit-small --patch-num 10 --cumulative --save-plot
mkdir ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10
echo "MORE: " > ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10/kstest-log_10.txt
python3 M31_signal_vs_entire_sky.py --b-limit-small --patch-num 10 --m31-ks-test --save-plot>> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10/kstest-log_10.txt
mv ./Results/*.png ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10
mv ../*.pkl ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10
echo "_________________________________________________________________________________"

echo "for 100:" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_vs_entire_sky.py --pickle --b-limit-small --patch-num 100 --m31-ks-test --rm-per-patch-hist --save-plot | grep "KS Test" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_density.py --b-limit-small --patch-num 100 --save-plot
python3 M31_signal_density.py --b-limit-small --patch-num 100 --cumulative --save-plot
mkdir ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_100
echo "MORE: " > ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_100/kstest-log_100.txt
python3 M31_signal_vs_entire_sky.py --b-limit-small --patch-num 100 --m31-ks-test --save-plot >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_100/kstest-log_100.txt
mv ./Results/*.png ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_100
mv ../*.pkl ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_100
echo "_________________________________________________________________________________"


echo "for 1000:" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_vs_entire_sky.py --pickle --b-limit-small --patch-num 1000 --m31-ks-test --rm-per-patch-hist --save-plot | grep "KS Test" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_density.py --b-limit-small --patch-num 1000 --save-plot
python3 M31_signal_density.py --b-limit-small --patch-num 1000 --cumulative --save-plot
mkdir ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_1000
echo "MORE: " > ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_1000/kstest-log_1000.txt
python3 M31_signal_vs_entire_sky.py --b-limit-small --patch-num 1000 --m31-ks-test --save-plot >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_1000/kstest-log_1000.txt
mv ./Results/*.png ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_1000
mv ../*.pkl ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_1000
echo "_________________________________________________________________________________"


echo "for 10K:" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_vs_entire_sky.py --pickle --b-limit-small --patch-num 10000 --m31-ks-test --rm-per-patch-hist --save-plot | grep "KS Test" >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/kstest-log.txt
python3 M31_signal_density.py --b-limit-small --patch-num 10000 --save-plot
python3 M31_signal_density.py --b-limit-small --patch-num 10000 --cumulative --save-plot
mkdir ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10K
echo "MORE: " > ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10K/kstest-log_10K.txt
python3 M31_signal_vs_entire_sky.py --b-limit-small --patch-num 10000 --m31-ks-test --save-plot >> ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10K/kstest-log_10K.txt
mv ./Results/*.png ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10K
mv ../*.pkl ../1e1_1e2_1e3_and_1e4_with_b_limit_small_imposed/Results_10K
echo "_________________________________________________________________________________"


