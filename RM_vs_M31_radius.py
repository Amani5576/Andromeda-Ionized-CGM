#Importing package alias'
from main import np, stats
import matplotlib.pyplot as plt

#Importing vairables
from main import m31_sep_Rvir, rm_m31, bin_num, d_bg, rm_bg, d_rm, R_vir, d_bin_centers, bin_means

#Using standard error of mean for error bars: 
bin_std, bin_edges, binnumber = stats.binned_statistic(m31_sep_Rvir, 
    rm_m31, 
    statistic = lambda rm_m31: stats.sem(rm_m31), #Standard Error ofMean
    bins = bin_num)  

plt.figure(figsize=(8,8))
plt.plot(d_rm, rm_m31, 'x', markersize=3)
plt.plot(d_bg, rm_bg, 'x', markersize=3)
plt.plot([R_vir.value, R_vir.value], [1000., -1000.], 'r--', linewidth=2)
plt.xlim(0., np.max(d_bg.value))
plt.ylim(-R_vir.value, R_vir.value) #within Virial radius
plt.errorbar(d_bin_centers, bin_means, yerr = bin_std, fmt = 'k.-')
plt.xlabel(r'$r_{M31}$ [kpc]', fontsize=10)
plt.ylabel(r'RM [rad/$m^2$]', fontsize=10)
plt.grid(True, which="both")
plt.show()