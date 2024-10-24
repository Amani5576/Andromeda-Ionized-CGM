from main import (
    
# Importing functions
dens_verses_rm,
characterize_densVersesRM,
assessing_boxScatter_with_HI,

#Importing package alias'
plt,

#Importing variables
rm_m31, m31_pos, m31_pa, rm_pos_icrs, err_m31, m31_condition, eq_pos,

)

#Creating a scatter plot of RM verses the Column Density
plt.figure(figsize=(7, 5))
plt.title("Within Virial Radius", fontsize=18, pad=24)
dens_verses_rm(ra=eq_pos.ra[m31_condition],
                dec=eq_pos.dec[m31_condition], 
                rm_errors= err_m31,
                above_zero=None, 
                rm_values=rm_m31)
characterize_densVersesRM()
plt.show()

# assessing_boxScatter_with_HI("Within Boxed filament Region")

# # If you uncomment this, you must change the information in the function to
# # reference the correct confining region.
# assessing_boxScatter_with_HI("Focusing on M31")