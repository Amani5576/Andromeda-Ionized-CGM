import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

#setting number of random points to generate
n_points = 10000

#sampling right ascension values uniformly between 0 and 360 degrees
ra_vals = np.random.uniform(0, 360, n_points)

#defining function for doing sin(dec) sampling and converting to galactic
def testing(val):
    #lower and upper bounds of dec in degrees
    dec_min_deg, dec_max_deg = -90+val, 90
    dec_min_rad, dec_max_rad = np.deg2rad(dec_min_deg), np.deg2rad(dec_max_deg) #converting to radians

    #sine of declination bounds
    sin_min, sin_max = np.sin(dec_min_rad), np.sin(dec_max_rad)
    sin_vals = np.random.uniform(sin_min, sin_max, n_points)
    
    dec_rad = np.arcsin(sin_vals) #declination in radians
    dec_vals = np.rad2deg(dec_rad) #to degrees

    #computing sin(dec) for plotting
    sin_vals = np.sin(np.deg2rad(dec_vals))

    #creating sky coordinate object ith sampled centre points
    coords = SkyCoord(ra=ra_vals * u.deg, dec=dec_vals * u.deg, frame="icrs")

    b_vals = coords.galactic.b.deg
    sin_b_vals = np.sin(np.deg2rad(b_vals)) #sin(b) for plotting

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    
    axs[0].hist(sin_vals, bins=30, edgecolor="black")
    axs[0].set_title("Uniform Sampling in sin(Dec)")
    axs[0].set_xlabel("sin(Dec)")
    axs[0].set_ylabel("Number of Points")
    
    axs[1].hist(sin_b_vals, bins=30, edgecolor="black")
    axs[1].set_title("Distribution in sin(b) (Galactic Latitude)")
    axs[1].set_xlabel("sin(b)")
    axs[1].set_ylabel("Number of Points")

    fig.suptitle(f"{dec_min_deg} < dec < {dec_max_deg}")
    plt.tight_layout()
    plt.show()

#running testing function for different values of declination lower bounds
for i in [0, 10, 20, 30, 40, 50, 60, 70] + list(range(70,91,2)):
    testing(i)
