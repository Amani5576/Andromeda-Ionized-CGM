 #Imorting alias's from attempt
from main import np, Table, u, SkyCoord, plt

def get_real_rm_data(): #get rm_data from catalog.dat
    #Read in catalog.dat file
    t = Table.read("catalog.dat", format = "ascii")

    t.rename_column("col9", "l")
    t.rename_column("col10", "b")
    t.rename_column("col17", "RM")
    t.rename_column("col18", "e_RM")

    t.keep_columns(["l", "b", "RM", "e_RM"])

    RM_unit = u.rad/u.m**2

    t["l"].unit = u.deg
    t["b"].unit = u.deg
    t["RM"].unit = RM_unit
    t["e_RM"].unit = RM_unit

    position = SkyCoord(t["l"], t["b"], unit = "degree", frame = "galactic")
    eq_pos = position.transform_to('icrs') #Equitorial positions
    
    return eq_pos, t["RM"], t["e_RM"]

def get_points(num=2000, ra_range=(0, 360), dec_range=(-90, 90), random=True,
               accounted_for_poles=True):

    def convert_to_cartesian(ra, dec):
        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)
        return x, y, z
    
    if random:
        ra = np.random.uniform(*ra_range, num)  # This is in degrees
        ra_rad = np.deg2rad(ra)  # Conversion to radians
        
        if accounted_for_poles:#Proper declination distribution on sphere
            dec_rad = np.arcsin(np.random.uniform(np.sin(np.deg2rad(dec_range[0])), 
                                                  np.sin(np.deg2rad(dec_range[1])), 
                                                  num))  # Gives back in radians
        else: #Improper declination distribution on sphere that collect more along the poles
            dec_rad = np.random.uniform(np.deg2rad(dec_range[0]), 
                        np.deg2rad(dec_range[1]), num)  # Gives back in radians
        
        return convert_to_cartesian(ra_rad, dec_rad)

    else: #Deal with real rm data
        eq_pos, _, _ = get_real_rm_data()
        return convert_to_cartesian(eq_pos.ra.rad, eq_pos.dec.rad)
        
# Function to calculate dot product with viewing direction
def calculate_dot_product(x, y, z, theta, phi):
    theta, phi = np.deg2rad(theta), np.deg2rad(phi)
    
    # Calculate the direction vector from theta (elevation) and phi (azimuth)
    x_dir = np.cos(phi) * np.cos(theta)
    y_dir = np.sin(phi) * np.cos(theta)
    z_dir = np.sin(theta)
    
    # Dot product between the direction vector and point (x, y, z)
    return x * x_dir + y * y_dir + z * z_dir

# Function to plot the sphere's surface
def plot_sphere_surface(ax):
    # Generate spherical coordinates
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]

    # Convert spherical coordinates to Cartesian coordinates
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    # Plot the surface of the sphere
    ax.plot_surface(x, y, z, color='k', alpha=0.1)  # Set alpha for transparency

# Function to plot the points on the sphere
def plot_points(ax, x, y, z, theta, phi, **kw):
    # Compute dot product with the specified direction
    dot_product = calculate_dot_product(x, y, z, theta, phi)

    # Classify points as visible (blue) or not visible (red) based on dot product
    visible = dot_product >= 0
    not_visible = dot_product < 0

    # Plot the points
    s_1, s_2, c_1, c_2 = (.1,.1,"b","r")
    
    if kw["labelling"]:

        #stuff for sphere's legend    
        ax.scatter([],[],marker='o', s=1, color = 'red')
        ax.scatter([],[],marker='o', s=1, color = 'blue')
        
        handles, labels = ax.get_legend_handles_labels()  # Get current handles and labels
        for handle in handles:
            handle.set_sizes([100])#Setting marker size for legend for clarity (above WCS acxis)
        kw["fig"].legend(fontsize = 25, loc = 'upper center', bbox_to_anchor = (.5, 1.1), 
                    framealpha = 0, ncols = 2)
        
    ax.scatter(x[visible], y[visible], z[visible], c=c_1, marker='o', s=s_1)
    ax.scatter(x[not_visible], y[not_visible], z[not_visible], c=c_2, 
               marker='o', s=s_2)

def get_circle(theta_center, phi_center, radius, num_points=1000):
    
    angles = np.linspace(0, 2*np.pi, num_points) #Generating angles to describe circle

    # Spherical coordinates of the circle (small circle on the surface of the sphere)
    theta_circle = theta_center + radius * np.cos(angles)  # Vary latitude
    phi_circle = phi_center + radius * np.sin(angles)  # Vary longitude

    # Convert the circle points to Cartesian coordinates
    x_circle = np.cos(theta_circle) * np.cos(phi_circle)
    y_circle = np.cos(theta_circle) * np.sin(phi_circle)
    z_circle = np.sin(theta_circle)

    return x_circle, y_circle, z_circle

def store_into_txt(filename, RA, DEC):

    with open(filename, "w") as file:
        file.write("RA (degrees)    DEC (degrees)\n")
        for ra, dec in zip(RA, DEC):
            file.write(f"{ra:<15} {dec}\n")
    
    print(f"Data sotred in {filename}")

# Set up the figure and 3D axis
fig, axes = plt.subplots(3,3, figsize=(15,15),
                         subplot_kw={'projection': '3d'})

axes = axes.flatten()

#IMPORTANT
np.random.seed(0)
#default has correct declination distribution otherwise add 
#"accounted_for_poles=False" to see collection of points near the poles
                #Dec is already in radians
x, y, z = get_points(random=True, num=200_000) 
    
# elevation angle for random points on the sky
elev_angles = [ 66, -45,   2,  77,  32,  29,  11, -42,  86] #in degrees
azim_angles = [161, 112, 262, 219, 123, 109,  77,  28, 225] #in degrees

labelling=True
for ax, elev, azim in list(zip(axes, elev_angles, azim_angles)):

    # Set up axis limits and make it symmetrical
    ax.set_aspect('equal')
    
    # Hide axis for better visualization
    ax.set_axis_off()
    
    # Create the surface of the sphere
    plot_sphere_surface(ax)
    
    plot_points(ax, x, y, z, elev, azim, labelling=labelling, fig=fig)
    labelling=False #Otherwise ther will be repetition in the legend
    
    lim = 0.8 #Zooming in parametre
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    
    # Adjust the viewing angle
    ax.view_init(elev, azim)  # Keep elevation fixed, rotate azimuthal angle
    ax.set_title("$(\\alpha,\\delta)=$" + f"({str(azim)}"+ "$^{\circ}$" + f",{str(elev)}"+ "$^{\circ})$", fontsize=20)

plt.tight_layout()
plt.show()

n_bins = 20
ra_bins = np.linspace(0, 360, n_bins + 1)  # RA bins from 0 to 360 degrees
dec_bins = np.linspace(np.cos(np.deg2rad(90)), 
                       np.cos(np.deg2rad(0)),
                       n_bins + 1)  # DEC bins from -90 to 90 degrees

# Convert back to spherical coordinates
r = np.sqrt(x**2 + y**2 + z**2)

dec_rad = np.arcsin(z/r)
DEC = np.rad2deg(dec_rad)
ra_rad = np.arctan2(y, x)

sin_Dec = np.sin(dec_rad)
RA = np.rad2deg(ra_rad)

# store_into_txt("Random_ra_dec_points.txt", RA, DEC)

# Sum counts per RA bin
ra_hist, _ = np.histogram(RA, bins=ra_bins)

# Sum counts per DEC bin
dec_hist, _ = np.histogram(sin_Dec, bins=dec_bins)

# Plot RA histogram
plt.figure(figsize=(10, 5))
plt.hist(np.rad2deg(ra_rad), bins=n_bins)
plt.xlabel('RA (degrees)')
plt.ylabel('Counts')
plt.grid(True)
plt.show()

# Plot histogram of cos(dec)
plt.figure(figsize=(10, 5))
plt.hist(np.rad2deg(dec_rad), bins=n_bins)
plt.xlabel('DEC (degrees)')
plt.ylabel('Counts')
plt.grid(True)
plt.show()

# Plot histogram of cos(dec)
plt.figure(figsize=(10, 5))
plt.hist(sin_Dec, bins=n_bins)
plt.xlabel('sin(DEC)')
plt.ylabel('Counts')
plt.grid(True)
plt.show()

