# Concept of Conversion: https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_coords.php
# Code used to make script: https://github.com/adrn/gala/compare/2dc48335da3d92ef3f494d4b6d7c0b5451557dce..9d5f9f6b9cbdf8b9c129dbf4efc447239e5d55e3

from astropy.coordinates.matrix_utilities import (rotation_matrix,
                                matrix_product, matrix_transpose)
from astropy.coordinates.baseframe import (frame_transform_graph, 
                                           BaseCoordinateFrame, 
                                           RepresentationMapping)
from astropy.coordinates.transformations import StaticMatrixTransform
from astropy.coordinates import representation as r
from astropy.coordinates import Galactic
from astropy.coordinates import ICRS

import numpy as np
from astropy.coordinates import SkyCoord

import astropy.units as u


class MagellanicStream(BaseCoordinateFrame):
    """
    A coordinate or frame aligned with the Magellanic Stream, 
    as defined by Nidever et al. (2008).
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'L'),
            RepresentationMapping('lat', 'B')
        ]
    }

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    _ngp = Galactic(l=188.5*u.deg, b=-7.5*u.deg)
    _lon0 = Galactic(l=280.47*u.deg, b=-32.75*u.deg)

@frame_transform_graph.transform(StaticMatrixTransform, 
                                 Galactic, MagellanicStream)
def gal_to_mag():
    mat1 = rotation_matrix(57.275785782128686*u.deg, 'z')
    mat2 = rotation_matrix(90*u.deg - MagellanicStream._ngp.b, 'y')
    mat3 = rotation_matrix(MagellanicStream._ngp.l, 'z')

    return matrix_product(mat1, mat2, mat3)

@frame_transform_graph.transform(StaticMatrixTransform, 
                                 MagellanicStream, Galactic)
def mag_to_gal():
    return matrix_transpose(gal_to_mag())


def ms_to_radec(Lambda_MS, Beta_MS, l_pole=188.5, epsilon=7.5):
    Lambda_MS = np.radians(Lambda_MS)
    Beta_MS = np.radians(Beta_MS)
    l_pole = np.radians(l_pole)
    epsilon = np.radians(epsilon)

    sin_b = np.sin(Beta_MS) * np.cos(epsilon) - np.cos(Beta_MS) * np.sin(epsilon) * np.sin(Lambda_MS)
    b = np.arcsin(sin_b)

    y = np.sin(Beta_MS) * np.sin(epsilon) + np.cos(Beta_MS) * np.cos(epsilon) * np.sin(Lambda_MS)
    x = np.cos(Beta_MS) * np.cos(Lambda_MS)
    l = np.arctan2(y, x) + l_pole

    l = np.mod(l, 2 * np.pi)

    l = np.degrees(l)
    b = np.degrees(b)

    galactic_coord = SkyCoord(l=l * u.deg, b=b * u.deg, frame='galactic')
    ra_dec_coord = galactic_coord.icrs
    return ra_dec_coord.ra.deg, ra_dec_coord.dec.deg

def convert_magellanic_stream_coords(l_ms, b_ms, frame="ICRS"):
    
    typ = ICRS if frame.lower() == "icrs" else Galactic
    
    return SkyCoord(L=l_ms, 
             B=b_ms, 
             unit=(u.deg, u.deg), 
             frame=MagellanicStream).transform_to(typ)


def convert_string_where_applicable(table):
    
    """This makes sure that the RA and DEC found in the McConnachie paper
        are utitlised instead of relying on Lehner et. al. Table 3. 
        However, Not all dwarf galaxies from Lehner et. al Table 3 are
        Listed in the McConnachie. Those missing stallites are then acquired 
        via conversion from lehner et. al. Table 3 Magellanic Stream to
        ICRS.
    """
    #Making SkyCoord object
    coord = SkyCoord(ra=table["RA"], dec=table["DEC"], unit=(u.hourangle, u.deg), frame='icrs')
    
    # No need to return since they are astropy tables are mutable
    table["RA"], table["DEC"] = coord.ra.deg, coord.dec.deg
    
    #If There wasnt a value then dont give o, but rather convert from l_ms and b_ms 
    table["RA"] = np.where(table["RA"] == 0, 
                           convert_magellanic_stream_coords(table["l_MS"], 
                                                            table["b_MS"], 
                                                            frame="ICRS").ra.deg, 
                           table["RA"])
    table["DEC"] = np.where(table["DEC"] == 0, 
                           convert_magellanic_stream_coords(table["l_MS"], 
                                                            table["b_MS"], 
                                                            frame="ICRS").dec.deg,
                           table["DEC"])
    
    


#Just storing this here because I really dont wanna see it in other scripts
alpha = lambda: round(np.random.uniform(.3,.6), 3)
def generate_marker_color_pairs(num_pairs):
    
    #Defining marker styles and colors
    markers = ['>', 's', '*', '^', 'o']
    colors = ['b', 'g', 'r', 'm', 'y', 'c']  # Added 'c' for cyan

    mca = [] #Marker Color Alpha
    
    #Generating pairs
    for i in range(num_pairs):
        marker = markers[i % len(markers)]  #Cycling through markers
        color = colors[i % len(colors)]      #Cycling through colors
        mca.append((marker, color, alpha())) #Adding tuple

    return mca


