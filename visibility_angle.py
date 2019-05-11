import numpy as np
import matplotlib.pyplot as plt
import scipy
import xarray
from .data import load_grt
from pathlib import *


#def calculate_vis_ang(geo_data, observer):


    #lats=geo_data.lat
    #longs=geo_data.lat
    #observer[0]
    #observer[1]

Rz=6.378e6
lato=16
lono=50
alto=1200
latp=17
lonp=51
altp=450

lonp=np.linspace(0,90)
latp=np.linspace(0,90)
alts=np.sin(lonp)
coords=(latp,lonp) #latp,lonp

coords_rad=np.radians(coords) #convert to radians

#print(coords_rad)



def geo_to_cart(lat,lon,alt):
    x=(Rz+alt)*np.sin(lon)*np.cos(lat)
    y=(Rz+alt)*np.cos(lon)*np.cos(lat)
    z=(Rz+alt)*np.sin(lat)
    #print(Rz+alt)
    return np.vstack((x,y,z))


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))



Z=geo_to_cart(coords_rad[0],coords_rad[1],alts)

O=geo_to_cart(np.radians(lono),np.radians(lato),alto)

OZ=Z-O
print(np.shape(O))
print(np.shape(OZ))
A=np.apply_along_axis(angle_between,1,arr=OZ.T,v2=O.T[0])
print(np.shape(A))
plt.plot(lonp,A)
plt.show()

project_path=Path.cwd().parent[1]
file_path=Path.joinpath(project_path, 'Data', 'GMRTv3.6','CR_SR.grd')

d=load_grt()
#print(A)