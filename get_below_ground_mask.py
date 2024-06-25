# get geopotential height from model levels at bottom level to give us orography
# then read geopotential height at requested time on pressure levels
import os.path as path
from get_era5_on_plevs import *
import iris
import iris.cube
import iris.analysis
from set_latlon_bounds import *
import datetime as dt
import numpy as np

def create_era5_orog():
    fname='/home/users/jcrook001/FORSEA/orog.nc'
    era5_dir='/badc/ecmwf-era5/data/oper/an_ml/2018/01/01/'
    filename='ecmwf-era5_oper_an_ml_201801010000.z.nc'
    z_cube=iris.load_cube(era5_dir+filename)
    orog=z_cube[0]/9.80665
    iris.save(orog, fname)
    return orog

def get_orography(intersection, for_um=False):
    if for_um:
        orog_file='/gws/nopw/j04/forsea/users/jcrook/um_fcst/20190111T0000Z_SEA4_km4p4_ra1tld_pa000.pp'
        print('reading orog from ', orog_file)
        orog=iris.load_cube(orog_file, 'surface_altitude')
    else:
        fname='/home/users/jcrook001/FORSEA/orog.nc'
        print('reading orography')
        if path.isfile(fname)==True:
            big_intersection= {'latitude': [-20,25], 'longitude': [70,160]} # use this to cover all cases
            orog=iris.load_cube(fname).intersection(**big_intersection)
            orog=orog.intersection(**intersection)
        else:
            orog=create_era5_orog()
    return orog

def get_era5_below_ground_mask(start_date, end_date, intersection, orog):
    geopot=get_era5_on_plevs(start_date, end_date, intersection, 'geopot', 'geopotential',dhours=6)
    geoh=geopot/9.80665
    pressures=geopot.coord('pressure_level').points
    return get_below_ground_mask(orog, geoh, pressures)
    
def get_below_ground_mask(orog, geoh, pressures):

    below_mask=np.ma.zeros_like(geoh.data)
    below_mask.mask=True
    shape=np.shape(below_mask)
    nt=0
    if len(shape)==4:
        nt=shape[0]
    if nt>0:
        for t in range(nt):
            this_below_mask=get_singlet_below_ground_mask(orog, geoh[t,:,:,:], pressures)
            below_mask[t,:,:,:]=this_below_mask
    else:
        below_mask=get_singlet_below_ground_mask(orog, geoh, pressures)

    return below_mask


def get_singlet_below_ground_mask(orog, geoh, pressures):
    nplev=len(pressures)
    below_mask=np.ma.zeros_like(geoh.data)+np.nan
    below_mask.mask=True
    for p in range(nplev):
        ix=np.where((geoh[p,:,:].data-orog.data)<0)
        if len(ix[0])>0:
            this_mask=np.ma.zeros_like(below_mask[p,:,:])+np.nan
            this_mask.mask=True
            this_mask[ix]=1
            this_mask[ix].mask=False
            below_mask[p,:,:]=this_mask
            #print(len(ix[0]), 'points below', pressures[p])
        #else:
            #print('all points above at', pressures[p])
    return below_mask

