# calculate surge data from 925 hPa winds - roughly model level 120 - and store in .nc file
# Constraint
import os.path as path
from get_era5_on_plevs import *
import iris
import iris.cube
import iris.analysis
from set_latlon_bounds import *
import datetime as dt
import numpy as np

def get_cross_surge(v925_cube):

    cross_intersection = {'latitude': [-5,0], 'longitude': [105,115]}
    v925_cross=v925_cube.intersection(**cross_intersection)
    set_latlon_bounds(v925_cross)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        weights_cross = iris.analysis.cartography.area_weights(v925_cross)
    mean_v_cross = v925_cross.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,
                          weights=weights_cross)

    # Find when cross-equatorial northerly surge occurs (v < -5)
    mean_v_cross.convert_units('m s-1')
    #print(mean_v_cross)
    cross_equatorial_northerly_surge = mean_v_cross.data < -5
    cross_equatorial_northerly_surge_coord = iris.coords.AuxCoord(
                                                 cross_equatorial_northerly_surge.astype(int),
                                                 var_name='cross_equatorial_northerly_surge')
    mean_v_cross.add_aux_coord(cross_equatorial_northerly_surge_coord, 0)  
    return mean_v_cross

def get_meridional_surge(v925_cube):
    merid_intersection={'latitude':[14.99,15.01], 'longitude':[110,117.5]}
    v925_merid=v925_cube.intersection(**merid_intersection)
    mean_v_merid = v925_merid.collapsed(['latitude','longitude'], iris.analysis.MEAN) # don't need weights as only 1 latitude

    # Find when meridional surge occurs (v < -8)
    mean_v_merid.convert_units('m s-1')
    #print(mean_v_merid)
    merid_surge = mean_v_merid.data < -8
    merid_surge_coord = iris.coords.AuxCoord(merid_surge.astype(int),
                                             var_name='meridional_surge')
    mean_v_merid.add_aux_coord(merid_surge_coord, 0)
    # Find when weak (-10 <= v < -8), moderate (-12 <= v < -10) and
    # strong (v < -12) meridional surges occur
    for strength, low, high in [('weak', -10, -8), ('moderate', -12, -10),
                                ('strong', -1e20, -12)]:
        meridional_surge = (low <= mean_v_merid.data) & (mean_v_merid.data < high)
        meridional_surge_coord = iris.coords.AuxCoord(
            meridional_surge.astype(int), var_name=strength+'_meridional_surge')
        mean_v_merid.add_aux_coord(meridional_surge_coord, 0)

    return mean_v_merid

def get_easterly_surge(u925_cube):

    set_latlon_bounds(u925_cube)
    mean_u = u925_cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN) # don't need weights as only 1 

    # Find when easterly surge occurs: 
    mean_u.convert_units('m s-1')
    #print(mean_u)
    easterly_surge = mean_u.data < -8
    easterly_surge_coord = iris.coords.AuxCoord(easterly_surge.astype(int),
                                                var_name='easterly_surge')
    mean_u.add_aux_coord(easterly_surge_coord, 0)
    # find when weak (-10 <= u < -8), moderate (-12 <= u < -10) and
    # strong (u < -12)
    for strength, low, high in [('weak', -10, -8), ('moderate', -12, -10),
                                ('strong', -1e20, -12)]:
        easterly_surge = (low <= mean_u.data) & (mean_u.data < high)
        easterly_surge_coord = iris.coords.AuxCoord(easterly_surge.astype(int),
                                                    var_name=strength+'_easterly_surge')
        mean_u.add_aux_coord(easterly_surge_coord, 0)

    return mean_u

def get_surge(year,  months=[1,2,3,10,11,12]):
    v_intersection={'latitude':[-5,16], 'longitude':[100,120]}
    easterly_intersection={'latitude':[7.5,15.0], 'longitude':[119.99,120.01]}
    outdir='/home/users/jcrook001/FORSEA/surge_data/'

    for m in months:
        this_start_date=dt.datetime(year,m,1)
        if m==12:
            this_end_date=dt.datetime(year+1,1,1)-dt.timedelta(days=1)
        else:
            this_end_date=dt.datetime(year,m+1,1)-dt.timedelta(days=1)
        cross_surge_filename=this_start_date.strftime(outdir+'cross_surge_%Y%m.nc')
        merid_surge_filename=this_start_date.strftime(outdir+'merid_surge_%Y%m.nc')
        if path.isfile(cross_surge_filename)==True and path.isfile(merid_surge_filename)==True:
            print('cross and merid surge already done', this_start_date.strftime('%Y%m'))
        else:
            v925month_cube=get_daily_mean_era5(this_start_date, this_end_date, 925, v_intersection, 'v', 'northward_wind')

            print('getting cross surge')
            mean_v_cross=get_cross_surge(v925month_cube)
            print(mean_v_cross)
            iris.save(mean_v_cross, cross_surge_filename)
            print('getting merid surge')
            mean_v_merid=get_meridional_surge(v925month_cube)
            print(mean_v_merid)
            iris.save(mean_v_merid, merid_surge_filename)
 
        east_surge_filename=this_start_date.strftime(outdir+'easterly_surge_%Y%m.nc')
        if path.isfile(east_surge_filename)==True:
            print('easterly surge already done', this_start_date.strftime('%Y%m'))
        else:
            u925month_cube=get_daily_mean_era5(this_start_date, this_end_date, 925, easterly_intersection, 'u', 'eastward_wind')
            print('getting easterly surge')
            mean_u=get_easterly_surge(u925month_cube)
            iris.save(mean_u, east_surge_filename)


