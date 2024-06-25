# get ERA5 data (u,v,omega, T or q) on pressure levels as stored in Terramaris GWS
# Note that this does not have the 925 hPa pressure level
import os.path as path
import iris
import iris.cube
import datetime as dt
import numpy as np

def get_data_dir(this_date, var_str):

    era5_dir='/gws/nopw/j04/terramaris/era5/'+var_str+'/'
    era5u201810_dir='/gws/nopw/j04/terramaris/era5_u/u/'
    start_special=dt.datetime(2018,10,5)
    end_special=dt.datetime(2018,10,31,23)
    this_dir=era5_dir
    if var_str=='u' and this_date>=start_special and this_date<=end_special:
        this_dir=era5u201810_dir
    data_dir=this_dir+this_date.strftime('%Y/%m/%d/')

    return data_dir

def in_subset_pressure(cell,subset_pressures):
    in_subset=np.asarray([cell==pressure for pressure in subset_pressures])
    return len(np.where(in_subset)[0])>0
    
def get_era5_on_plevs(start_date, end_date, intersection, var_str, var_name, dhours=6, pressure=None, subset_pressures=None):

    # data is in 6 hourly files
    str_p=''
    level_constraint=None
    if pressure!=None:
        level_constraint=iris.Constraint(pressure_level=pressure)
        str_p=' on pressure level {p:d} '.format(p=pressure)
    elif subset_pressures!=None:
        level_constraint=iris.Constraint(pressure_level=lambda cell: in_subset_pressure(cell,subset_pressures))
        str_p=' on pressure levels '
        for p in subset_pressures:
            str_p=str_p+'{p:d} '.format(p=p)
        str_p=str_p+' hPa'
        
    these_cubes=iris.cube.CubeList()
    this_date=start_date
    print('reading',var_str, start_date, 'to', end_date, str_p)
    while(this_date<=end_date):
        data_dir=get_data_dir(this_date, var_str)
        # Read in data
        filename='era5.'+var_str+this_date.strftime('.%Y%m%d_%H00.nc')
        print('reading',data_dir+filename)
        this_cube=iris.load_cube(data_dir+filename,var_name).intersection(**intersection)
        if level_constraint!=None:
            this_cube=this_cube.extract(level_constraint)
        these_cubes.append(this_cube)
        this_date=this_date+dt.timedelta(hours=dhours)
    if len(these_cubes)>1:
        # concatenate times
        iris.util.unify_time_units(these_cubes)
        iris.util.equalise_attributes(these_cubes)
        these_cubes=these_cubes.merge_cube()
    elif len(these_cubes)==1:
        these_cubes=these_cubes[0]
    return these_cubes

# note we only have this for Jan-March and Oct-Dec each year
def get_era5_925(start_date, end_date, intersection, var_str, var_name, dhours=1):

    era5_dir='/gws/nopw/j04/terramaris/era5_925/'+var_str+'/'
    # data is in hourly files
    these_cubes=iris.cube.CubeList()
    this_date=start_date
    print('reading 925',var_str, start_date, 'to', end_date)
    while(this_date<=end_date):
        data_dir=era5_dir+this_date.strftime('%Y/%m/%d/')
        # Read in data
        filename='era5.'+var_str+this_date.strftime('.%Y%m%d_%H00.nc')
        this_cube=iris.load_cube(data_dir+filename,var_name).intersection(**intersection)
        these_cubes.append(this_cube)
        this_date=this_date+dt.timedelta(hours=dhours)
    if len(these_cubes)>1:
        # concatenate times
        iris.util.unify_time_units(these_cubes)
        iris.util.equalise_attributes(these_cubes)
        these_cubes=these_cubes.merge_cube()
    elif len(these_cubes)==1:
        these_cubes=these_cubes[0]

    return these_cubes

def get_daily_mean_era5(start_date, end_date, pressure, intersection, var_str, var_name,dhours=1):

    # data is in hourly files
    these_cubes=iris.cube.CubeList()
    this_date=start_date
    nhours=int(24/dhours)
    print('reading',var_str, start_date, 'to', end_date)
    level_constraint=iris.Constraint(pressure_level=pressure)
    while(this_date<=end_date):
        if pressure==925:
            data_dir='/gws/nopw/j04/terramaris/era5_925/'+var_str+'/'+this_date.strftime('%Y/%m/%d/')
        else:
            data_dir=get_data_dir(this_date, var_str)
        # Read in data
        this_hour=this_date
        era5_day_cubes=iris.cube.CubeList()
        for h in range(nhours):
            filename='era5.'+var_str+this_hour.strftime('.%Y%m%d_%H00.nc')
            print('reading',filename)
            this_cube=iris.load_cube(data_dir+filename,var_name).intersection(**intersection)
            era5_day_cubes.append(this_cube.extract(level_constraint))
            this_hour=this_hour+dt.timedelta(hours=dhours)
        # concatenate times for this day
        iris.util.unify_time_units(era5_day_cubes)
        iris.util.equalise_attributes(era5_day_cubes)
        era5_day_cubes=era5_day_cubes.merge_cube()
        # get daily mean
        era5_day_cubes=era5_day_cubes.collapsed('time', iris.analysis.MEAN)
        these_cubes.append(era5_day_cubes)
        this_date=this_date+dt.timedelta(days=1)
    # concatenate times
    iris.util.unify_time_units(these_cubes)
    iris.util.equalise_attributes(these_cubes)
    these_cubes=these_cubes.merge_cube()

    return these_cubes

