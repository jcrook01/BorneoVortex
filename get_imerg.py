# get GPM INERG data as stored in BADC CEDA archive
import os.path as path
import iris
import iris.cube
import datetime as dt
import numpy as np
import warnings
from pathlib import Path
import datetime
import numpy as np
import iris
import h5py

# function to read the actual HDF5 30 minute files
def read_gpm(filenames=None, times=None, var_name='precipitationCal',slice_lat=slice(None), slice_lon=slice(None), intersection=None,constraint=None, missing='error'):
    """Reads in a variable from a GPM file and returns an iris Cube.
    Note that this function uses the h5py library to read in GPM files
    and, therefore, does NOT support lazy loading.

    Exactly one of *filenames* (string(s) or pathlib.Path(s)) or *times*
    (datetime object(s)) must be given.

    Available *var_name*s: precipitationCal, precipitationUncal,
    randomError, HQprecipitation, HQprecipSource, HQobservationTime,
    IRprecipitation, IRkalmanFilterWeight,
    probabilityLiquidPrecipitation, precipitationQualityIndex.

    Optional *slice_lat* and *slice_lon* are slice objects.

    Optional *intersection* and *constraint* are applied, in that order,
    to each Cube as follows:
        cube = cube.intersection(**intersection)
        cube = cube.extract(constraint)

        *missing* determines how to handle missing files:
        'error': raise Exception
        'warning': issue warning and move on
        'ignore': move on silently
    """
    # Check arguments
    if [filenames, times].count(None) != 1:
        raise ValueError('Exactly one of filenames and times must be given')
    if times is not None:
        if isinstance(times[0], datetime.datetime):
#           times = [times]
            filenames = []
            for time in times:
                #print(time)
                t1 = time.strftime('%Y%m%d-S%H%M%S')
                time2 = time + datetime.timedelta(minutes=30)
                time2 -= datetime.timedelta(seconds=1)
                t2 = time2.strftime('E%H%M%S')
                mins = time.hour*60 + time.minute
                this_dir=time.strftime('/badc/gpm/data/GPM-IMERG-v6/%Y/%j/')
                filename = this_dir+'3B-HHR.MS.MRG.3IMERG.{t1}-{t2}.{mins:04d}.V06B.HDF5'.format(t1=t1,t2=t2, mins=mins)
                filenames.append(filename)
    else:
        if isinstance(filenames, (str, Path)):
            filenames = [filenames]

    # Make sure slices are valid
    if not all((isinstance(slice_lat, slice), isinstance(slice_lon, slice))):
        raise TypeError('*slice_lat* and *slice_lon* must be slice objects')

    # Iterate for each file
    return_cubes = iris.cube.CubeList()
    time_fmt = '%Y%m%d-S%H%M%S'
    time_units_fmt = 'hours since %Y-%m-%d %H:%M'
    slices = (slice(None), slice_lon, slice_lat)  # HDF files have lon then lat
    for filename in filenames:

        # Get data
        try:
            print('opening ',filename)
            group = h5py.File(filename, 'r')['Grid']
        except OSError as err:
            if missing == 'error':
                print(filename)
                raise err
            elif missing == 'warning':
                print(filename)
                warnings.warn('Missing file {filename}'.format(filename=filename), UserWarning)
                continue
            elif missing == 'ignore':
                continue
        var = group[var_name]
        atts = dict(var.attrs)
        try:
            var_arr = np.ma.MaskedArray(var[slices],mask=var[slices]==atts['_FillValue'])
        except KeyError:
            var_arr = var[slices]
        cube = iris.cube.Cube(var_arr, units=atts.get('units').decode())
        if cube.units == 'mm/hr':
            cube.units = 'mm hr-1'

        # Get dimensions
        dim_names = dict(var.attrs)['DimensionNames'].decode().split(',')
        lat_name = [d for d in dim_names if d.lower().startswith('lat')][0]
        lon_name = [d for d in dim_names if d.lower().startswith('lon')][0]
        for ii, dim_name in enumerate(dim_names):
            file_dim = group[dim_name]
            atts = dict(file_dim.attrs)
            units = atts['units'].decode()
            standard_name = atts['standard_name'].decode()
            coord_slice = {lat_name: slice_lat, lon_name: slice_lon}.get(
            dim_name, slice(None))
            coord = iris.coords.DimCoord(file_dim[coord_slice], units=units,standard_name=standard_name)
            cube.add_dim_coord(coord, ii)

        # Set metadata
        if var_name in ['precipitationCal', 'precipitationUncal','HQprecipitation', 'IRprecipitation']:
            cube.var_name = 'pcp'
            cube.standard_name = 'lwe_precipitation_rate'
            cube.long_name = var_name
        else:
            cube.var_name = var_name

        # Constrain
        if intersection is not None:
            cube = cube.intersection(**intersection)
        if constraint is not None:
            cube = cube.extract(constraint)

        # Append to list
        return_cubes.append(cube)

    # Quick attempt to concatenate
    iris.util.unify_time_units(return_cubes)
    return_cubes = return_cubes.concatenate()

    # Transpose axes
    for c in return_cubes:
        c.transpose(tuple(c.coord_dims(c.coord(axis=ax))[0] for ax in 'TYX'))

    # Return
    if len(return_cubes) == 1:
        return return_cubes[0]
    else:
        return return_cubes

# read the 30 minute files but take mean over 1 hour periods
def get_hourly_imerg(start_date, end_date, intersection):

    # data is in half hourly files so for every 6 hour period read all the files and take mean
    this_imerg_cubes=iris.cube.CubeList()
    this_date=start_date
    while(this_date<=end_date):
        # Read in data and average each hour of day
        dates=[this_date-dt.timedelta(minutes=30),this_date]
        imerg=read_gpm(times=dates, intersection=intersection)
        # take mean of these times
        imerg=imerg.collapsed('time', iris.analysis.MEAN)
        this_imerg_cubes.append(imerg)
        this_date=this_date+dt.timedelta(hours=1)
        
    this_imerg_cubes=this_imerg_cubes.merge_cube()
    return this_imerg_cubes

# read the 30 minute files but take mean over 6 hour periods
def get_6hourly_imerg(start_date, end_date, intersection):

    # data is in half hourly files so for every 6 hour period read all the files and take mean
    this_imerg_cubes=iris.cube.CubeList()
    this_date=start_date
    while(this_date<=end_date):
        # Read in data and average each hour of day
        dates=[this_date-dt.timedelta(minutes=180),this_date-dt.timedelta(minutes=150), this_date-dt.timedelta(minutes=120),this_date-dt.timedelta(minutes=90),this_date-dt.timedelta(minutes=60),this_date-dt.timedelta(minutes=30),this_date,this_date+dt.timedelta(minutes=30),this_date+dt.timedelta(minutes=60),this_date+dt.timedelta(minutes=90), this_date+dt.timedelta(minutes=120), this_date+dt.timedelta(minutes=150)]
        imerg=read_gpm(times=dates, intersection=intersection)
        # take mean of these times
        imerg=imerg.collapsed('time', iris.analysis.MEAN)
        this_imerg_cubes.append(imerg)
        this_date=this_date+dt.timedelta(hours=6)
        
    this_imerg_cubes=this_imerg_cubes.merge_cube()
    return this_imerg_cubes

def get_daily_imerg(start_date, end_date, intersection):

    # data is in half hourly files so for every 24 hour period read all the files and take mean
    this_day=start_date.replace(hour=0,minute=0)
    this_imerg_cubes=iris.cube.CubeList()
    while(this_day<=end_date):
        # Read in data and average each hour of day
        dates=[]
        this_date=this_day
        for h in range(24):
           dates.append(this_date)
           dates.append(this_date+dt.timedelta(minutes=30))
           this_date=this_date+dt.timedelta(hours=1)
        
        imerg=read_gpm(times=dates, intersection=intersection)
        # take mean of these times
        imerg=imerg.collapsed('time', iris.analysis.MEAN)
        this_imerg_cubes.append(imerg)
        this_day=this_day+dt.timedelta(days=1)
        
    this_imerg_cubes=this_imerg_cubes.merge_cube()
    return this_imerg_cubes
