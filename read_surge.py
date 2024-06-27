import datetime as dt
import iris
import iris.cube
import os.path as path
import numpy as np
import pdb

# lowercase means weak easterly, uppercase means moderate, uppercase italic means strong
Eweak='e'
Emod='E'
Estrong='$\mathbf{E}$'
Mweak='m'
Mmod='M'
Mstrong='$\mathbf{M}$'

surge_text=['N','X',Eweak,Eweak+'X',Emod,Emod+'X',Estrong,Estrong+'X', 
            Mweak,Mweak+'X',Mweak+Eweak,Mweak+Eweak+'X',Mweak+Emod,Mweak+Emod+'X',Mweak+Estrong, Mweak+Estrong+'X', 
            Mmod,Mmod+'X',Mmod+Eweak,Mmod+Eweak+'X',Mmod+Emod,Mmod+Emod+'X',Mmod+Estrong, Mmod+Estrong+'X',
            Mstrong,Mstrong+'X',Mstrong+Eweak,Mstrong+Eweak+'X',Mstrong+Emod,Mstrong+Emod+'X',Mstrong+Estrong, Mstrong+Estrong+'X']

def get_surge_indices(is_cross_surge,is_merid_surge,is_easterly_surge):
    surge_ix=np.asarray(is_cross_surge+(2*is_easterly_surge)+(8*is_merid_surge))
    return surge_ix

def read_surge_file(surge_file, coord_name, get_mod_strong, check_months, start_date, end_date):
    if path.isfile(surge_file)==True:
        this_surge=iris.load_cube(surge_file)
        time_coord = this_surge.coord('time')
        time_unit = time_coord.units
        this_surge_times=time_unit.num2date(time_coord.points)
        this_points=this_surge.coord(coord_name).points
        if get_mod_strong:
            moderate=this_surge.coord('moderate_'+coord_name).points
            strong=this_surge.coord('strong_'+coord_name).points
            ix=np.where(moderate==1)
            this_points[ix]=2
            ix=np.where(strong==1)
            this_points[ix]=3
        if check_months:
            months=np.asarray([date.month for date in this_surge_times])
            ix=np.where(((months<4) | (months>=10)) & ((this_surge_times>=start_date) & (this_surge_times<=end_date)))
            this_surge_times=this_surge_times[ix]
            this_points=this_points[ix]
        return this_points, this_surge_times
    else:
        print(surge_file, 'does not exist')
        return [], []
    
# Simon Peatman created yearly files of surge data in terramaris gws up to 2018
# I created surge data in monthly files thereafter
# We are only interested in Jan-March and Oct-Dec
# If start_date and end_date are not within one season there will be a gap in the
# dates as we only have Jan-March and Oct-Dec
# read_simon will for 2018 dates force to read Simon's data instead of mine - this is just for checking it matches
def read_surge(start_date, end_date, read_simon=False):
    indir='/gws/nopw/j04/forsea/users/jcrook/BorneoVortex/surge_data/'
    simon_dir='/gws/nopw/j04/terramaris/peatman/'
    simon_cross_dir=simon_dir+'cross_equatorial_northerly_surge/'
    simon_merid_dir=simon_dir+'meridional_surge/'
    simon_east_dir=simon_dir+'easterly_surge/'
    cross_coord='cross_equatorial_northerly_surge'
    merid_coord='meridional_surge'
    simon_merid_coord='cold_surge'
    east_coord='easterly_surge'
    first_month=dt.datetime(start_date.year,start_date.month,1)
    last_month=dt.datetime(end_date.year,end_date.month,1)
    this_date=first_month
    is_cross=[]
    is_merid=[]
    is_easterly=[]
    surge_times=[]
    while this_date <= last_month:
        if this_date.year<2018 or (this_date.year==2018 and read_simon):
            # read Simon's data
            cross_surge_filename = this_date.strftime(simon_cross_dir+'era5.v925.%Y.dmean.cross_equatorial_northerly_surge.nc')
            merid_surge_filename = this_date.strftime(simon_merid_dir+'ERA5_%Y_v_925hPa.daily.cold_surge.nc')
            east_surge_filename = this_date.strftime(simon_east_dir+'era5.u925.%Y.dmean.easterly_surge.nc')
            print('reading surge for ',this_date.year)
            check_months=True # we need to remove Apr-Sep
            this_merid_coord=simon_merid_coord
        else:
            cross_surge_filename=this_date.strftime(indir+'cross_surge_%Y%m.nc')
            merid_surge_filename=this_date.strftime(indir+'merid_surge_%Y%m.nc')
            east_surge_filename=this_date.strftime(indir+'easterly_surge_%Y%m.nc')
            print('reading surge for ',this_date.year, this_date.month)
            check_months=False
            this_merid_coord=merid_coord
        cross_surge, this_surge_times=read_surge_file(cross_surge_filename,cross_coord, False,check_months, first_month, end_date)
        merid_surge, this_surge_times=read_surge_file(merid_surge_filename,this_merid_coord,True,check_months,first_month, end_date)
        east_surge, this_surge_times=read_surge_file(east_surge_filename,east_coord,True,check_months,first_month, end_date)
        if len(surge_times)==0:
            is_cross=cross_surge
            surge_times=this_surge_times
            is_merid=merid_surge
            is_easterly=east_surge
        else:
            is_cross=np.append(is_cross, cross_surge)
            surge_times=np.append(surge_times, this_surge_times)
            is_merid=np.append(is_merid,merid_surge)
            is_easterly=np.append(is_easterly,east_surge)

        this_date=this_surge_times[-1].replace(hour=0)+dt.timedelta(days=1)
        if this_date.month==4: # skip to October
            this_date=dt.datetime(this_date.year,10,1)
        
    return is_cross, is_merid, is_easterly, surge_times

def read_surge_indices(start_date, end_date):
    is_cross, is_merid, is_easterly, surge_times=read_surge(start_date, end_date)
    surge_indices=get_surge_indices(is_cross,is_merid,is_easterly)
            
    return surge_indices, surge_times