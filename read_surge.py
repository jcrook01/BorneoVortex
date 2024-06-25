import datetime as dt
import iris
import iris.cube
import os.path as path
import numpy as np

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

# if start_date and end_date are not within one season there will be a gap in the dates as we only have Jan-March and Oct-Dec
def read_surge(start_date, end_date):
    indir='/home/users/jcrook001/FORSEA/surge_data/'
    first_month=dt.datetime(start_date.year,start_date.month,1)
    last_month=dt.datetime(end_date.year,end_date.month,1)
    this_date=first_month
    is_cross=[]
    is_merid=[]
    is_easterly=[]
    surge_times=[]
    while this_date <= last_month:
        cross_surge_filename=this_date.strftime(indir+'cross_surge_%Y%m.nc')
        merid_surge_filename=this_date.strftime(indir+'merid_surge_%Y%m.nc')
        east_surge_filename=this_date.strftime(indir+'easterly_surge_%Y%m.nc')
        print('reading surge for ',this_date.year, this_date.month)
        if path.isfile(cross_surge_filename)==True:
            cross_surge=iris.load_cube(cross_surge_filename)
            this_is_cross=cross_surge.coord('cross_equatorial_northerly_surge').points
            time_coord = cross_surge.coord('time')
            time_unit = time_coord.units
            this_surge_times=time_unit.num2date(time_coord.points)
            if len(is_cross)==0:
                is_cross=this_is_cross
                surge_times=this_surge_times
            else:
                is_cross=np.append(is_cross, this_is_cross)
                surge_times=np.append(surge_times, this_surge_times)
        if path.isfile(merid_surge_filename)==True:
            merid_surge=iris.load_cube(merid_surge_filename)
            this_is_merid=merid_surge.coord('meridional_surge').points
            if len(is_merid)==0:
                is_merid=this_is_merid
            else:
                is_merid=np.append(is_merid, this_is_merid)
            moderate=merid_surge.coord('moderate_meridional_surge').points
            ix=np.where(moderate==1)
            is_merid[ix]=2
            strong=merid_surge.coord('strong_meridional_surge').points
            ix=np.where(strong==1)
            is_merid[ix]=3
        if path.isfile(east_surge_filename)==True:
            easterly_surge=iris.load_cube(east_surge_filename)
            this_is_easterly=easterly_surge.coord('easterly_surge').points
            if len(is_easterly)==0:
                is_easterly=this_is_easterly
            else:
                is_easterly=np.append(is_easterly,this_is_easterly)
            moderate=easterly_surge.coord('moderate_easterly_surge').points
            ix=np.where(moderate==1)
            is_easterly[ix]=2
            strong=easterly_surge.coord('strong_easterly_surge').points
            ix=np.where(strong==1)
            is_easterly[ix]=3
        this_date=this_surge_times[-1].replace(hour=0)+dt.timedelta(days=1)
        if this_date.month==4: # skip to October
            this_date=dt.datetime(this_date.year,10,1)
        
    return is_cross, is_merid, is_easterly, surge_times

def read_surge_indices(start_date, end_date):
    is_cross, is_merid, is_easterly, surge_times=read_surge(start_date, end_date)
    surge_indices=get_surge_indices(is_cross,is_merid,is_easterly)
            
    return surge_indices, surge_times