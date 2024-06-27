import warnings
import numpy as np
import datetime
import string
import pdb

INVALID_LAT_LON=-1000
# this class holds the basic data for a tropcal cyclone
class Tropical_Cyclones():
    def __init__(self):
        self.name=''
        self.basins=[]
        self.subbasins=[]
        self.dates=[]
        self.lats=[]
        self.lons=[]
        self.usa_sshs=[]  # saffir simpson scale

    def set(self, name, basins, subbasins, dates, lats, lons, usa_sshs):
        self.name=name
        self.basins=basins
        self.subbasins=subbasins
        self.dates=dates
        self.lats=lats
        self.lons=lons
        self.usa_sshs=usa_sshs

    def get_description(self):
        description=self.name+' max cat='+str(np.amax(self.usa_sshs))+' '+self.dates[0].strftime('%Y%m%d %H:%M')+' to '+self.dates[-1].strftime('%Y%m%d %H:%M')
        return description

    # get TC centre where dates match date
    # if there is no date matching but date is within 3 hours of 2 tc dates find the centre interpolated between the 2 dates
    # if tc is never near date return INVALID_LAT_LON
    def get_tc_centre_at_date(self, date):

        tc_lat=INVALID_LAT_LON
        tc_lon=INVALID_LAT_LON
        ix=np.where(self.dates==date)
        if len(ix[0]>0):
            tc_lat=self.lats[ix[0][0]]
            tc_lon=self.lons[ix[0][0]]
        else:
            three_hours=3
            delta_t=self.dates-date
            delta_t_hours=np.asarray([dt.days*24+dt.seconds/3600 for dt in delta_t])
            ix=np.where(abs(delta_t_hours)<three_hours)
            nt=len(ix[0])
            if nt==2:
                if delta_t_hours[ix[0][0]]<0 and delta_t_hours[ix[0][1]]>0:
                    first_tc_lat=self.lats[ix[0][0]]
                    first_tc_lon=self.lons[ix[0][0]]
                    first_date=self.dates[ix[0][0]]
                    second_tc_lat=self.lats[ix[0][1]]
                    second_tc_lon=self.lons[ix[0][1]]
                    second_date=self.dates[ix[0][1]]
                    tc_delta_hours=(second_date-first_date).days*24+(second_date-first_date).seconds/3600
                    tc_lat=first_tc_lat+(second_tc_lat-first_tc_lat)*-1*delta_t_hours[ix[0][0]]/tc_delta_hours
                    tc_lon=first_tc_lon+(second_tc_lon-first_tc_lon)*-1*delta_t_hours[ix[0][0]]/tc_delta_hours
                    #print(first_date, first_tc_lat, first_tc_lon, second_date, second_tc_lat, second_tc_lon, date, tc_lat, tc_lon)
                #else:
                    #print(self.name, 'at delta hours',delta_t_hours[ix[0][0]], delta_t_hours[ix[0][1]], date)


        return tc_lat, tc_lon


if __name__ == '__main__':
    main()

