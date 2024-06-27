import numpy as np
import datetime as dt
import pdb
from tropical_cyclones import *
import iris
import iris.analysis
import numpy.ma as ma

# this reads the Borneo Vortex track files that Kevin has produced for us with the tracks from era5

MIN_LON_TO_PLOT=75
MAX_LON_TO_PLOT=145
class track_info:

    def __init__(self):
        self.track_id=-1
        self.nt=0

    def create(self, track_id, nt):
        self.track_id=track_id
        self.is_tc=np.zeros(nt,int)-1
        self.nt=nt
        self.track_times=np.asarray([dt.datetime(1900,1,1) for i in range(nt)])
        # all following data has longitude, latitude and value
        self.vort_data=np.zeros((3,nt))
        self.mslp_data=np.zeros((4,nt)) # use 4th element to indicate minimum in MSLP was not found so vort centre used
        self.max_windspeed_data=np.zeros((3,nt))
        self.max_windspeed10m_data=np.zeros((3,nt))
        self.imerg=np.zeros(self.nt)-1

    # update the current_lon_range and current_lat_range with the range of lons and lats that this BV covers
    # if this BV goes outside current ranges 
    def get_lon_lat_range(self, current_lon_range, current_lat_range):
        min_lon=np.amin(self.vort_data[0,:])
        max_lon=np.amax(self.vort_data[0,:])
        min_lat=np.amin(self.vort_data[1,:])
        max_lat=np.amax(self.vort_data[1,:])
        #print(self.track_id, min_lon, max_lon, min_lat, max_lat)
        if current_lon_range[0]>min_lon:
            if min_lon<MIN_LON_TO_PLOT:
                #print('setting min_lon to MIN_LON_TO_PLOT')
                min_lon=MIN_LON_TO_PLOT
            current_lon_range[0]=min_lon-0.5
        if current_lon_range[1]<max_lon:
            if max_lon>MAX_LON_TO_PLOT:
                #print('setting max_lon to MAX_LON_TO_PLOT')
                max_lon=MAX_LON_TO_PLOT
            current_lon_range[1]=max_lon+0.5
        if current_lat_range[0]>min_lat:
            current_lat_range[0]=min_lat-0.5
            #print('setting min lat to', min_lat)
        if current_lat_range[1]<max_lat:
            current_lat_range[1]=max_lat+0.5

    #----------------------------------------------------
    # methods to check for existence in region and time
    #-----------------------------------------------------
    def exists_in_region(self, lon_range, lat_range):
        lons=self.vort_data[0,:]
        lats=self.vort_data[1,:]

        ix=np.where((lons>=lon_range[0]) & (lons<=lon_range[1]) & (lats>=lat_range[0]) & (lats<=lat_range[1]))
        return ix

    def exists_in_region_mostly(self, lon_range, lat_range, fraction):
        lons=self.vort_data[0,:]
        lats=self.vort_data[1,:]

        ix=np.where((lons>=lon_range[0]) & (lons<=lon_range[1]) & (lats>=lat_range[0]) & (lats<=lat_range[1]))
        nin_reg=len(ix[0])
        return nin_reg>=self.nt*fraction

    def exists_after_time(self, this_date):
        after=np.asarray([date>this_date for date in self.track_times])
        
        ix=np.where(after)
        return len(ix[0])>0

    def exists_in_month(self, year, month):
        my_years=np.asarray([date.year for date in self.track_times])
        my_months=np.asarray([date.month for date in self.track_times])
        
        ix=np.where((my_years==year) & (my_months==month))
        return len(ix[0])>0

    def exists_in_month_in_region(self, year, month, lon_range, lat_range):
        my_years=np.asarray([date.year for date in self.track_times])
        my_months=np.asarray([date.month for date in self.track_times])
        
        ix=np.where((my_years==year) & (my_months==month))
        this_lons=self.vort_data[0,ix[0]]
        this_lats=self.vort_data[1,ix[0]]
        ix=np.where((this_lons>=lon_range[0]) & (this_lons<=lon_range[1]) & (this_lats>=lat_range[0]) & (this_lats<=lat_range[1]))
        return len(ix[0])>0

    # days is an array of dates where we just need to check day and month
    def exists_on_days(self, days):
        my_months=np.asarray([date.month for date in self.track_times])
        my_days=np.asarray([date.day for date in self.track_times])
        exists=False
        for d in days:
            ix=np.where((my_days==d.day) & (my_months==d.month))
            if len(ix[0])>0:
                exists=True

        return exists

    # days is an array of dates where we need to check if the track existed at that time in the region and return the indices to the track when true
    def get_matching_times_in_region(self, days, lon_range, lat_range):

        my_years=np.asarray([date.year for date in self.track_times])
        my_months=np.asarray([date.month for date in self.track_times])
        my_days=np.asarray([date.day for date in self.track_times])
        exists=np.zeros(self.nt, int)
        for d in days:
            ix=np.where((my_days==d.day) & (my_months==d.month) & (my_years==d.year))
            if len(ix[0])>0:
                my_lons=self.vort_data[0,ix[0]]
                my_lats=self.vort_data[1,ix[0]]
                ix_reg=np.where((my_lons-lon_range[0]>=0) & (lon_range[1]-my_lons>=0) & (my_lats-lat_range[0]>=0) & (lat_range[1]-my_lats>=0))
                if len(ix_reg[0])>0:
                    exists[ix[0][ix_reg]]=True

        return np.where(exists==True)

    # days is an array of dates where we just need to check day and month, unless check_year=True
    def exists_on_days_in_region(self, days, lon_range, lat_range, check_year=False):

        if check_year==True:
            my_years=np.asarray([date.year for date in self.track_times])
        my_months=np.asarray([date.month for date in self.track_times])
        my_days=np.asarray([date.day for date in self.track_times])
        exists=False
        for d in days:
            if check_year==True:
                ix=np.where((my_days==d.day) & (my_months==d.month) & (my_years==d.year))
            else:
                ix=np.where((my_days==d.day) & (my_months==d.month))
            if len(ix[0])>0:
                my_lons=self.vort_data[0,ix[0]]
                my_lats=self.vort_data[1,ix[0]]
                ix_reg=np.where((my_lons-lon_range[0]>=0) & (lon_range[1]-my_lons>=0) & (my_lats-lat_range[0]>=0) & (lat_range[1]-my_lats>=0))
                if len(ix_reg[0])>0:
                    exists=True

        return exists

    def reaches_vorticity_threshold(self,threshold):
        ix=np.where(self.vort_data[2,:]>=threshold)
        return len(ix[0])>0

    # find whether there exists a TC that occurs at the same time and close to this BV
    # tc is a list of Tropical_Cyclones objects
    def matches_tc(self,tc):
        is_tc=np.zeros(self.nt, int)
        hours_start=np.asarray([(date-tc.dates[0]).days*24+(date-tc.dates[0]).seconds/3600 for date in self.track_times])
        hours_end=np.asarray([(date-tc.dates[-1]).days*24+(date-tc.dates[-1]).seconds/3600 for date in self.track_times])
        ix=np.where((hours_start>=0) & (hours_end<=0))
        if len(ix[0])>0:
            #print('track times',self.track_times[0], self.track_times[-1])
            #print(tc.name, 'tc times',tc.dates[0], tc.dates[-1])

            # this track does overlap with times of TC now check position
            lat_diff=np.zeros(len(ix[0]))
            lon_diff=np.zeros(len(ix[0]))
            lats=self.vort_data[1,ix[0]]
            lons=self.vort_data[0,ix[0]]
            for i in range(len(ix[0])):
                tc_lat, tc_lon=tc.get_tc_centre_at_date(self.track_times[ix[0][i]])
                lat_diff[i]=abs(lats[i]-tc_lat)
                lon_diff[i]=abs(lons[i]-tc_lon)
            pos_thresh=3 # BVs at 850 hPa can be a little distance away from surface position of TC in IBTrACS
            #pdb.set_trace()
            if np.mean(lat_diff)<pos_thresh and np.mean(lon_diff)<pos_thresh:
               is_tc[ix[0][i]]=1
            if np.mean(lat_diff)<2 and np.mean(lon_diff)<2:
                print('TC does exist at same time', self.track_times[ix[0][0]],self.track_times[ix[0][-1]], tc.name, 'lon diff',np.mean(lon_diff), 'lat diff',np.mean(lat_diff), is_tc)
        return is_tc

    def get_tc(self):
        tc=-1
        ix=np.where(self.is_tc>=0)
        if len(ix[0])>0:
            tc=self.is_tc[ix[0][0]]
        return tc
        
    # get data averaged over a circle of radius dlonlat centred on the track centre at time tix
    def get_mean_centred_data(self,data,lons,lats,dlonlat,tix):
        cenlon=self.vort_data[0,tix]
        cenlat=self.vort_data[1,tix]
        if cenlon<np.amin(lons) or cenlon>np.amax(lons) or cenlat<np.amin(lats) or cenlat>np.amax(lats):
            print('centre outside data region')
            pdb.set_trace()
            mean_data=None
        else:
            X, Y = np.meshgrid(lons, lats)
            dist = np.hypot(X-cenlon, Y-cenlat)
            #pdb.set_trace()
            mask=np.zeros_like(dist)
            ix=np.where(dist>dlonlat)
            mask[ix]=True
            if hasattr(data, 'mask'):
                data = ma.masked_array(data.data, mask=mask)
            else:
                data = ma.masked_array(data, mask=mask)
            mean_data=np.mean(data)
        return mean_data

    # vorticity is 6 hourly data at a single pressure level that matches the BV track level
    # this picks out the relevant times for this vortex track and takes the mean over a +/- 4 degree circle
    # centred on the track centre
    def add_vorticity(self, vorticity, lons, lats, dates, shift_minutes=0):
        dlonlat=4
        if hasattr(self, 'mean_vort')==False:
            self.mean_vort=np.zeros(self.nt)-1

        for t in range(self.nt):
            tix=np.where(dates==self.track_times[t]-dt.timedelta(minutes=shift_minutes))
            if len(tix[0])>0:
                this_vort=vorticity[tix[0][0]]
                vort_mean_data=self.get_mean_centred_data(this_vort,lons,lats,dlonlat,t)
                self.mean_vort[t]=vort_mean_data
            else:
                print('no matching time in vorticity')
                pdb.set_trace()
                
    # imerg_cube is a cube of hourly mean imerg every hour
    # this picks out the relevant times for this vortex track and takes the mean over a +/- 4 degree circle
    # centred on the track centre
    # for imerg we get the mean from half hourly data centred on the hours of the track (although the time comes
    # out as 15 minutes after the track
    # but for UM simulations we have hourly mean data that is not centered on the track times so times are 30 minutes out
    def add_rain(self, imerg_cube, shift_minutes=15):
        dlonlat=4
        lons_p=imerg_cube.coord('longitude').points
        lats_p=imerg_cube.coord('latitude').points
        imerg_t_coord=imerg_cube.coord('time')
        imerg_dates = imerg_t_coord.units.num2date(imerg_t_coord.points)
        imerg_start=imerg_dates[0]+dt.timedelta(minutes=shift_minutes)
        imerg_end=imerg_dates[-1]+dt.timedelta(minutes=shift_minutes)
        hours_start=np.asarray([(date-imerg_start).days*24+(date-imerg_start).seconds/3600 for date in self.track_times])
        hours_end=np.asarray([(date-imerg_end).days*24+(date-imerg_end).seconds/3600 for date in self.track_times])
        ix=np.where((hours_start>=0) & (hours_end<=0))
        if len(ix[0])>0:
          for t in range(self.nt):
            ix_imerg=np.where(imerg_dates==self.track_times[t]-dt.timedelta(minutes=shift_minutes)) # imerg dates are 15 mins out
            if len(ix_imerg[0])>0:
                this_imerg=imerg_cube[ix_imerg[0][0]]
                imerg_mean_data=self.get_mean_centred_data(this_imerg.data,lons_p,lats_p,dlonlat,t)
                #print(self.track_id, 'imerg data added at tix', imerg_mean.data, t)
                self.imerg[t]=imerg_mean_data
        else:
            print(self.track_id, 'no matching dates to imerg')
            pdb.set_trace()

# read all tracks from start_date onwards
def read_kevin_bv_tracks(track_file, start_date, my_tcs=[]):

    print("reading "+track_file)
    file = open(track_file, "r")
    l=0 # line number
    found_track=0
    nt=0
    this_id=-1
    this_tix=-1
    bv_tracks=[]
    ntracks=0
    ntcs=len(my_tcs)
    for line in file:
        words=line.split()
        if words[0] == 'TRACK_ID':
            start_time=dt.datetime.strptime(words[3],"%Y%m%d%H")
            if (start_time>=start_date):
                found_track=l # this is the line number of track id
                this_id=int(words[1])
                #print('line', l, ' found track', this_id, start_time)
                bv_tracks.append(track_info())
                ntracks=ntracks+1

        # next line has number of observations
        if found_track > 0 and l == found_track+1:
            words=line.split()
            nwords=len(words)
            if words[0] == 'POINT_NUM':
                nt=int(words[nwords-1])
                this_tix=0
                #print('line', l,' creating track', this_id, nt)
                bv_tracks[ntracks-1].create(this_id,nt)

       # next line has first useful stuff on it
        if found_track > 0 and l > found_track+1 and l <= found_track+1+nt:
            words=line.split()
            nwords=len(words)
            if nwords < 22:
                pdb.set_trace()
            bv_tracks[ntracks-1].track_times[this_tix]=dt.datetime.strptime(words[0],"%Y%m%d%H")
            #print('line',l,' adding data to track', bv_tracks[ntracks-1].track_id, bv_tracks[ntracks-1].track_times[this_tix])
            bv_tracks[ntracks-1].vort_data[0,this_tix]=float(words[1])
            bv_tracks[ntracks-1].vort_data[1,this_tix]=float(words[2])
            bv_tracks[ntracks-1].vort_data[2,this_tix]=float(words[3])
            bv_tracks[ntracks-1].mslp_data[0,this_tix]=float(words[5])
            bv_tracks[ntracks-1].mslp_data[1,this_tix]=float(words[7])
            bv_tracks[ntracks-1].mslp_data[2,this_tix]=float(words[9])
            if bv_tracks[ntracks-1].mslp_data[0,this_tix]>=1.0e20:
                bv_tracks[ntracks-1].mslp_data[0,this_tix]=bv_tracks[ntracks-1].vort_data[0,this_tix]
                bv_tracks[ntracks-1].mslp_data[3,this_tix]=1 # set flag to indicate no mslp min found
            if bv_tracks[ntracks-1].mslp_data[1,this_tix]>=1.0e20:
                bv_tracks[ntracks-1].mslp_data[1,this_tix]=bv_tracks[ntracks-1].vort_data[1,this_tix]
                bv_tracks[ntracks-1].mslp_data[3,this_tix]=1 # set flag to indicate no mslp min found

            bv_tracks[ntracks-1].max_windspeed_data[0,this_tix]=float(words[11])
            bv_tracks[ntracks-1].max_windspeed_data[1,this_tix]=float(words[13])
            bv_tracks[ntracks-1].max_windspeed_data[2,this_tix]=float(words[15])
            bv_tracks[ntracks-1].max_windspeed10m_data[0,this_tix]=float(words[17])
            bv_tracks[ntracks-1].max_windspeed10m_data[1,this_tix]=float(words[19])
            bv_tracks[ntracks-1].max_windspeed10m_data[2,this_tix]=float(words[21])
            this_tix=this_tix+1
            if this_tix==nt:
                #print('last time found for track', bv_tracks[ntracks-1].track_id)
                for t in range(ntcs):
                    if bv_tracks[ntracks-1].matches_tc(my_tcs[t]):
                        bv_tracks[ntracks-1].is_tc=t
                        break
                found_track=0
                nt=0
                this_tix=-1
			
        l=l+1

    return np.asarray(bv_tracks)

