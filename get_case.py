from read_kevin_bv_tracks import *
import datetime as dt
from get_era5_on_plevs import *
from get_below_ground_mask import *
from get_imerg import *
from set_latlon_bounds import *
from pathlib import Path
import pdb

# get the matching 925 track at each time and also find lat and lon range of 850 and 925 tracks
# track_850 is the bv_track object for the case we are interested in
# tracks_925 is a list of all the 925 hPa bv_track objects that matched (i.e. overlapped) with the 850 track at some point
def get_case_track_data(track_850, tracks_925):

    track_lats850=track_850.vort_data[1,:]
    track_lons850=track_850.vort_data[0,:]
    min_lat=-10
    max_lat=15
    min_lon=95
    max_lon=125
    # use extended area for u and v so we can determine vortex height and tilt
    ext_min_lat=min_lat
    ext_max_lat=max_lat
    ext_min_lon=min_lon
    ext_max_lon=max_lon
    min_lat850=np.amin(track_lats850)
    if min_lat850<min_lat:
        min_lat=min_lat850-1
        print('min lat extended by 850 to ', min_lat)
    if min_lat850-5<ext_min_lat:
        ext_min_lat=min_lat850-5
    max_lat850=np.amax(track_lats850)
    if max_lat850>max_lat:
        max_lat=max_lat850+1
        print('max lat extended by 850 to ', max_lat)
    if max_lat850+5>ext_max_lat:
        ext_max_lat=max_lat850+5
    min_lon850=np.amin(track_lons850)
    if min_lon850<min_lon:
        min_lon=min_lon850-1
        print('min lon extended by 850 to ', min_lon)
    if min_lon850-5<ext_min_lon:
        ext_min_lon=min_lon850-5
    max_lon850=np.amax(track_lons850)
    if max_lon850>max_lon:
        max_lon=max_lon850+1
        print('max lon extended by 850 to ', max_lon)
    if max_lon850+5>ext_max_lon:
        ext_max_lon=max_lon850+5

    # find matching 925 track at each time
    ntracks925=len(tracks_925)
    best_match=np.zeros(track_850.nt, int)-1
    tix925=np.zeros(track_850.nt, int)-1
    track_lats925=np.zeros(track_850.nt)
    track_lons925=np.zeros(track_850.nt)
    for t in range(track_850.nt):
        min_lon_diff=2
        min_lat_diff=2
        for j in range(ntracks925):
            exists=np.asarray([time == track_850.track_times[t] for time in tracks_925[j].track_times])
            ix_match=np.where(exists)
            if len(ix_match[0])>0:
                if len(ix_match[0])>1:
                    pdb.set_trace()
                lon925=tracks_925[j].vort_data[0,ix_match[0][0]]
                lat925=tracks_925[j].vort_data[1,ix_match[0][0]]
                lon_diff=abs(track_lons850[t]-lon925)
                lat_diff=abs(track_lats850[t]-lat925)
                if lon_diff<min_lon_diff and lat_diff<min_lat_diff:
                    # this is the best match so far
                    best_match[t]=j
                    min_lon_diff=lon_diff
                    min_lat_diff=lat_diff
                    track_lats925[t]=lat925
                    track_lons925[t]=lon925
                    tix925[t]=ix_match[0][0]
    ix=np.where(best_match>=0)
    max_lat925=np.amax(track_lats925[ix[0]])
    if max_lat925>max_lat:
        max_lat=max_lat925
        print('max lat extended by 925 to ', max_lat)
    max_lon925=np.amax(track_lons925[ix[0]])
    if max_lon925>max_lon:
        max_lon=max_lon925
        print('max lon extended by 925 to ', max_lon)
    min_lon925=np.amin(track_lons925[ix[0]])
    if min_lon925<min_lon:
        min_lon=min_lon925
        print('min lon extended by 925 to ', min_lon)
    if max_lat>20:
        max_lat=20
    print('lon lat range', min_lon, max_lon, min_lat, max_lat)
    return track_lats925, track_lons925, best_match, tix925, min_lat, max_lat, min_lon, max_lon, ext_min_lat, ext_max_lat, ext_min_lon, ext_max_lon
    
def get_case_track925_rain_vort(nt, tracks925, best_match, tix925):
    rain925=np.zeros(nt)-1
    vort925=np.zeros(nt)-1
    for t in range(nt):
        j=best_match[t]
        if j>=0:
            rain925[t]=tracks925[j].imerg[tix925[t]]
            vort925[t]=tracks925[j].vort_data[2,tix925[t]]

    ix=np.where(best_match>=0)
    bad_ix=np.where(rain925[ix]<0)
    if len(bad_ix[0])>0:
        print('925 track rain', rain925[ix])
    return rain925, vort925

def get_common_stuff(bv_tracks850, bv_tracks925, case_id850):
    if case_id850==9636:
        case_ids925=[9239,10051]
    elif case_id850==215: # not in top since 2018 BVs as too far north for true BV
        case_ids925=[108, 530]
    elif case_id850==8371:
        case_ids925=[6537,8871]
    elif case_id850==7732:
        case_ids925=[7165, 7913]
    elif case_id850==5655:
        case_ids925=[5260]
    elif case_id850==1025:
        case_ids925=[1029]
    elif case_id850==5020:  # not in our top since 2018 BVs
        case_ids925=[4717,5554]
    elif case_id850==10475:
        case_ids925=[11484,11656]
    elif case_id850==10158:
        case_ids925=[9511]
    elif case_id850==6970:
        case_ids925=[5915] 
    elif case_id850==4966:
        case_ids925=[4905,4955]
    elif case_id850==8200:
        case_ids925=[7426]
    elif case_id850==9138: # Natasha's case
        case_ids925=[8423,8588]
    elif case_id850==6032:
        case_ids925=[5260, 5915]
    elif case_id850==7032:
        case_ids925=[6537]
    elif case_id850==16941:
        case_ids925=[16382, 16947]
    elif case_id850==10616:  # not in top since 2018 BVs
        case_ids925=[9990, 10859]
    elif case_id850==4551:
        case_ids925=[4266] 
    elif case_id850==10832:
        case_ids925=[10074, 10129] 
    elif case_id850==7012:
        case_ids925=[6516, 6637]
    elif case_id850==6328:
        case_ids925=[5944,6407,7354]
    elif case_id850==1706: # Sam's case
        case_ids925=[1612]
    elif case_id850==5863:
        case_ids925=[5450,5642,6161]
    elif case_id850==15722:
        case_ids925=[15051]
    elif case_id850==11557:
        case_ids925=[10505,11557]
    elif case_id850==4646: # not in top since 2018 BVs
        case_ids925=[4447,4516,4392] 
    elif case_id850==7716:
        case_ids925=[6637]
    else:
        print('Invalid case')
        pdb.set_trace()
        
    ids850=np.asarray([track.track_id for track in bv_tracks850])
    ix850=np.where(ids850==case_id850)
    track850=bv_tracks850[ix850[0][0]]
    start_date=track850.track_times[0]
    end_date=track850.track_times[-1]

    ids925=np.asarray([track.track_id for track in bv_tracks925])
    tracks925=[]
    for id in case_ids925:
        ix925=np.where(ids925==id)
        if len(ix925[0])==0:
            print(id, '925 track not found')
        else:
            if len(ix925[0])>1:
                print('too many matching 925 tracks')
                pdb.set_trace()
            tracks925.append(bv_tracks925[ix925[0][0]])
        
    track_lats925, track_lons925, best_match, tix925, min_lat, max_lat, min_lon, max_lon, ext_min_lat, ext_max_lat, ext_min_lon, ext_max_lon=get_case_track_data(track850, tracks925)

    outdir='/home/users/jcrook001/FORSEA/Cases/'+start_date.strftime('%Y%m%d-')+end_date.strftime('%Y%m%d/')
    Path(outdir).mkdir(parents=True, exist_ok=True)

    print('case id', case_id850, outdir)

    return track850, track_lats925, track_lons925, best_match, tix925, min_lat, max_lat, min_lon, max_lon, ext_min_lat, ext_max_lat, ext_min_lon, ext_max_lon, outdir
    
from get_theta_e import *
def get_case_for_composite(bv_tracks850, bv_tracks925, case_id850):
    track850, track_lats925, track_lons925, best_match, tix925, min_lat, max_lat, min_lon, max_lon, ext_min_lat, ext_max_lat, ext_min_lon, ext_max_lon, outdir=get_common_stuff(bv_tracks850, bv_tracks925, case_id850)
    intersection = {'latitude': [min_lat,max_lat], 'longitude': [min_lon,max_lon]}
    ext_intersection = {'latitude': [min_lat-10,max_lat+10], 'longitude': [min_lon-10,max_lon+10]}
    start_date=track850.track_times[0]
    end_date=track850.track_times[-1]
    track_lons=track850.vort_data[0,:]
    track_lats=track850.vort_data[1,:]

    composite_ew_fname=outdir+'composite_ew.npy'
    composite_ns_fname=outdir+'composite_ns.npy'
    composite_dlons_fname=outdir+'composite_dlons.npy'
    composite_dlats_fname=outdir+'composite_dlats.npy'
    if path.isfile(composite_ew_fname)==True and path.isfile(composite_ns_fname)==True:
        print('reading', composite_ew_fname)
        composites_ew=np.load(composite_ew_fname)
        print('reading', composite_ns_fname)
        composites_ns=np.load(composite_ns_fname)
        delta_lons=np.load(composite_dlons_fname)
        delta_lats=np.load(composite_dlats_fname)
        era5_u=get_era5_on_plevs(start_date, start_date, ext_intersection, 'u', 'eastward_wind', dhours=6)
        pressure = era5_u.coord('pressure_level').points

        return composites_ew, composites_ns, delta_lons, delta_lats, pressure, outdir
        
    u_ix=0
    v_ix=1
    w_ix=2
    q_ix=3
    thetae_ix=4
    ncomp_ew=0
    ncomp_ns=0
    dlonlat=8
    wave_fbase=outdir+'wave_hovmoller_'
    kelvin_w2=np.load(wave_fbase+'kelvin.npy',allow_pickle=True)
    wmrg_w2=np.load(wave_fbase+'wmrg.npy',allow_pickle=True)
    r1_w2=np.load(wave_fbase+'r1.npy',allow_pickle=True)
    lons_wave=np.load(wave_fbase+'wave_lons.npy',allow_pickle=True)
    nt_wave=len(kelvin_w2[:,0])
    nwaves=3
    nwave_phases=5
    for t in range(track850.nt):
        this_date=track850.track_times[t]
        era5_u=get_era5_on_plevs(this_date, this_date, ext_intersection, 'u', 'eastward_wind', dhours=6)
        era5_v=get_era5_on_plevs(this_date, this_date, ext_intersection, 'v', 'northward_wind', dhours=6)
        era5_w=get_era5_on_plevs(this_date, this_date, ext_intersection, 'omega', 'lagrangian_tendency_of_air_pressure', dhours=6)
        era5_q=get_era5_on_plevs(this_date, this_date, ext_intersection, 'q', 'specific_humidity', dhours=6)
        era5_q.convert_units('g/kg')
        era5_T=get_era5_on_plevs(this_date, this_date, ext_intersection, 'T', 'air_temperature', dhours=6)

        era5_thetae=get_equivalent_potential_temperature(era5_T, era5_q)
        if t==0:
            lons_u = era5_u.coord('longitude').points
            lats_u = era5_u.coord('latitude').points
            lons = era5_q.coord('longitude').points
            lats = era5_q.coord('latitude').points
            ix=np.where(abs(lons_u-lons)>0.001)
            if len(ix[0])>0:
                pdb.set_trace()
            ix=np.where(abs(lats_u-lats)>0.001)
            if len(ix[0])>0:
                pdb.set_trace()
            pressure = era5_u.coord('pressure_level').points
            nplev=len(pressure)
            print('lon lat range read:', np.amin(lons), np.amax(lons), np.amin(lats), np.amax(lats))
        
        # get cross sections centred on the track centre
        ix_lat=np.where(abs(lats-track_lats[t])<=1)
        ix_lon=np.where(abs(lons-track_lons[t])<=dlonlat)
        this_delta_lons=lons[ix_lon[0]]-track_lons[t]
        if t==0:
            ndlon=len(this_delta_lons)
            delta_lons=np.zeros((ndlon,track850.nt))+np.nan
            composites_ew=np.zeros((nplev,ndlon,5))
            composites_ew_by_wave=np.zeros((nplev,ndlon,5,nwaves,nwave_phases))
        if len(ix_lon[0]) != ndlon:
            print('wrong number of longitudes so missing this time', track_lons[t])
            pdb.set_trace()
        else:
            delta_lons[:,t]=this_delta_lons
            this_v_cross_ew=np.mean(era5_v[:,ix_lat[0],:].data, axis=1)[:,ix_lon[0]]
            composites_ew[:,:,v_ix]=composites_ew[:,:,v_ix]+this_v_cross_ew
            this_u_cross_ew=np.mean(era5_u[:,ix_lat[0],:].data, axis=1)[:,ix_lon[0]]
            composites_ew[:,:,u_ix]=composites_ew[:,:,u_ix]+this_u_cross_ew
            this_w_cross_ew=np.mean(era5_w[:,ix_lat[0],:].data, axis=1)[:,ix_lon[0]]
            composites_ew[:,:,w_ix]=composites_ew[:,:,w_ix]+this_w_cross_ew
            this_thetae_cross_ew=np.mean(era5_thetae[:,ix_lat[0],:].data, axis=1)[:,ix_lon[0]]
            composites_ew[:,:,thetae_ix]=composites_ew[:,:,thetae_ix]+this_thetae_cross_ew
            this_q_cross_ew=np.mean(era5_q[:,ix_lat[0],:].data, axis=1)[:,ix_lon[0]]
            composites_ew[:,:,q_ix]=composites_ew[:,:,q_ix]+this_q_cross_ew
            ncomp_ew=ncomp_ew+1
            
        ix_lon=np.where(abs(lons-track_lons[t])<=1)
        ix_lat=np.where(abs(lats-track_lats[t])<=dlonlat)
        this_delta_lats=lats[ix_lat[0]]-track_lats[t]
        if t==0:
            ndlat=len(this_delta_lats)
            delta_lats=np.zeros((ndlat,track850.nt))+np.nan
            composites_ns=np.zeros((nplev,ndlat,5))
        if len(ix_lat[0]) != ndlat:
            print('wrong number of latitudes so missing this time', track_lats[t])
            pdb.set_trace()
        else:
            delta_lats[:,t]=this_delta_lats

            this_u_cross_ns=np.mean(era5_u[:,:,ix_lon[0]].data, axis=2)[:,ix_lat[0]]
            composites_ns[:,:,u_ix]=composites_ns[:,:,u_ix]+this_u_cross_ns
            this_v_cross_ns=np.mean(era5_v[:,:,ix_lon[0]].data, axis=2)[:,ix_lat[0]]
            composites_ns[:,:,v_ix]=composites_ns[:,:,v_ix]+this_v_cross_ns
            this_w_cross_ns=np.mean(era5_w[:,:,ix_lon[0]].data, axis=2)[:,ix_lat[0]]
            composites_ns[:,:,w_ix]=composites_ns[:,:,w_ix]+this_w_cross_ns
            this_thetae_cross_ns=np.mean(era5_thetae[:,:,ix_lon[0]].data, axis=2)[:,ix_lat[0]]
            composites_ns[:,:,thetae_ix]=composites_ns[:,:,thetae_ix]+this_thetae_cross_ns
            this_q_cross_ns=np.mean(era5_q[:,:,ix_lon[0]].data, axis=2)[:,ix_lat[0]]
            composites_ns[:,:,q_ix]=composites_ns[:,:,q_ix]+this_q_cross_ns
            ncomp_ns=ncomp_ns+1
            
    composites_ew=composites_ew/ncomp_ew
    composites_ns=composites_ns/ncomp_ns
    np.save(composite_ew_fname, composites_ew)
    np.save(composite_ns_fname, composites_ns)
    delta_lons=np.nanmean(delta_lons,axis=1)
    delta_lats=np.nanmean(delta_lats,axis=1)
    np.save(composite_dlons_fname, delta_lons)
    np.save(composite_dlats_fname, delta_lats)
    return  composites_ew, composites_ns, delta_lons, delta_lats, pressure, outdir
    
def get_case(bv_tracks850, bv_tracks925, case_id850):
        
    track850, track_lats925, track_lons925, best_match, tix925, min_lat, max_lat, min_lon, max_lon, ext_min_lat, ext_max_lat, ext_min_lon, ext_max_lon, outdir=get_common_stuff(bv_tracks850, bv_tracks925, case_id850)

    intersection = {'latitude': [min_lat,max_lat], 'longitude': [min_lon,max_lon]}
    ext_intersection = {'latitude': [ext_min_lat,ext_max_lat], 'longitude': [ext_min_lon,ext_max_lon]}
    if max_lat>20:
        print('max lat', max_lat, 'cannot get orog above 20N')
        pdb.set_trace()


    fname=outdir+'imerg_hourly.nc'
    if path.isfile(fname)==True:
        print('reading', fname)
        imerg=iris.load_cube(fname,'pcp')
    else:
        imerg=get_hourly_imerg(start_date, end_date, intersection)
        set_latlon_bounds(imerg) # need to set this before doing mean over area in add_rain
        iris.save(imerg, fname)

    print('adding rain to tracks')
    track850.add_rain(imerg)
    #for track in tracks925:
    #    print(track.track_id, 'adding rain')
    #    track.add_rain(imerg)
    #rain925, vort925=get_case_track925_rain_vort(track850.nt, tracks925, best_match, tix925)
    orog=get_orography(intersection)    
    below_mask=get_era5_below_ground_mask(start_date, end_date, intersection, orog)
    era5_u=get_era5_on_plevs(start_date, end_date, ext_intersection, 'u', 'eastward_wind', dhours=6)
    era5_v=get_era5_on_plevs(start_date, end_date, ext_intersection, 'v', 'northward_wind', dhours=6)
    era5_w=get_era5_on_plevs(start_date, end_date, intersection, 'omega', 'lagrangian_tendency_of_air_pressure', dhours=6)
    era5_u925=get_era5_925(start_date, end_date, intersection, 'u', 'eastward_wind', dhours=6)
    era5_v925=get_era5_925(start_date, end_date, intersection, 'v', 'northward_wind', dhours=6)
    era5_w925=get_era5_925(start_date, end_date, intersection, 'omega', 'lagrangian_tendency_of_air_pressure', dhours=6)
    era5_q=None #get_era5_on_plevs(start_date, end_date, intersection, 'q', 'specific_humidity', dhours=6)
    #era5_q.convert_units('g/kg')

    era5_T=get_era5_on_plevs(start_date, end_date, intersection, 'T', 'air_temperature', dhours=6)

    return  track850, track_lats925, track_lons925, best_match, min_lat, max_lat, min_lon, max_lon, era5_u, era5_v, era5_w, era5_u925, era5_v925, era5_w925, era5_q, era5_T, below_mask, imerg, outdir
