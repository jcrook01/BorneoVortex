import warnings
import os.path as path
import numpy as np
import datetime as dt
import pdb
from read_kevin_bv_tracks import *

# This is the code to read TrueBVs850.npy or TrueBVs850_since_2000_with_rain.npy which were
# created in Leeds as follows:
# Vortex must go through the Borneo Vortex region for at least 24 hours
# Vortex has to have max_vort reaching 3.5
# Vortex has to also be seen in 925 hPa tracks for some of its life
#
# The minimum criteria is also that it shouldn't be a TC but as I wanted to look at
# which BVs became TCs I did not remove the TCs at this point.
# Therefore once the BVs have been read we also must remove BVs that are TCs
# We often also specify it should be init in region

bv_region_lon_range=[100,117]
bv_region_lat_range=[-2.5,10]

def get_non_tcs(bv_tracks850, allow_becomes_tc):

    nbvs=len(bv_tracks850)
    is_tc=np.zeros(nbvs, bool)
    for i in range(nbvs):
        became_tc850=False
        is_tc850=False
        if bv_tracks850[i].get_tc()>=0:
            is_tc850=True
            if allow_becomes_tc:
                # has been matched to a TC at some point but did this happen after it was in BV region or before?
                tix_in_reg=bv_tracks850[i].exists_in_region(bv_region_lon_range, bv_region_lat_range)
                tix_tc=np.where(bv_tracks850[i].is_tc>=0)
                if tix_tc[0][0] > tix_in_reg[0][-1]:
                    print(bv_tracks850[i].track_id, 'becomes TC ', bv_tracks850[i].is_tc[tix_tc[0][0]])
                    became_tc850=True

        is_tc[i]=(is_tc850 & (became_tc850==False))

    ix=np.where(is_tc==False)
    bv_tracks850=bv_tracks850[ix]
    print(len(bv_tracks850), 'tracks not TCs')
    return bv_tracks850

def get_tracks_in_region(bv_tracks):
    in_reg=np.asarray([len(track.exists_in_region(bv_region_lon_range, bv_region_lat_range)[0])>3 for track in bv_tracks])
    ix=np.where(in_reg)
    print(len(ix[0]), 'in reg')
    bv_tracks=bv_tracks[ix]

    return bv_tracks

def get_subset_init_in_reg(bv_tracks):
    # check bv was not formed out of this region
    lon_range=[100,125]
    lat_range=[-5,12]
    out_reg=np.asarray([((bv.vort_data[0,0]<lon_range[0]) | (bv.vort_data[0,0]>lon_range[1]) | (bv.vort_data[1,0]<lat_range[0]) | (bv.vort_data[1,0]>lat_range[1])) for bv in bv_tracks])
    
    ix=np.where(out_reg==False)
    print(len(ix[0]), 'init in reg')

    return ix

def get_tracks_init_in_region(bv_tracks):

    ix=get_subset_init_in_reg(bv_tracks)
    bv_tracks=bv_tracks[ix]

    return bv_tracks

def get_high_vort_tracks(bv_tracks):
    vort_threshold=3.5
    is_bv=np.asarray([track.reaches_vorticity_threshold(vort_threshold) for track in bv_tracks])
    ix_bv=np.where(is_bv)
    bv_tracks=bv_tracks[ix_bv]
    print(len(bv_tracks),'high vort')
    return bv_tracks

def_start=dt.datetime(1900,1,1)
def get_true_bvs(start_date=def_start, init_in_reg=False):

    basedir='/gws/nopw/j04/forsea/users/jcrook/BorneoVortex/'
    trackdir=basedir+'BV_tracks/'
    true_file=trackdir+'TrueBVs850.npy'
    #if path.isfile(true_file)==True:
    try:
        bv_tracks850=np.load(true_file, allow_pickle=True)
    except OSError as err:
        print(true_file)
        raise err        


    if start_date>def_start:
        ntot_bvs=len(bv_tracks850)
        after=np.asarray([track.exists_after_time(start_date) for track in bv_tracks850])
        ix=np.where(after)
        bv_tracks850=bv_tracks850[ix]
        nbvs=len(bv_tracks850)
        print(ntot_bvs, nbvs, 'total bvs and bvs since ', start_date)
    if init_in_reg:
        bv_tracks850=get_tracks_init_in_region(bv_tracks850)

    return bv_tracks850

def get_true_bvs_with_rain(year, init_in_reg=False):

    # we have rain from oct 2000 so year must be no earlier than 2000 
    # this gets the true BV tracks with rain since year

    basedir='/gws/nopw/j04/forsea/users/jcrook/BorneoVortex/'
    trackdir=basedir+'BV_tracks/'
    track_file=trackdir+'TrueBVs850_since_2000_with_rain.npy'
    if year < 2000:
        print('invalid year')
        return []


    try:
        bv_tracks=np.load(track_file, allow_pickle=True)
    except OSError as err:
        print(track_file)
        raise err        


    if year> 2000:
        # remove earlier ones
        years=np.asarray([track.track_times[0].year for track in bv_tracks])
        ix=np.where(years>=year)
        # print(len(bv_tracks), len(ix[0]))
        bv_tracks=bv_tracks[ix]

    if init_in_reg:
        bv_tracks=get_tracks_init_in_region(bv_tracks)

    print(len(bv_tracks), 'true bvs with rain')

    return bv_tracks
