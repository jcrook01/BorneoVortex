import warnings
import numpy as np
import pdb

# to determine the height of a vortex we will look for the location where winds within 1 degree either side of
# a centre point are of opposite sign. This could be a little bit away from the max vorticity centre -
# allow up to 3 degrees away.
# The centre of rotation may shift as we go up due to tilt.
# We will start at lower levels using the max vorticity centre from the tracking and work upwards allowing
# up to 3 degrees shift at each pressure level.
# If we cannot find a centre within 3 degrees of given centre where winds are of opposite sign then we assume
# we have reached the top of the vortex and don't try any further up.
# The highest pressure level where we have found a centre gives us the vortex height and the shift in the centre
# with height gives us the tilt.
# This is done using the cross sections south-north and west-east using meridional and zonal winds. 
# So we get 2 heights and 2 tilts

ddeg=3  # look for changes in wind sign over this degrees either side of centre - allows for tilt but still don't want to be too far away

def get_vortex_centre_lat(cenlat,lats,u_ns):
    vortex_ns_found=False
    vortex_lat=np.nan
    # look at u +/- ddeg dgrees from current centre and find where we have positive to the south and negative to the north
    ix_centre_lat=np.where((lats>=cenlat-ddeg) & (lats<=cenlat+ddeg))[0]
    centre_u=u_ns[ix_centre_lat]
    centre_lats=lats[ix_centre_lat]
    ix_neg=np.where(centre_u<0)
    ix_pos=np.where(centre_u>0)
    if len(ix_pos[0])>0 and len(ix_neg[0])>0:
        # there are some positive and negative values
        pos_val=np.zeros(len(centre_u), int)
        pos_val[ix_pos]=1
        # as lats go from big to small we need to look for negative to positive change
        pos_diff=pos_val[1:]-pos_val[:-1]
        ix_change=np.where(pos_diff==1)
        change_lats=centre_lats[ix_change[0]]
        if len(change_lats)>0:
            # find the change closest to the original centre
            ix_sort=np.argsort(abs(change_lats-cenlat))
            change_lats=change_lats[ix_sort]
            vortex_lat=change_lats[0]
            vortex_ns_found=True
            #print('NS vortex found: centre at ', vortex_lat)
        #else:
            #print('NS vortex not found: no change negative to positive')
    #else:
        #print('NS vortex not found: no positive and negative u near centre, len pos, len neg=', len(ix_pos[0]), len(ix_neg[0]))
            
    return vortex_lat, vortex_ns_found
    
def get_vortex_centre_lon(cenlon,lons,v_ew):
    vortex_ew_found=False
    vortex_lon=np.nan
    # look at v +/- ddeg dgrees from current centre and find where we have positive to the east and negative to the west
    ix_centre_lon=np.where((lons>=cenlon-ddeg) & (lons<=cenlon+ddeg))[0]
    centre_lons=lons[ix_centre_lon]
    centre_v=v_ew[ix_centre_lon]
    ix_neg=np.where(centre_v<0)
    ix_pos=np.where(centre_v>0)
    if len(ix_pos[0])>0 and len(ix_neg[0])>0:
        # there are some positive and negative values
        pos_val=np.zeros(len(centre_v), int)
        pos_val[ix_pos]=1
        pos_diff=pos_val[1:]-pos_val[:-1]
        # lons increase so we need to look for negative to positive change
        ix_change=np.where(pos_diff==1)
        change_lons=centre_lons[ix_change[0]]
        if len(change_lons)>0:
            # find the change closest to the original centre
            ix_sort=np.argsort(abs(change_lons-cenlon))
            change_lons=change_lons[ix_sort]
            vortex_lon=change_lons[0]
            vortex_ew_found=True
            #print('EW vortex found: centre at ', vortex_lon)
        #else:
            #print('EW vortex not found: no change negative to positive')
    #else:
        #print('EW vortex not found: no positive and negative v near centre, len pos, len neg=', len(ix_pos[0]), len(ix_neg[0]))

    return vortex_lon,vortex_ew_found

def get_vortex_tilt(vortex_centre_lats, vortex_centre_lons, pressures, pix850, pix500):
    tilt_north=0
    tilt_east=0

    # use 850 to 500 hPa centre to determine tilt as centre above 500 hPa is not always too clear
    if np.isnan(vortex_centre_lats[pix850])==False:
      if np.isnan(vortex_centre_lats[pix500])==False:
        tilt_north=(vortex_centre_lats[pix850]-vortex_centre_lats[pix500])/(500-850)
      else:
        # vortex did not get as high as 500 so use a lower value if there is one
        ix_top=np.where(np.isnan(vortex_centre_lats)==False)
        if len(ix_top[0])>0:
            poss_pressures=pressures[ix_top[0]]
            pix=np.argmin(poss_pressures)
            tilt_north=(vortex_centre_lats[pix850]-vortex_centre_lats[pix])/(pressures[pix]-850)
        
    if np.isnan(vortex_centre_lons[pix850])==False:
      if np.isnan(vortex_centre_lons[pix500])==False:
        tilt_east=(vortex_centre_lons[pix850]-vortex_centre_lons[pix500])/(500-850)
      else:
        # vortex did not get as high as 500 so use a lower value if there is one
        ix_top=np.where(np.isnan(vortex_centre_lons)==False)
        if len(ix_top[0])>0:
            poss_pressures=pressures[ix_top[0]]
            pix=np.argmin(poss_pressures)
            tilt_east=(vortex_centre_lons[pix850]-vortex_centre_lons[pix])/(pressures[pix]-850)
    print('tilt north and east', tilt_north, tilt_east)
    return tilt_north, tilt_east

def find_vortex_height_and_tilt(track_lat850, track_lon850, this_u_cross_ns, this_v_cross_ew, lons, lats, pressure, pix850, pix500):
    # find the vortex centres up to 200 hPa
    ix=np.where(pressure>=200)
    press_upto200=pressure[ix[0]]
    nplev=len(ix[0])
    vortex_lats=np.zeros(nplev)+np.nan
    vortex_lons=np.zeros(nplev)+np.nan
    this_vortex_height_ns=1000
    this_vortex_height_ew=1000
    this_tilt_north=0
    this_tilt_east=0
    #print('calculating height and tilt')
    cenlat=track_lat850
    cenlon=track_lon850
    ix=np.where((pressure>=200) & (pressure<850))
    press_above850=np.sort(pressure[ix[0]])[::-1] # go from 850 to 200
    ix=np.where(pressure>=850)
    press_below850=np.sort(pressure[ix[0]]) # go from 850 to 1000
    for p in press_below850:
        pix=np.where(pressure==p)
        pixto200=np.where(press_upto200==p)
        #print('pressure', pressure[pix], 'hPa: orig centre lon lat', cenlon,cenlat)
        # find the centre at each pressure level 
        this_vortex_lat, vortex_ns_found=get_vortex_centre_lat(cenlat,lats,this_u_cross_ns[pix[0][0],:])
        if vortex_ns_found:
            vortex_lats[pixto200]=this_vortex_lat

        this_vortex_lon, vortex_ew_found=get_vortex_centre_lon(cenlon,lons,this_v_cross_ew[pix[0][0],:])
        if vortex_ew_found:
            vortex_lons[pixto200]=this_vortex_lon

    cenlat=track_lat850
    cenlon=track_lon850
    vortex_ns_found=True
    vortex_ew_found=True
    for p in press_above850:
        pix=np.where(pressure==p)[0][0]
        pixto200=np.where(press_upto200==p)
        #print('pressure', pressure[pix], 'hPa: orig centre lon lat', cenlon,cenlat)
        # find the centre at each pressure level until we cant find it
        if pix==pix850 or vortex_ns_found:
            this_vortex_lat, vortex_ns_found=get_vortex_centre_lat(cenlat,lats,this_u_cross_ns[pix,:])
            if vortex_ns_found:
                vortex_lats[pixto200]=this_vortex_lat
                cenlat=this_vortex_lat
                this_vortex_height_ns=p

        if pix==pix850 or vortex_ew_found:
            this_vortex_lon, vortex_ew_found=get_vortex_centre_lon(cenlon,lons,this_v_cross_ew[pix,:])
            if vortex_ew_found:
                vortex_lons[pixto200]=this_vortex_lon
                cenlon=this_vortex_lon
                this_vortex_height_ew=p
                
        if vortex_ns_found==False and vortex_ew_found==False:
            break # neither centre found so stop loop

    if this_vortex_height_ns<850 and this_vortex_height_ew<850:        
        this_tilt_north, this_tilt_east=get_vortex_tilt(vortex_lats, vortex_lons, pressure, pix850, pix500) 
        
    print('this vortex height', this_vortex_height_ns, this_vortex_height_ew)

    return vortex_lats, vortex_lons, this_vortex_height_ns, this_vortex_height_ew, this_tilt_north, this_tilt_east