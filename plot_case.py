import warnings
import os.path as path
import numpy as np
import iris
from map_formatter import *
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import pdb
from read_kevin_bv_tracks import *
from read_surge import *
from format_clabel import *
from get_vorticity import *
from get_vortex_latlon_tilt import *
from get_case import *
from plot_case_qt import *
from plotting_functions import *

#------------------------------------------------------------------------------
# plot wind field uvw maps at 400, 600, 850 and 925 hPa
# then plot wind field uvw cross sections in NS and EW direction centred on the BV track centre
# inputs:
#    this_date - the date for which the data is given - to be added to the top of the figure
#    surge_str - the type of surge at this time as a string to be displayed along with the date
#    min_lon,max_lon, min_lat, max_lat -  the range for lon and lat coordinates to plot
#    u,v,w - data for pressure levels except 925, lats and lons
#    u925, v925, w925 - data jsut for 925hPa as this is read from a different file
#    below_mask indicates at which pressure levels cells become under ground
#    lons_ext, lats_ext - lons and lats of u and v
#    lons, lats - lons and lats of w
#    pressure, pix<xxx> - pressure levels of u,v,w and the indices at which to find differet pressure levels
#    track_lon850, track_lat850 - bv 850hPa location
#    track_lon925, track_lat925, best_match - bv 925hPa location - only vald of best_match==1
#    outdir - where to store the plot
#-----------------------------------------------------------------------------------
def plot_uvw(this_date, surge_str, min_lon, max_lon, min_lat, max_lat, u, v, w, u925, v925, w925, below_mask, lons_ext, lats_ext, lons, lats, pressure, pix950, pix900, pix850, pix600, pix400, pix500, track_lon850, track_lat850, track_lon925, track_lat925, best_match, outdir):
    
    # plot map of 400 hPa winds 
    # map of 600 hPa winds
    # map of 850 hPa winds showing BV track current centre
    # map of 925 hPa winds showing BV track current centre
    # longitude height profile of winds showing both track centres
    # latitude height profile of winds showing both track centres

    # u and v are on extended grid, w and 925 winds are not

    cb_top_pos=[0.98, 0.65, 0.02, 0.24]
    cb_mid_pos=[0.98, 0.38, 0.02, 0.24]
    cb_bottom_pos=[0.98, 0.1, 0.02, 0.24]
    nrows=3
    ncols=2

    u400=u[pix400,:,:]
    v400=v[pix400,:,:]
    w400=w[pix400,:,:]
    u600=u[pix600,:,:]
    v600=v[pix600,:,:]
    w600=w[pix600,:,:]
    u850=u[pix850,:,:]
    v850=v[pix850,:,:]
    w850=w[pix850,:,:]

    print('plotting uvw maps')

    fig = plt.figure(figsize = (9,13))
    title='(a) 400 hPa wind'
    con, sp=plot_uvw_map(fig, nrows,ncols, 1, title, min_lon, max_lon, min_lat, max_lat, u400, v400, w400, below_mask[pix400,:,:], lons_ext, lats_ext, lons, lats, track_lon850, track_lat850, '')
    plt.text(max_lon, max_lat*1.3, this_date.strftime('%d/%m/%Y %H:%M ')+surge_str, fontsize=12)

    title='(b) 600 hPa wind'
    con, sp=plot_uvw_map(fig, nrows,ncols, 2, title, min_lon, max_lon, min_lat, max_lat, u600, v600, w600, below_mask[pix600,:,:], lons_ext, lats_ext, lons, lats, track_lon850, track_lat850, '')
    cb = fig.add_axes(cb_top_pos)
    fig.colorbar(con, cax=cb,orientation = 'vertical')
    cb.set_ylabel(w_label, fontsize=10)

    title='(c) 850 hPa wind'
    con,sp=plot_uvw_map(fig, nrows,ncols, 3, title, min_lon, max_lon, min_lat, max_lat, u850, v850, w850, below_mask[pix850,:,:], lons_ext, lats_ext, lons, lats, track_lon850, track_lat850, col850)

    title='(d) 925 hPa wind'
    track_col=''
    if best_match>=0:
        track_col=col925 # put dot at current centre
    con, sp=plot_uvw_map(fig, nrows,ncols, 4, title, min_lon, max_lon, min_lat, max_lat, u925, v925, w925, below_mask[pix950,:,:], lons, lats, lons, lats, track_lon925, track_lat925, track_col)
    cb = fig.add_axes(cb_mid_pos)
    fig.colorbar(sp.lines, cax=cb,orientation = 'vertical')
    cb.set_ylabel(wspeed_label, fontsize=10)

    # plot wind cross sections
    ix_lat_ext=np.where(abs(lats_ext-track_lat850)<=1)
    ix_lat=np.where(abs(lats-track_lat850)<=1)
    this_v_cross_ew=np.mean(v[:,ix_lat_ext[0],:].data, axis=1)
    this_w_cross_ew=np.mean(w[:,ix_lat[0],:].data, axis=1)
    ix_lon_ext=np.where(abs(lons_ext-track_lon850)<=1)
    ix_lon=np.where(abs(lons-track_lon850)<=1)
    this_u_cross_ns=np.mean(u[:,:,ix_lon_ext[0]].data, axis=2)
    this_w_cross_ns=np.mean(w[:,:,ix_lon[0]].data, axis=2)
    vortex_lats, vortex_lons, this_vortex_height_ns, this_vortex_height_ew, this_tilt_north, this_tilt_east = find_vortex_height_and_tilt(track_lat850, track_lon850, this_u_cross_ns, this_v_cross_ew, lons_ext, lats_ext, pressure, pix850, pix500)

    track_lat_str='{l:2.2f}'.format(l=abs(track_lat850))
    if track_lat850>=0:
        track_lat_str=track_lat_str+'N'
    else:
        track_lat_str=track_lat_str+'S'
    title='(e) meridional wind at '+track_lat_str
    this_below_mask=np.nanmean(below_mask[:,ix_lat[0],:], axis=1)
    con_v=plot_uvw_ew_cross_section(fig, nrows,ncols, 5, title, min_lon, max_lon, this_v_cross_ew, this_w_cross_ew, this_below_mask, lons_ext, lons, pressure, track_lon850, track_lon925, best_match, vortex_lons, this_vortex_height_ew<1000)

    track_lon_str='{l:2.2f}E'.format(l=track_lon850)
    title='(f) zonal wind at '+track_lon_str
    this_below_mask=np.nanmean(below_mask[:,:,ix_lon[0]], axis=2)
    con_u=plot_uvw_ns_cross_section(fig, nrows,ncols, 6, title, min_lat, max_lat, this_u_cross_ns, this_w_cross_ns, this_below_mask, lats_ext, lats, pressure, track_lat850, track_lat925, best_match, vortex_lats, this_vortex_height_ns<1000)

    cb = fig.add_axes(cb_bottom_pos)
    fig.colorbar(con_v, cax=cb,orientation = 'vertical')
    cb.set_ylabel(uv_label, fontsize=10)

    outfile=outdir+this_date.strftime('uvw_%Y%m%d_%H.%M.png')
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()

    return this_vortex_height_ns, this_vortex_height_ew, this_tilt_north, this_tilt_east

#------------------------------------------------------------------------------
# plot relative vorticity at 850 hPa and 925 hPa with horizontal winds as streamlines,
# plot rain with horizontal winds as streamlines
# plot PV at 850 hPa
# plot PV cross sections in NS and EW direction centred on the BV track centre
# inputs:
#    this_date - the date for which the data is given - to be added to the top of the figure
#    surge_str - the type of surge at this time as a string to be displayed along with the date
#    min_lon,max_lon, min_lat, max_lat -  the range for lon and lat coordinates to plot
#    u850,v850 - [nlat,nlon] horizontal wind at 850 hPa
#    u925, v925 - [nlat,nlon] horizontal wind at 925hPa
#    rain - [nlat,nlon]
#    vort850, vort925 - [nlat,nlon] relative vorticity at 850 and 925
#    pv - [nlat,nlon] potential vorticity
#    below_mask - [nlat,nlon] indicates at which pressure levels cells are under ground
#    lons_ext, lats_ext - lons and lats of u's and v's
#    lons, lats - lons and lats of 
#    lons_p, lats_p - lons and lats of rain
#    pressure, pix<xxx> - pressure levels of u,v,w and the indices at which to find differet pressure levels
#    track_lon850, track_lat850 - bv 850hPa location
#    track_lon925, track_lat925, best_match - bv 925hPa location - only vald of best_match==1
#    outdir - where to store the plot
#-----------------------------------------------------------------------------------
def plot_rain_vort_pv(this_date, surge_str, min_lon,max_lon,min_lat,max_lat, u850, v850, u925, v925, rain, vort850, vort925, pv, below_mask, lons_ext, lats_ext, lons, lats, pressure, pix850, pix900, pix950, lons_p, lats_p, track_lon850, track_lat850, track_lon925, track_lat925, best_match, outdir):

    #-----------------------------------------------------------
    # plot 850 hPa vorticity
    # 925 hPa vorticity
    # imerg rain
    # PV at 850 hPa
    # longitude height profile of PV with both track centres
    # latitude height profile of PV with both track centres
    #-----------------------------------------------------------
    cb_top1_pos=[0.92, 0.65, 0.015, 0.24]
    cb_top2_pos=[0.999, 0.65, 0.015, 0.24]
    cb_mid_pos=[0.98, 0.38, 0.02, 0.24]
    cb_bottom_pos=[0.98, 0.1, 0.02, 0.24]
    nrows=3
    ncols=2

    fig = plt.figure(figsize = (9,13))
    title='(a) 850 hPa \u03B6 and streamlines'
    this_below_mask=below_mask[pix850,:,:]
    con_v, sp=plot_vort_map(fig, nrows,ncols, 1, title, min_lon,max_lon,min_lat,max_lat, u850, v850, vort850, this_below_mask, lons_ext, lats_ext, track_lon850, track_lat850, col850)
    plt.text(max_lon, max_lat*1.3, this_date.strftime('%d/%m/%Y %H:%M ')+surge_str, fontsize=12)
   
    title='(b) 925 hPa \u03B6 and streamlines'
    this_below_mask=below_mask[pix950,:,:]
    track_col=''
    if best_match>=0:
        track_col=col925
    con_v, sp=plot_vort_map(fig, nrows,ncols, 2, title, min_lon,max_lon,min_lat,max_lat, u925, v925, vort925, this_below_mask, lons, lats, track_lon925, track_lat925, col925)
    cb = fig.add_axes(cb_top1_pos)
    fig.colorbar(con_v, cax=cb,orientation = 'vertical')
    cb.set_ylabel(vort_label, fontsize=8)
    cb = fig.add_axes(cb_top2_pos)
    fig.colorbar(sp.lines, cax=cb,orientation = 'vertical')
    cb.set_ylabel(wspeed_label, fontsize=8)

    title='(c) IMERG and 850 hPa streamlines'
    this_below_mask=below_mask[pix850,:,:]
    con_r, sp=plot_rain_map(fig, nrows,ncols, 3, title, min_lon,max_lon,min_lat,max_lat, u850, v850, rain, this_below_mask, lons_ext, lats_ext, lons_p, lats_p, track_lon850, track_lat850, col850)
    
    cb = fig.add_axes(cb_mid_pos)
    fig.colorbar(con_r, cax=cb,orientation = 'vertical')
    cb.set_ylabel(rainrate_label, fontsize=8)
        
    pv850=pv[pix850,:,:]
    title='(d) 850 hPa PV and streamlines'
    plot_pv_map(fig, nrows,ncols, 4, title, min_lon,max_lon,min_lat,max_lat, u850, v850, pv850, this_below_mask, lons_ext, lats_ext, lons, lats, track_lon850, track_lat850, col850)

    # plot potential vorticity cross sections
    track_lat_str='{l:2.2f}'.format(l=abs(track_lat850))
    if track_lat850>=0:
        track_lat_str=track_lat_str+'N'
    else:
        track_lat_str=track_lat_str+'S'
    title='(e) PV at '+track_lat_str
    ix_lat=np.where(abs(lats-track_lat850)<=1)
    this_pv=np.mean(pv[:,ix_lat[0],:], axis=1)
    this_below_mask=np.nanmean(below_mask[:,ix_lat[0],:], axis=1)
    con_pv=plot_pv_ew_cross_section(fig, nrows,ncols, 5, title, min_lon,max_lon, this_pv, this_below_mask, lons, pressure, track_lon850, track_lon925, best_match)

    track_lon_str='{l:2.2f}E'.format(l=track_lon850)
    title='(f) PV at '+track_lon_str
    ix_lon=np.where(abs(lons-track_lon850)<=1)
    this_pv=np.mean(pv[:,:,ix_lon[0]], axis=2)
    this_below_mask=np.nanmean(below_mask[:,:,ix_lon[0]], axis=2)
    con_pv=plot_pv_ns_cross_section(fig, nrows,ncols, 6, title, min_lat,max_lat, this_pv, this_below_mask, lats, pressure, track_lat850, track_lat925, best_match)
    cb = fig.add_axes(cb_bottom_pos)
    fig.colorbar(con_pv, cax=cb,orientation = 'vertical')
    cb.set_ylabel(pv_label, fontsize=8)

    outfile=outdir+this_date.strftime('imerg_pv_%Y%m%d_%H.%M.png')
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()


#------------------------------------------------------------------------------
# Given the data for a specific BV case do the following:
# 1. work out potential temperature
# 2. work out potential vorticity
# 3. work out vorticity at 850 and 925 hPa and add the track centred mean 850 vorticity to the BV track object
# 4. Then for each time plot uvw maps and cross sections, plot rain vorticity and pv maps and cross sections
#    calculate vortex height and tilt
# 5. finally print the timeseries of mean rain and vort, vortex height and tilt
#------------------------------------------------------------------------------
def plot_case(track_850, track_lats925, track_lons925, best_match, min_lat, max_lat, min_lon, max_lon, surge_indices, surge_times, era5_u, era5_v, era5_w, era5_u925, era5_v925, era5_w925, era5_q, era5_T, era5_below_mask, rain, outdir):

    era5_theta=get_potential_temperature(era5_T)
    era5_pv=get_potential_vorticity(era5_u, era5_v, era5_theta)
    # every 6hours for requested case plot

    start_date=track_850.track_times[0]
    end_date=track_850.track_times[-1]
    track_lons850=track_850.vort_data[0,:]
    track_lats850=track_850.vort_data[1,:]
    
    lons_ext = era5_u.coord('longitude').points
    lats_ext = era5_u.coord('latitude').points
    lons = era5_w.coord('longitude').points
    lats = era5_w.coord('latitude').points
    pressure = era5_u.coord('pressure_level').points
    # use these pressure levels to get below_mask maps
    pix900=np.where(pressure==900)[0][0]
    pix950=np.where(pressure==950)[0][0]
    pix850=np.where(pressure==850)[0][0]
    # use these pressure levels in wind plots and height and tilt measurement
    pix600=np.where(pressure==600)[0][0]
    pix400=np.where(pressure==400)[0][0]
    pix500=np.where(pressure==500)[0][0]
    wind_t_coord=era5_u.coord('time')
    wind_dates = wind_t_coord.units.num2date(wind_t_coord.points)
    
    u850=era5_u[:,pix850,:,:]
    v850=era5_v[:,pix850,:,:]
    print('getting vorticity')
    vort850,lons_vort,lats_vort=get_vorticity(u850,v850)
    vort925,lons_vort925,lats_vort925=get_vorticity(era5_u925,era5_v925)
    # get storm centred average vorticity at 850 hPa
    track_850.add_vorticity(vort850*vort_factor, lons_ext, lats_ext, wind_dates)
    
    rain_t_coord=rain.coord('time')
    rain_dates = rain_t_coord.units.num2date(rain_t_coord.points)
    lons_p = rain.coord('longitude').points
    lats_p = rain.coord('latitude').points

    this_date=start_date

    surge_years=np.asarray([surge_time.year for surge_time in surge_times])        
    surge_months=np.asarray([surge_time.month for surge_time in surge_times])
    surge_days=np.asarray([surge_time.day for surge_time in surge_times])
    for t in range(track_850.nt):
        this_date=track_850.track_times[t]
        print(this_date)
        
        six=np.where((surge_years==this_date.year) & (surge_months==this_date.month) & (surge_days==this_date.day))
        if len(six[0])==0:
            print('no surge',  this_date)
        this_surge=surge_indices[six[0][0]]
        surge_str=' surge='+surge_text[this_surge]
        
        ix_wind=np.where(wind_dates==track_850.track_times[t])
        if len(ix_wind[0])==0:
            print('no wind')
            pdb.set_trace()
        u=era5_u[ix_wind[0][0]]
        v=era5_v[ix_wind[0][0]]
        w=era5_w[ix_wind[0][0]]
        u925=era5_u925[ix_wind[0][0]]
        v925=era5_v925[ix_wind[0][0]]
        w925=era5_w925[ix_wind[0][0]]
        pv=era5_pv[ix_wind[0][0]]
        below_mask=era5_below_mask[ix_wind[0][0]]
        
        this_vortex_height_ns, this_vortex_height_ew, this_tilt_north, this_tilt_east=plot_uvw(this_date, surge_str, min_lon, max_lon, min_lat, max_lat, u, v, w, u925, v925, w925, below_mask, lons_ext, lats_ext, lons, lats, pressure, pix950, pix900, pix850, pix600, pix400, pix500, track_lons850[t], track_lats850[t], track_lons925[t], track_lats925[t], best_match[t], outdir)
        if t==0:
            vortex_height_ns=[this_vortex_height_ns]
            vortex_height_ew=[this_vortex_height_ew]
            vortex_height=[np.amax([this_vortex_height_ns,this_vortex_height_ew])]
            vortex_tilt_east=[this_tilt_east]
            vortex_tilt_north=[this_tilt_north]
        else:
            vortex_height_ns.append(this_vortex_height_ns)
            vortex_height_ew.append(this_vortex_height_ew)
            vortex_height.append(np.amax([this_vortex_height_ns,this_vortex_height_ew]))
            vortex_tilt_east.append(this_tilt_east)
            vortex_tilt_north.append(this_tilt_north)
        
        ix_rain=np.where(rain_dates+dt.timedelta(minutes=15)==track_850.track_times[t])
        if len(ix_rain[0])==0:
            print('no rain')
            pdb.set_trace()
        this_rain=rain[ix_rain[0][0],:,:]
        u850=u[pix850,:,:]
        v850=v[pix850,:,:]
        plot_rain_vort_pv(this_date, surge_str, min_lon,max_lon,min_lat,max_lat, u850, v850, u925, v925, this_rain, vort850[ix_wind[0][0]], vort925[ix_wind[0][0]],pv, below_mask, lons_ext, lats_ext, lons, lats, pressure, pix850, pix900, pix950, lons_p, lats_p, track_lons850[t], track_lats850[t], track_lons925[t], track_lats925[t], best_match[t], outdir)

    #   plot accumulated rainfall map
    acc_rain_levels=np.arange(10)*20+20
    total_rain = rain.collapsed('time', iris.analysis.SUM)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    print('max total rain=',np.amax(total_rain.data))
    con = ax.contourf(lons_p, lats_p, total_rain.data,levels=acc_rain_levels,extend='max')
    cb = plt.colorbar(con, orientation = 'vertical')
    cb.set_label('accumulated rain (mm)')
    ax.plot(track_lons850,track_lats850, color=col850, zorder=2, label='850 hPa')
    track_hours=np.asarray([time.hour for time in track_850.track_times])
    ix=np.where(track_hours==0)
    ax.scatter(track_lons850[ix],track_lats850[ix], marker='o',s=8,color=col850, zorder=2)
    tix=np.where(best_match>=0)
    ax.plot(track_lons925[tix],track_lats925[tix], color=col925, zorder=2, label='925 hPa')
    ax.legend()
    map_formatter(ax,tick_base_x=tick_base_x,tick_base_y=tick_base_y,labelsize=8,res='50m')
    ax.set_ylim(min_lat,max_lat)
    ax.set_xlim(min_lon,max_lon)
    outfile=outdir+'rain_accumulation.png'
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()

    # plot whole timeseries data for requested case:
    fig = plt.figure()
    nrows=1
    ncols=1
    hours, hours850, date_labels=plot_timeseries_rain_vort(fig,nrows,ncols,1,track_850, surge_years, surge_months, surge_days, surge_indices,'')
    outfile=outdir+'rain_vort_timeseries.png'
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()
    
    fig = plt.figure()
    hours, hours850, date_labels=plot_timeseries_height_tilt(fig,nrows,ncols,1,track_850, surge_years, surge_months, surge_days, surge_indices, vortex_height, vortex_tilt_north, vortex_tilt_east,'')

    outfile=outdir+'height_tilt_timeseries.png'
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(hours850, vortex_height_ns, linestyle='dotted',color='k', label='NS')
    ax.plot(hours850, vortex_height_ew, linestyle='dashed',color='k', label='EW')
    ax.invert_yaxis()
    ax.legend()    
    ylim=ax.get_ylim()
    ndays=int(len(hours850)/4)
    offset_hour=2 # this is hour=12 on 6 hourly data
    for d in range(ndays):
        this_surge=surge_indices[d+first_surge_tix]
        ax.text(hours[d*4]+offset_hour, ylim[0]+0.01*(ylim[1]-ylim[0]), surge_text[this_surge], color='red', horizontalalignment='center')
    ax.set_xlabel('date')
    # plot labels only on the hour=0 dates
    ax.set_xticks(hours[::4])
    ax.set_xticklabels(date_labels, rotation=45)
    ax.set_ylabel('vortex height (hPa)')
    ax.tick_params(axis='y', colors='k')
    outfile=outdir+'vortex_height_ns_ew.png'
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()

def plot_composites(composites_ew, composites_ns, delta_lons, delta_lats, pressure, outdir):
    u_ix=0
    v_ix=1
    w_ix=2
    q_ix=3
    thetae_ix=4
    nrows=2
    ncols=2
    fig = plt.figure(figsize=(8,8))
    plot_composite_cross_sections(fig, nrows,ncols, 1, 'EW', composites_ew[:,:,v_ix], composites_ew[:,:,w_ix], composites_ew[:,:,q_ix], composites_ew[:,:,thetae_ix], delta_lons, pressure, 'longitude from centre')
    plot_composite_cross_sections(fig, nrows,ncols, 3, 'NS', composites_ns[:,:,u_ix], composites_ns[:,:,w_ix], composites_ns[:,:,q_ix], composites_ns[:,:,thetae_ix], delta_lats, pressure, 'latitude from centre')
    outfile=outdir+'composites.png'
    plt.savefig(outfile,dpi=100,bbox_inches='tight')
    plt.close()
