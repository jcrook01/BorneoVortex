import warnings
import os.path as path
import numpy as np
import iris
from map_formatter import *
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from read_kevin_bv_tracks import *
from read_surge import *
import pdb

wind_levels=[-20,-15,-10,-5,-1,1,5,10,15,20]
w_levels=np.arange(10)*0.4-1.8
wprof_levels=np.arange(10)*0.4-1.8
w_linestyles=['dotted','dotted','dotted','dotted','dotted','solid','solid','solid','solid','solid']
rain_levels=np.arange(10)*2+2
vort_factor=1.0e5
vort_levels=np.arange(12)*2-9
pv_levels=np.arange(10)*0.2-0.9
q_levels=np.arange(14)*0.2-1.1
thetae_levels=np.arange(14)-6.5

w_label='vertical wind (Pa s$^{-1}$)'
uv_label='horizontal wind (m s$^{-1}$)'
wspeed_label='horizontal windspeed (m s$^{-1}$)'
rainrate_label='rain (mm h$^{-1}$)'
vort_label='\u03B6$ (10$^{-5}$ s$^{-1}$)'
Ra_label='R$_{a}$ (mm h$^{-1}$)'  # rain associated with BV (area averaged)
centre_vort_label='\u03B6$_{c}$ (10$^{-5}$ s$^{-1}$)' # the BV tracked vorticity
mean_vort_label='\u03B6$_{a}$ (10$^{-5}$ s$^{-1}$)'  # the BV area averaged vorticity
pv_label='PV (PVU)'
tick_base_x=5
tick_base_y=5
# constants for streamlines
str_density=2
arr_sz=0.5
linewidth=1

col925='cyan'
col850='magenta'

wind_cmap = plt.get_cmap('bwr')
pv_cmap = plt.get_cmap('bwr')
q_cmap = plt.get_cmap('viridis')
delta_q_cmap = plt.get_cmap('bwr_r')
mask_cmap=plt.get_cmap('binary',2)
min_wspeed=0
max_wspeed=12
stream_cmap=plt.get_cmap('gray_r',max_wspeed)
stream_norm=plt.Normalize(min_wspeed, max_wspeed)




#-----------------------------------------------------------------------------------------------------
# plot a vertical cross section of data_cross shaded 
# also plots the position of the BV and how it tilts if required
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_x,max_x - defines the x coordinate extent to plot 
#     data_cross - data to plot with dimensions[np,nx]
#     x - the x coordinates of the data
#     pressure - the y coordinate
#     shade_levels - the contour levels in the shaded field
#     shade_cmap - the cmap to use for filling contours for shaded field
#     x_label - the label for the x axis (longitude or latitude)
#     track_x850 - this is the x location of the BV in 850 hPa tracking to plot - if INVALID_LAT_LON dont plot
#     track_x925 - this is the x location of the BV in 925 hPa tracking to plot - if INVALID_LAT_LON dont plot
#     bv_x - the x location of the BV as estimated at each pressure level up to 200 hPa
#     plot_bv - whether to plot the bv_x
# returns:
#     ax - the axes object
#     con - the contour plot object so colour bar can be added to the plot
#-----------------------------------------------------------------------------------------------------
def plot_cross_section(fig, nrows,ncols, nplot, title, min_x, max_x, data_cross, x, pressure, shade_levels, shade_cmap, x_label, track_x850, track_x925, bv_x, plot_bv):

    ax = fig.add_subplot(nrows,ncols,nplot)
    plt.title(title, loc='left')
    con=ax.contourf(x_shade, pressure, shade_cross, levels=shade_levels, extend='both',cmap=shade_cmap)
    ax.set_xlim(min_x,max_x)
    ax.set_ylim(1000, 200)
    ax.set_xlabel(x_label)
    ax.set_ylabel('pressure (hPa)')
    if track_x850!=INVALID_LAT_LON:
        ax.scatter(track_x850, 850, color=col850,marker='*', zorder=4) # put marker at current 850 BV centre
    if track_x925!=INVALID_LAT_LON:
        ax.scatter(track_x925, 925, color=col925,marker='*', zorder=4) # put marker at current 925 BV centre
    if plot_bv==1:
        ix=np.where(pressure>=200)
        bv_pressures=pressure[ix[0]]
        ax.plot(bv_x, bv_pressures, color='k')
    return ax,con

#-------------------------------------------------------------------------------------------------------------------
# plot a vertical cross section of w as contour lines and horizontal wind as shading and another plot with anomaly
# in q_cross as shading and anomaly in thetea_cross as contour lines
# inputs:
#     fig - the matplotlib figure to add these subplots to
#     nrows, ncols, first_nplot - first_nplot is which plot out of the nrows x ncols for the 1st plot
#     title - base title of this plot
#     horiz_wind_cross, w_cross, q_cross, thetae_cross - data to plot with dimensions [np,nx]
#     delta_x - the x coordinates of the data for x axis
#     pressure - the pressure coordinates of the data for y axis
#     xlabel - label to display on x axis
#-------------------------------------------------------------------------------------------------------------------
def plot_composite_cross_sections(fig, nrows,ncols, first_nplot, title, horiz_wind_cross, w_cross, q_cross, thetae_cross, delta_x, pressure, xlabel):
 
    ax,con_horiz=plot_cross_section(fig, nrows,ncols, first_nplot+1, title+' horizontal winds and omega', delta_x[0], delta_x[-1], horiz_wind_cross, delta_x, pressure, wind_levels, wind_cmap, xlabel, INVALID_LAT_LON, INVALID_LAT_LON, [], 0)
    cb=plt.colorbar(con_horiz, orientation = 'horizontal')
    cb.set_label(uv_label, fontsize=10)
    con_w=ax.contour(delta_x, pressure, w_cross, levels=wprof_levels,linestyles=w_linestyles, extend='both')
    ax.clabel(con_w, con_w.levels, inline=True, fontsize=10)
    ax.axvline(x=0)
    
    mean_q=np.mean(q_cross, axis=1)
    delta_q_cross=np.asarray([q_cross[p,:]-mean_q[p] for p in range(len(pressure))])
    mean_thetae=np.mean(thetae_cross, axis=1)
    delta_thetae_cross=np.asarray([thetae_cross[p,:]-mean_thetae[p] for p in range(len(pressure))])
    
    ax,con_q=plot_cross_section(fig, nrows,ncols, first_nplot+1, title+' delta q and theta_e', delta_x[0], delta_x[-1], delta_q_cross, delta_x, pressure, q_levels, delta_q_cmap, xlabel, INVALID_LAT_LON, INVALID_LAT_LON, [], 0)

    cb=plt.colorbar(con_q, orientation = 'horizontal')
    cb.set_label('specific humidity (g kg$^{-1}$)', fontsize=10)
    con_c=ax.contour(lons_w, pressure, delta_thetae_cross, levels=thetae_levels, extend='both')
    ax.clabel(con_th, con_th.levels, inline=True, fontsize=10)
    ax.axvline(x=0)


#-----------------------------------------------------------------------------------------------------
# plot an EW vertical cross section of horizontal wind v shaded and w as contour lines 
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lon,max_lon - defines the x coordinate extent to plot 
#     this_v_cross,this_w_cross, below_mask - data[np,nlon] already selected at the right latitude and time
#            below_mask indicates when a cell is below the surface so we can mask out data
#     lons_v,lons_w - the coordinates of the data - v may have an extended grid compared to w
#     pressure - the y coordinates of the data
#     track_lon850 - this is the longitude of the BV in 850 hPa tracking to plot 
#     track_lon925 - this is the longitude of the BV in 925 hPa tracking to plot 
#     best_match - if 1 then track_lon925 is valid (there was a matching 925 BV at this time)
#     bv_lons - the longitude of the BV as estimated at each pressure level up to 200 hPa
#     plot_bv_lons - whether to plot the bv_lons
# returns:
#     con_v - the contour plot object for this_v_cross_ew so colour bar can be added to the plot
#-----------------------------------------------------------------------------------------------------
def plot_uvw_ew_cross_section(fig, nrows,ncols, nplot, title, min_lon, max_lon, this_v_cross_ew, this_w_cross_ew, below_mask, lons_v, lons_w, pressure, track_lon850, track_lon925, best_match, bv_lons, plot_bv_lons):
    # plot ew wind cross section
    track_x925=INVALID_LAT_LON
    if best_match:
       track_x925=track_lon925
    ax,con_v=plot_cross_section(fig, nrows,ncols, nplot, title, min_lon, max_lon, this_v_cross_ew, lons_v, pressure, wind_levels, wind_cmap, 'longitude', track_lon850, track_x925, bv_lons, plot_bv_lons)
    con_c=ax.contour(lons_w, pressure, this_w_cross_ew, levels=wprof_levels,linestyles=w_linestyles, extend='both')
    con_m=ax.contourf(lons_w, pressure, below_mask, levels=[-1,0,1],cmap=mask_cmap)


    return con_v

#-----------------------------------------------------------------------------------------------------
# plot a NS vertical cross section of horizontal wind v shaded and w as contour lines 
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lat,max_lat - defines the x coordinate extent to plot 
#     this_u_cross,this_w_cross, below_mask - data[np,nlat] already selected at the right longitude and time
#            below_mask indicates when a cell is below the surface so we can mask out data
#     lats_u,lats_w - the x coordinates of the data - u may have an extended grid compared to w
#     pressure - the y coordinates of the data
#     track_lat850 - this is the latitude of the BV in 850 hPa tracking to plot 
#     track_lat925 - this is the latitude of the BV in 925 hPa tracking to plot 
#     best_match - if 1 then track_lat925 is valid (there was a matching 925 BV at this time)
#     bv_lats - the latitude of the BV as estimated at each pressure level up to 200 hPa
#     plot_bv_lats - whether to plot the bv_lats
# returns:
#     con_u - the contour plot object for this_u_cross_ns so colour bar can be added to the plot
#-----------------------------------------------------------------------------------------------------
def plot_uvw_ns_cross_section(fig, nrows,ncols, nplot, title, min_lat, max_lat, this_u_cross_ns, this_w_cross_ns, below_mask, lats_uv, lats_w, pressure, track_lat850, track_lat925, best_match, bv_lats, plot_bv_lats):
 
    track_x925=INVALID_LAT_LON
    if best_match:
       track_x925=track_lat925
    ax,con_u=plot_cross_section(fig, nrows,ncols, nplot, title, min_lon, max_lon, this_u_cross_ns, lats_uv, pressure, wind_levels, wind_cmap, 'latitude', track_lat850, track_x925, bv_lats, plot_bv_lats)
    con_c=ax.contour(lats_w, pressure, this_w_cross_ns, levels=wprof_levels,linestyles=w_linestyles, extend='both')
    con_m=ax.contourf(lats_w, pressure, below_mask, levels=[-1,0,1],cmap=mask_cmap)
 
    return con_u

#-----------------------------------------------------------------------------------------------------
# plot an EW vertical cross section of pv shaded and showing track position 
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lon,max_lon - defines the x coordinate extent to plot 
#     pv_cross_ew - data[np,nlon] already selected at the right latitude and time
#     below_mask - indicates when a cell is below the surface so we can mask out data
#     lons - the x coordinates of the data
#     pressure - the y coordinates of the data
#     track_lon850 - this is the longitude of the BV in 850 hPa tracking to plot 
#     track_lon925 - this is the longitude of the BV in 925 hPa tracking to plot 
#     best_match - if 1 then track_lon925 is valid (there was a matching 925 BV at this time)
# returns:
#     con_pv - the contour plot object for pv_cross_ew so colour bar can be added to the plot
#-----------------------------------------------------------------------------------------------------
def plot_pv_ew_cross_section(fig, nrows,ncols, nplot, title, min_lon,max_lon, pv_cross_ew, below_mask, lons, pressure, track_lon850, track_lon925, best_match):

    track_x925=INVALID_LAT_LON
    if best_match:
       track_x925=track_lon925
       
    ax,con_pv=plot_cross_section(fig, nrows,ncols, nplot, title, min_lon, max_lon, pv_cross_ew, lons, pressure, pv_levels, pv_cmap, 'longitude', track_lon850, track_x925, [], 0)
    con_m=ax.contourf(lons, pressure, below_mask, levels=[-1,0,1],cmap=mask_cmap)

    return con_pv

#-----------------------------------------------------------------------------------------------------
# plot a NS vertical cross section of pv shaded and showing track position
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lat,max_lat - defines the x coordinate extent to plot 
#     pv_cross_ns - data[np,nlat] already selected at the right latitude and time
#     below_mask - indicates when a cell is below the surface so we can mask out data
#     lats - the x coordinates of the data
#     pressure - the y coordinates of the data
#     track_lat850 - this is the latitude of the BV in 850 hPa tracking to plot 
#     track_lat925 - this is the latitude of the BV in 925 hPa tracking to plot 
#     best_match - if 1 then track_lon925 is valid (there was a matching 925 BV at this time)
# returns:
#     con_pv - the contour plot object for pv_cross_ew so colour bar can be added to the plot
#-----------------------------------------------------------------------------------------------------
def plot_pv_ns_cross_section(fig, nrows,ncols, nplot, title, min_lat, max_lat, pv_cross_ns, below_mask, lats, pressure, track_lat850, track_lat925, best_match):

    track_x925=INVALID_LAT_LON
    if best_match:
       track_x925=track_lat925

    ax,con_pv=plot_cross_section(fig, nrows,ncols, nplot, title, min_lat, max_lat, pv_cross_ns, lats, pressure, pv_levels, pv_cmap, 'latitude', track_lat850, track_x925, bv_lons, plot_bv_lons)
    con_m=ax.contourf(lats, pressure, below_mask, levels=[-1,0,1],cmap=mask_cmap)
    return con_pv


#-------------------------------------------------------------------------------------------------------------------
# plot a map of w as contours and u and v as streamlines
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lon,max_lon,min_lat,max_lat - defines coordinate extent to plot 
#     u,v,w - data[nlat,nlon] already selected at the right pressure level and time
#     below_mask[nlat,nlon] - indicates when a cell is below the surface for this pressure levels so we can mask out data
#     lons_uv,lats_uv,lons_w,lats_w - the coordinates of the data - u and v may have an extended grid compared to w
#                                     so lons and lats for both are passed in
#     track_lons, track_lats - the lats and lons of the BV track
#     track_tix - the index into track_lats/lons at this time
#     track_col - the colour in which to plot the position of the BV at this time (if '' dont plot the BV)
#     plot_track - whether to plot the whole BV track or not
# returns:
#     con_w - the contour plot object for w so colour bar can be added to the plot
#     sp - the streamline plot object for uv so colour bar can be added to the plot
#----------------------------------------------------------------------------------------------------------------------
def plot_uvw_map(fig, nrows,ncols, nplot, title, min_lon, max_lon, min_lat, max_lat, u, v, w, below_mask, lons_uv, lats_uv, lons_w, lats_w, track_lons, track_lats, track_tix, track_col, plot_track=False):

    ax = fig.add_subplot(nrows,ncols,nplot,projection=ccrs.PlateCarree())
    ax.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
    plt.title(title, loc='left')
    con_w = ax.contourf(lons_w, lats_w, w.data,levels=w_levels,extend='both', cmap=wind_cmap)
    wspeed=(u.data**2+v.data**2)**0.5
    sp = ax.streamplot(lons_uv, lats_uv, u.data, v.data,linewidth=linewidth,arrowsize = arr_sz, density=str_density, color=wspeed, cmap=stream_cmap, norm=stream_norm)
    ix_bel=np.where(below_mask)
    if len(ix_bel[0])>0:
        con_m = ax.contourf(lons_w, lats_w, below_mask,levels=[-1,0,1],cmap=mask_cmap)
    if track_col!='':
        ax.scatter(track_lons[track_tix], track_lats[track_tix], color=track_col,marker='*', zorder=4) # put marker at current centre
        if plot_track:
            ax.plot(track_lons, track_lats, color=track_col,zorder=4) # plor whole track
    map_formatter(ax,tick_base_x=tick_base_x,tick_base_y=tick_base_y,labelsize=8,res='50m')
    ax.set_ylim(min_lat,max_lat)
    ax.set_xlim(min_lon,max_lon)
    return con_w, sp

#-------------------------------------------------------------------------------------------------------------------
# plot a map of vorticity shaded and u and v as streamlines
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lon,max_lon,min_lat,max_lat - defines coordinate extent to plot 
#     u,v,vort, below_mask - data[nlat,nlon] already selected at the right pressure level and time
#            below_mask indicates when a cell is below the surface for this pressure levels so we can mask out data
#     lons_uv,lats_uv,lons_w,lats_w - the coordinates of the data - u and v may have an extended grid compared to w
#                                     so lons and lats for both are passed in
#     track_lons, track_lats - the lats and lons of the BV track
#     track_tix - the index into track_lats/lons at this time
#     track_col - the colour in which to plot the position of the BV at this time (if '' dont plot the BV)
#     plot_track - whether to plot the whole BV track or not
#----------------------------------------------------------------------------------------------------------------------
def plot_vort_map(fig, nrows,ncols, nplot, title, min_lon,max_lon,min_lat,max_lat, u, v, vort, below_mask, lons_uv, lats_uv, track_lons, track_lats, track_tix, track_col, plot_track=False):

    ax = fig.add_subplot(nrows,ncols,nplot,projection=ccrs.PlateCarree())
    ax.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
    plt.title(title, loc='left')
    con_v = ax.contourf(lons_uv, lats_uv, vort*vort_factor,levels=vort_levels,extend='both',cmap=pv_cmap)
    ax.scatter(track_lons[track_tix], track_lats[track_tix], color=track_col,marker='*', zorder=4) # put marker at current centre
    if plot_track:
        ax.plot(track_lons, track_lats, color=track_col,zorder=4) # plot whole track
    wspeed=(u.data**2+v.data**2)**0.5
    sp = ax.streamplot(lons_uv, lats_uv, u.data, v.data,linewidth=linewidth,arrowsize = arr_sz, density=str_density, color=wspeed, cmap=stream_cmap, norm=stream_norm)
    map_formatter(ax,tick_base_x=tick_base_x,tick_base_y=tick_base_y,labelsize=8,res='50m')
    ax.set_ylim(min_lat,max_lat)
    ax.set_xlim(min_lon,max_lon)
    return con_v, sp

#-------------------------------------------------------------------------------------------------------------------
# plot a map of rain shaded and u and v as streamlines
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lon,max_lon,min_lat,max_lat - defines coordinate extent to plot 
#     u,v,rain - data[nlat,nlon] already selected at the right pressure level and time
#     lons_uv,lats_uv,lons_p,lats_p - the coordinates of the data - u and v have an different grid compared to rain
#                                     so lons and lats for both are passed in
#     track_lons, track_lats - the lats and lons of the BV track
#     track_tix - the index into track_lats/lons at this time
#     track_col - the colour in which to plot the position of the BV at this time (if '' dont plot the BV)
#     plot_track - whether to plot the whole BV track or not
#----------------------------------------------------------------------------------------------------------------------
def plot_rain_map(fig, nrows,ncols, nplot, title, min_lon,max_lon,min_lat,max_lat, u, v, rain, lons_uv, lats_uv, lons_p, lats_p, track_lons, track_lats, track_tix,track_col, plot_track=False):
    ax = fig.add_subplot(nrows,ncols,nplot,projection=ccrs.PlateCarree())
    ax.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
    plt.title(title, loc='left')
    con_r = ax.contourf(lons_p, lats_p, rain.data,levels=rain_levels,extend='max')
    if plot_track:
        ax.plot(track_lons, track_lats, color=track_col,zorder=4) # plot whole track
    ax.scatter(track_lons[track_tix], track_lats[track_tix], color=track_col,marker='*', s=30, zorder=4) # put marker at current centre
    wspeed=(u.data**2+v.data**2)**0.5
    sp = ax.streamplot(lons_uv, lats_uv, u.data, v.data,linewidth=linewidth,arrowsize = arr_sz, density=str_density, color=wspeed, cmap=stream_cmap, norm=stream_norm)
    map_formatter(ax,tick_base_x=tick_base_x,tick_base_y=tick_base_y,labelsize=8,res='50m')
    ax.set_ylim(min_lat,max_lat)
    ax.set_xlim(min_lon,max_lon)
    return con_r, sp
        
#-------------------------------------------------------------------------------------------------------------------
# plot a map of pv shaded and u and v as streamlines
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     title - title of this plot
#     min_lon,max_lon,min_lat,max_lat - defines coordinate extent to plot 
#     u,v,pv - data[nlat,nlon] already selected at the right pressure level and time
#     lons_uv,lats_uv,lons_p,lats_p - the coordinates of the data - u and v have an different grid compared to rain
#                                     so lons and lats for both are passed in
#     track_lons, track_lats - the lats and lons of the BV track
#     track_tix - the index into track_lats/lons at this time
#     track_col - the colour in which to plot the position of the BV at this time (if '' dont plot the BV)
#     plot_track - whether to plot the whole BV track or not
#----------------------------------------------------------------------------------------------------------------------
def plot_pv_map(fig, nrows,ncols, nplot, title, min_lon,max_lon,min_lat,max_lat, u, v, pv, below_mask, lons_uv, lats_uv, lons_pv, lats_pv, track_lons, track_lats, track_tix, track_col, plot_track=False):
    ax = fig.add_subplot(nrows,ncols,nplot,projection=ccrs.PlateCarree())
    ax.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
    plt.title(title, loc='left')
    con_pv = ax.contourf(lons_pv, lats_pv, pv,levels=pv_levels,extend='both', cmap=pv_cmap)
    wspeed=(u.data**2+v.data**2)**0.5
    sp = ax.streamplot(lons_uv, lats_uv, u.data, v.data,linewidth=linewidth,arrowsize = arr_sz, density=str_density, color=wspeed, cmap=stream_cmap, norm=stream_norm)
    ax.scatter(track_lons[track_tix], track_lats[track_tix], color=track_col,marker='*', zorder=4) # put marker at current centre
    if plot_track:
        ax.plot(track_lons, track_lats, color=track_col, zorder=4) # plot whole track
    map_formatter(ax,tick_base_x=tick_base_x,tick_base_y=tick_base_y,labelsize=8,res='50m')
    ax.set_ylim(min_lat,max_lat)
    ax.set_xlim(min_lon,max_lon)
    return con_pv, sp
    

#-----------------------------------------------------------------------------------------------------
# plot time series of height and tilt - indicate when BV is in NH or SH
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     track850 - the bv_track object for this BV
#     surge_months, surge_days, surge_indices - the surge times and conditions
#     vortex_height, vortex_title_north, vortex_title_east - timeseries to plot
#     title - title of plot
#     flip - used to switch x and y axes round
#-----------------------------------------------------------------------------------------------------
def plot_timeseries_height_tilt(fig,nrows,ncols,nplot,track_850, surge_years, surge_months, surge_days, surge_indices, vortex_height, vortex_tilt_north, vortex_tilt_east,title, flip=False):

    # plot whole timeseries data for requested case:
    #   the hourly averaged rainfall and vorticity averaged over a +/-dlonlat degree circle centred on BV track as a timeseries with the surge on each day
    # if flip is true x axis is height/titlt and yaxis is time

    track_lats=track_850.vort_data[1,:]
    track_lons=track_850.vort_data[0,:]
    ix_n=np.where(track_lats>0)
    ix_s=np.where(track_lats<0)
    start_day=track_850.track_times[0].replace(hour=0,minute=0)
    end_day=track_850.track_times[-1].replace(hour=0,minute=0)
    ndays=(end_day-start_day).days+1
    hours=np.arange(ndays*4)*6 # this is 6 hourly
    dt850=(track_850.track_times[0]-start_day)
    first_hour850=dt850.days*24+dt850.seconds/3600
    first_hour850=int(first_hour850/6)
    hours850=hours[first_hour850:first_hour850+track_850.nt]
    offset_hour=12 
    six=np.where((surge_years==start_day.year) & (surge_months==start_day.month) & (surge_days==start_day.day))
    first_surge_tix=six[0][0]
    date_labels=[]
    ax = fig.add_subplot(nrows,ncols,nplot)
    plt.title(title, loc='left')
    if flip==True:
        ax.plot(vortex_height,hours850, color='k')
        ax.invert_xaxis()
        xlim=ax.get_xlim()
        for t in ix_n[0]:
            ax.text(xlim[1], hours850[t], 'N', color='red',verticalalignment='center')
        for t in ix_s[0]:
            ax.text(xlim[1], hours850[t], 'S', color='blue',verticalalignment='center')
        ax2 = ax.twiny()
        ax2.plot(vortex_tilt_north*100, hours850, linestyle='dotted',color='g', label='Tilt North')
        ax2.plot(vortex_tilt_east*100, hours850, linestyle='dashed', color='g', label='Tilt East')
        x=xlim[0]-0.01*(xlim[0]-xlim[1])
        for d in range(ndays):
            date_labels.append((start_day+dt.timedelta(days=d)).strftime('%d'))
            this_surge=surge_indices[d+first_surge_tix]
            ax.text(x, hours[d*4]+offset_hour, surge_text[this_surge], color='red', verticalalignment='center')
        ax.set_ylim(hours[0],hours[-1]+6)
        ax.set_ylabel('date')
        # plot labels only on the hour=0 dates
        ax.set_yticks(hours[::4])
        ax.set_yticklabels(date_labels)
        ax.set_xlabel('vortex height (hPa)')
        ax.tick_params(axis='x', colors='k')
        ax2.set_xlabel('vortex tilt (deg/Pa)', color='g')
        ax2.tick_params(axis='x', colors='g')
        ax2.set_xlim([-3,3])
    else:
        ax.plot(hours850, vortex_height, color='k')
        ax.invert_yaxis()
        ylim=ax.get_ylim()
        for t in ix_n[0]:
            ax.text(hours850[t], ylim[1], 'N', color='red',horizontalalignment='center')
        for t in ix_s[0]:
            ax.text(hours850[t], ylim[1], 'S', color='blue',horizontalalignment='center')
        ax.set_ylim(ylim)
        ax2 = ax.twinx()
        ax2.plot(hours850, vortex_tilt_north*100, linestyle='dotted',color='g', label='Tilt North')
        ax2.plot(hours850, vortex_tilt_east*100, linestyle='dashed', color='g', label='Tilt East')
        y=ylim[0]-0.01*(ylim[0]-ylim[1])
        for d in range(ndays):
            date_labels.append((start_day+dt.timedelta(days=d)).strftime('%d'))
            this_surge=surge_indices[d+first_surge_tix]
            ax.text(hours[d*4]+offset_hour, y, surge_text[this_surge], color='red', horizontalalignment='center')
        ax.set_xlabel('date')
        # plot labels only on the hour=0 dates
        ax.set_xticks(hours[::4])
        ax.set_xticklabels(date_labels)#, rotation=45)
        ax.set_ylabel('vortex height (hPa)')
        ax.tick_params(axis='y', colors='k')
        ax2.set_ylabel('vortex tilt (deg/Pa)', color='g')
        ax2.tick_params(axis='y', colors='g')
        ax2.set_ylim([-3,3])

    #ax2.legend()    
    return hours, hours850, date_labels

#-----------------------------------------------------------------------------------------------------
# plot time series of mean track rain and vorticity, indicate when BV is in NH or SH
# inputs:
#     fig - the matplotlib figure to add this subplot to
#     nrows, ncols, nplot - nplot is which plot out of the nrows x ncols
#     track850 - the bv_track object for this BV which contains the mean track rain and vorticity
#     surge_months, surge_days, surge_indices - the surge times and conditions
#     title - title of plot
#     flip - used to switch x and y axes round
#-----------------------------------------------------------------------------------------------------
def plot_timeseries_rain_vort(fig,nrows,ncols,nplot,track_850, surge_years, surge_months, surge_days, surge_indices,title, flip=False):
    track_lats=track_850.vort_data[1,:]
    track_lons=track_850.vort_data[0,:]
    ix_n=np.where(track_lats>0)
    ix_s=np.where(track_lats<0)
    start_day=track_850.track_times[0].replace(hour=0,minute=0)
    end_day=track_850.track_times[-1].replace(hour=0,minute=0)
    ndays=(end_day-start_day).days+1
    print(ndays, 'days')
    hours=np.arange(ndays*4)*6 # this is 6 hourly
    dt850=(track_850.track_times[0]-start_day)
    first_hour850=dt850.days*24+dt850.seconds/3600
    first_hour850=int(first_hour850/6)
    hours850=hours[first_hour850:first_hour850+track_850.nt]
    offset_hour=12
    six=np.where((surge_years==start_day.year) & (surge_months==start_day.month) & (surge_days==start_day.day))
    first_surge_tix=six[0][0]
    date_labels=[]
    ax = fig.add_subplot(nrows,ncols,nplot)
    plt.title(title,loc='left')
    if flip==True:
        ax.plot(track_850.imerg, hours850, color='blue')
        xlim=ax.get_xlim()
        ax.set_xlim(0,xlim[1])
        xlim=ax.get_xlim()
        for t in ix_n[0]:
            ax.text(xlim[1], hours850[t], 'N', color='red')
        for t in ix_s[0]:
            ax.text(xlim[1], hours850[t], 'S', color='blue')

        ax2 = ax.twiny()
        ax2.plot(track_850.mean_vort, hours850, color='k')
        x=xlim[0]+0.01*(xlim[1]-xlim[0])

        for d in range(ndays):
            date_labels.append((start_day+dt.timedelta(days=d)).strftime('%d'))
            this_surge=surge_indices[d+first_surge_tix]
            ax.text(x, hours[d*4]+offset_hour, surge_text[this_surge], color='red', horizontalalignment='center')
        ax.set_ylim(hours[0],hours[-1]+6)
        print('hours min max', hours[0], hours[-1], hours850[0], hours850[-1], 'ylim=',ax.get_ylim())

        ax.set_ylabel('date')
        # plot labels only on the hour=0 dates
        ax.set_yticks(hours[::4])
        ax.set_yticklabels(date_labels)
        ax.set_xlabel(Ra_label, color='blue')
        ax.tick_params(axis='x', colors='blue')
        ax2.set_xlabel(mean_vort_label, color='k')
        ax2.tick_params(axis='x', colors='k')
    else:
        ax.plot(hours850, track_850.imerg, color='blue')
        ylim=ax.get_ylim()
        for t in ix_n[0]:
            ax.text(hours850[t], ylim[1], 'N', color='red')
        for t in ix_s[0]:
            ax.text(hours850[t], ylim[1], 'S', color='blue')
        ax.set_ylim(ylim)
        ax2 = ax.twinx()
        ax2.plot(hours850, track_850.mean_vort, color='k')
        
        for d in range(ndays):
            date_labels.append((start_day+dt.timedelta(days=d)).strftime('%d'))
            this_surge=surge_indices[d+first_surge_tix]
            ax.text(hours[d*4]+offset_hour, ylim[0]+0.01*(ylim[1]-ylim[0]), surge_text[this_surge], color='red', horizontalalignment='center')

        ax.set_xlabel('date')
        # plot labels only on the hour=0 dates
        ax.set_xticks(hours[::4])
        ax.set_xticklabels(date_labels)#, rotation=45)
        ax.set_ylabel(Ra_label, color='blue')
        ax.tick_params(axis='y', colors='blue')
        ax2.set_ylabel(mean_vort_label, color='k')
        ax2.tick_params(axis='y', colors='k')
    
    return hours, hours850, date_labels
