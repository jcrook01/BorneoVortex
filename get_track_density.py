import warnings
import numpy as np
from read_kevin_bv_tracks import *
from map_formatter import *
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
import pdb

# code to calculate and plot track density
density_levels=np.arange(40)*0.02
density_title='tracks per 30 day'

def get_track_density(bv_tracks, lons, lats,ntimes):

    lon_range=[np.amin(lons),np.amax(lons)]
    lat_range=[np.amin(lats),np.amax(lats)]

    nlon=len(lons)
    nlat=len(lats)
    density=np.zeros((nlat,nlon))
    bvs_in_reg=0  # number of BVs
    for i in range(len(bv_tracks)):
        tix=bv_tracks[i].exists_in_region(lon_range, lat_range)
        if len(tix[0])>0:            
            bvs_in_reg=bvs_in_reg+1
        for t in tix[0]:
            this_lon=bv_tracks[i].vort_data[0,t]
            this_lat=bv_tracks[i].vort_data[1,t]
            ix=np.where(abs(this_lon-lons)<=0.5)
            iy=np.where(abs(this_lat-lats)<=0.5)
            if len(ix[0])>0 and len(iy[0])>0:
                #print(this_lon, this_lat, lons[ix[0][0]],lats[iy[0][0]])
                density[iy[0][0],ix[0][0]]=density[iy[0][0],ix[0][0]]+1

    ix=np.where(density>0)
    ix=np.where(density==0)
    density[ix]=np.nan

    return density/ntimes, bvs_in_reg

# plots density as a map - can add wind vectors too if u and v passed in
def plot_density(fig, nrows, ncols, plot_num, data, levels, cmap, title, lons, lats, extend='max', u=[], v=[]):

    min_lon=85
    max_lon=135
    min_lat=-10
    max_lat=20
    ax = fig.add_subplot(nrows,ncols,plot_num,projection=ccrs.PlateCarree())
    ax.set_extent([min_lon,max_lon,min_lat,max_lat],ccrs.PlateCarree())
    con=ax.contourf(lons,lats,data,cmap=cmap, levels=levels, extend=extend)
    if len(u)!=0 and len(v)!=0:
        wind_key_x=0.85
        wind_key_y=0.96
        wind_key_u=1
        wind_key_txt=' = 1 m/s'
        q1=ax.quiver(lons[::2],lats[::2],u[::2,::2], v[::2,::2])
        ax.quiverkey(q1, X=wind_key_x, Y=wind_key_y, U=wind_key_u,label=wind_key_txt, labelpos='E')
    
    map_formatter(ax,tick_base_x=10,tick_base_y=5,labelsize=10,res='50m')
    plt.title(title, loc='left')
    ax.set_ylim(min_lat,max_lat)
    ax.set_xlim(min_lon,max_lon)

    return con, ax


if __name__ == '__main__':
    main()
