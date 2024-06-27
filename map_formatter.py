########################################################################
# Selection of plotting tools specific to my needs. Contents below.    #
#                                              John Ashcroft, Oct 2017 #
########################################################################

import cartopy.crs as ccrs
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm

import iris
from matplotlib.colors import from_levels_and_colors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#import matplotlib.colors as colors


#### Contents:   #####
#
## map_formatter()
# Draws gridlines, coastlines, ticks and formats labels on axes. 
# tick_base defines the interval in degrees between points on the grid.
#
def map_formatter(ax,tick_base_x=15.0, tick_base_y=15.0, labelsize=20,top_label=False, bottom_label=True, right_label=False,left_label=True,res='10m', coast_col='k'):

    ax.coastlines(resolution=res, color=coast_col)

    x_ticks=ax.get_xticks()
    y_ticks=ax.get_yticks()
    x_lim=ax.get_xlim()
    if (x_ticks[1]-x_ticks[0])!= tick_base_x:
        x_ticks[0]=int(x_ticks[0]/tick_base_x)*tick_base_x
        x_ticks=np.arange(x_ticks[0], x_lim[-1]+1, tick_base_x)
    y_lim=ax.get_ylim()
    if (y_ticks[1]-y_ticks[0])!= tick_base_y:
        y_ticks[0]=int(y_ticks[0]/tick_base_y)*tick_base_y
        y_ticks=np.arange(y_ticks[0], y_lim[-1]+1, tick_base_y)
    ax.set_xticks(x_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(y_ticks, crs=ccrs.PlateCarree())
    #gl=ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.75, color='k', linestyle=':')
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    if top_label:
        ax.xaxis.set_tick_params(labeltop=True,labelsize=labelsize)
    else:
        ax.xaxis.set_tick_params(labeltop=False)
    if right_label:
        ax.yaxis.set_tick_params(labelright=True,labelsize=labelsize)
    else:
        ax.yaxis.set_tick_params(labelright=False)
    if left_label:
        ax.yaxis.set_tick_params(labelleft=True,labelsize=labelsize)
    else:
        ax.yaxis.set_tick_params(labelleft=False)
    if bottom_label:
        ax.xaxis.set_tick_params(labelbottom=True,labelsize=labelsize)
    else:
        ax.xaxis.set_tick_params(labelbottom=False)
    
import matplotlib.ticker as mticker
def map_formatter_old(ax,tick_base_x=15.0, tick_base_y=15.0, labelsize=20,top_label=False, bottom_label=True, right_label=False,left_label=True,central_longitude=0.0,res='10m', coast_col='k'):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='k', linestyle=':')
    gl.xlabels_top = top_label
    gl.ylabels_right = right_label
    gl.ylabels_left = left_label
    gl.xlabels_bottom = bottom_label        
    gl.xlocator = mticker.MultipleLocator(base=tick_base_x)
    gl.ylocator = mticker.MultipleLocator(base=tick_base_y)
    ax.coastlines(resolution=res, color=coast_col, linewidth=1)
    gl.xlabel_style= {'size':labelsize}
    gl.ylabel_style= {'size':labelsize}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
        

