import warnings
import numpy as np
import pdb
import math

R_earth = 6371200.

#-----------------------
# get_3d_lat_lon_array
#--------------------------
def get_3d_lev_lat_lon_array(levs,lats,lons):
    nlev=len(levs)
    nlat=len(lats)
    nlon=len(lons)
    levarr = np.zeros((nlev,nlat,nlon))
    dlevarr = np.zeros((nlev,nlat,nlon))
    dlevs = np.zeros_like(levs)
    dlevs[1:-1] = (levs[2:] - levs[:-2])/2
    dlevs[0] = (levs[1] - levs[0])
    dlevs[-1] = (levs[-1] - levs[-2])
    latarr = np.zeros((nlev,nlat,nlon))
    dlatarr = np.zeros((nlev,nlat,nlon))
    dlats = np.zeros(nlat)
    dlats[1:-1] = (lats[2:] - lats[:-2])/2  # use central differences
    dlats[0] = (lats[1] - lats[0])  # use 1st differences at boundary
    dlats[-1] = (lats[-2] - lats[-1])
    lonarr = np.zeros((nlev,nlat,nlon))
    dlonarr = np.zeros((nlev,nlat,nlon))
    dlons = np.zeros(nlon)
    dlons[1:-1] = (lons[2:] - lons[:-2])/2  # use central differences
    dlons[0] = (lons[1] - lons[0])  # use 1st differences at boundary
    dlons[-1] = (lons[-2] - lons[-1])
    for jj in range(0,nlat):
        for ii in range(0,nlon):
            levarr[:,jj,ii] = levs
            dlevarr[:,jj,ii] = dlevs
            latarr[:,jj,ii] = lats[jj]
            dlatarr[:,jj,ii] = dlats[jj]
            lonarr[:,jj,ii] = lons[ii]
            dlonarr[:,jj,ii] = dlons[ii]

    return levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr

#-----------------------
# get_2d_lat_lon_array
#--------------------------
def get_2d_lat_lon_array(lats,lons):
    nlat=len(lats)
    nlon=len(lons)
    latarr = np.zeros((nlat,nlon))
    dlatarr = np.zeros((nlat,nlon))
    dlats = np.zeros(nlat)
    dlats[1:-1] = (lats[2:] - lats[:-2])/2  # use central differences
    dlats[0] = (lats[1] - lats[0])  # use 1st differences at boundary
    dlats[-1] = (lats[-2] - lats[-1])
    lonarr = np.zeros((nlat,nlon))
    dlonarr = np.zeros((nlat,nlon))
    dlons = np.zeros(nlon)
    dlons[1:-1] = (lons[2:] - lons[:-2])/2  # use central differences
    dlons[0] = (lons[1] - lons[0])  # use 1st differences at boundary
    dlons[-1] = (lons[-2] - lons[-1])
    for jj in range(0,nlat):
        for ii in range(0,nlon):
            latarr[jj,ii] = lats[jj]
            dlatarr[jj,ii] = dlats[jj]
            lonarr[jj,ii] = lons[ii]
            dlonarr[jj,ii] = dlons[ii]

    return latarr, dlatarr, lonarr, dlonarr


#-----------------------
    """
    Return the gradient of a 3-dimensional array on a sphere given a pressure, latitude
    and longitude arrays.
    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.
    Parameters
    ----------
    f : A 3-dimensional array containing samples of a scalar function on pressure, lat and lons.
    levsarr, dlevarr, latsarr, dlatarr, lonsarr, dlonarr: computed in get_3d_lev_lat_lon_array
          from pressure vector, latitude vector and longitude vector
    axis - if set only compute the gradient on the given axis, otherwise compute in 3 D

    Returns
    -------
    g : dfdz,dfdy,dfdx arrays of the same shape as `f` giving the derivative of `f` with
        respect to each dimension.
    """
#-----------------------
def gradient_sphere(f, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=None):

    dlatsrad = np.deg2rad(dlatarr)
    dlonsrad = np.deg2rad(dlonarr)
    latrad = np.deg2rad(latarr)
    dy = R_earth * np.cos(latrad) * dlatsrad
    dx = R_earth * np.cos(latrad) * dlonsrad
    if axis==0:
       df = np.gradient(f,axis=axis)/dlevarr
    elif axis==1:
       df = np.gradient(f*np.cos(latrad),axis=axis)/dy
    elif axis==2:
       df = np.gradient(f,axis=axis)/dx
    else:
       dfdz = np.gradient(f,axis=0)/dlevarr
       dfdy = np.gradient(f*np.cos(latrad),axis=1)/dy
       dfdx = np.gradient(f,axis=2)/dx
       df= [dfdz,dfdy,dfdx]

    return df

def get_horizontal_divergence(u_cube, v_cube, latarr, dlatarr, lonarr, dlonarr):

    dlatsrad = np.deg2rad(dlatarr)
    dlonsrad = np.deg2rad(dlonarr)
    latrad = np.deg2rad(latarr)
    dy = R_earth * np.cos(latrad) * dlatsrad
    dx = R_earth * np.cos(latrad) * dlonsrad
    shape=np.shape(u_cube.data)
    if len(shape)==2:
        xaxis=1
        yaxis=0
    elif len(shape)==3:
        xaxis=2
        yaxis=1
    dvdy = np.gradient(v_cube.data*np.cos(latrad),axis=yaxis)/dy
    dudx = np.gradient(u_cube.data,axis=xaxis)/dx

    return dudx+dvdy
#-----------------------
# get_vorticity
#
# inputs:
#    u_cube, v_cube with dimensions [nz,nlat,nlon] - the horizontal wind arrays
#     where nz can be time array or pressure levels
#
# returns:
#    vertical relative vorticity calculated in spherical corrdinates
#-----------------------
def get_vorticity(u_cube,v_cube, plev_name='pressure_level'):

    plevs=u_cube.coord(plev_name).points
    lats=u_cube.coord('latitude').points
    lons=u_cube.coord('longitude').points
    times=u_cube.coord('time').points
    nt=len(times)
    nplev=len(plevs)
    shape=np.shape(u_cube.data)
    dim_unknown=True
    if len(shape)==2:
        # add third dimension
        udata=np.zeros((1,shape[0],shape[1]))
        udata[0,:,:]=u_cube.data
        vdata=np.zeros((1,shape[0],shape[1]))
        vdata[0,:,:]=v_cube.data
        levs=np.arange(3)
        dim_unknown=False
        levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr=get_3d_lev_lat_lon_array(levs,lats,lons)
        dudy=gradient_sphere(udata, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=1)
        dvdx=gradient_sphere(vdata, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=2)
        vort=dvdx[0,:,:]-dudy[0,:,:]
    else:
        if len(shape)==3:
            if nplev>1 and nplev==shape[0]:
                levs=plevs*100
                dim_unknown=False
            elif nt>1 and nt==shape[0]:
                levs=times
                dim_unknown=False
            if dim_unknown==False:
                udata=u_cube.data
                vdata=v_cube.data   
                levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr=get_3d_lev_lat_lon_array(levs,lats,lons)
                dudy=gradient_sphere(udata, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=1)
                dvdx=gradient_sphere(vdata, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=2)
                vort=dvdx-dudy
        elif len(shape)==4:
            vort=np.zeros(shape)
            if nplev>1 and nplev==shape[1] and nt==shape[0]:    
                levs=plevs*100
                levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr=get_3d_lev_lat_lon_array(levs,lats,lons)
                udata=u_cube.data
                vdata=v_cube.data
                dim_unknown=False
                # need to loop round all times as 1st index
                for t in range(nt):
                    dudy=gradient_sphere(udata, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=1)
                    dvdx=gradient_sphere(vdata, levarr, dlevarr, latarr, dlatarr, lonarr, dlonarr, axis=2)
                    vort[t,:,:,:]=dvdx-dudy
                      
    if dim_unknown==True:
        print('which dimensions to use?')
        print(u_cube)
        pdb.set_trace()
    return vort,lons,lats

#-----------------------
# get_potential_vorticity
#-----------------------
def get_potential_vorticity(u_cube, v_cube, theta_cube, plev_name='pressure_level'):

    plevs=u_cube.coord(plev_name).points
    lats=u_cube.coord('latitude').points
    lons=u_cube.coord('longitude').points
    shape=np.shape(u_cube.data)
    if len(shape)==4:
        nt=shape[0]
    else:
        nt=0
    plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr=get_3d_lev_lat_lon_array(plevs*100,lats,lons)
    # Coriolis parameter
    f = 2*(7.292e-05)*np.sin(np.deg2rad(latarr))
    g=9.81
    pv=np.zeros_like(u_cube.data)
    if nt==0:
        dtheta=gradient_sphere(theta_cube.data, plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr)
        du=gradient_sphere(u_cube.data, plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr)
        dv=gradient_sphere(v_cube.data, plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr)
        pv=-10**6*g*(-dv[0]*dtheta[2]+du[0]*dtheta[1]+dv[2]*dtheta[0]-du[1]*dtheta[0] +f*dtheta[0])
    else:
      for t in range(nt):
        dtheta=gradient_sphere(theta_cube[t,:,:,:].data, plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr)
        du=gradient_sphere(u_cube[t,:,:,:].data, plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr)
        dv=gradient_sphere(v_cube[t,:,:,:].data, plevarr, dplevarr, latarr, dlatarr, lonarr, dlonarr)
        pv[t,:,:,:]=-10**6*g*(-dv[0]*dtheta[2]+du[0]*dtheta[1]+dv[2]*dtheta[0]-du[1]*dtheta[0] +f*dtheta[0])

    return pv
