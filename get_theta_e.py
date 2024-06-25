# convert temperature to potential temp and equivalent potential temp
import numpy as np
import iris

def get_potential_temperature(T_cube, plev_name='pressure_level'):

    print('getting potential temp')
    R_over_Cp=0.286
    P_zero=1000.0

    p_coord=T_cube.coord(plev_name)
    pressure=p_coord.points
    nplev=len(pressure)
    t_coord=T_cube.coord('time')
    shape=np.shape(T_cube.data)
    if len(shape)==4:
        nt=shape[0]
    else:
        nt=0
    theta=np.zeros_like(T_cube.data)
    for pix in np.arange(nplev):
        if nt==0:
            this_temp=T_cube[pix,:,:].data
        else:
            this_temp=T_cube[:,pix,:,:].data
        if np.any(this_temp ==np.inf):
            pdb.set_trace()
        this_theta=this_temp*(P_zero/pressure[pix])**R_over_Cp
        if nt==0:
            theta[pix,:,:]=this_theta
        else:
            theta[:,pix,:,:]=this_theta

    lat_coord=T_cube.coord('latitude')
    lon_coord=T_cube.coord('longitude')
    if nt==0:
        dim_coords_and_dims=[(p_coord, 0), (lat_coord, 1), (lon_coord,2)]
    else:
        dim_coords_and_dims=[(t_coord, 0), (p_coord, 1), (lat_coord, 2), (lon_coord,3)]
    theta_cube=iris.cube.Cube(theta, var_name='theta', long_name='potential temperature', units=T_cube
.units, dim_coords_and_dims=dim_coords_and_dims)
 
    return theta_cube
  
def get_equivalent_potential_temperature(T_cube, q_cube, plev_name='pressure_level'):
    print('getting equivalent potential temp')
    L_v = 2.25e6  # latent heat of vaporization of water, J kg-1
    cp_d = 1004  # heat capacity of dry air at constant pressure, J kg-1 K-1
    p_coord=T_cube.coord(plev_name)

    theta_cube=get_potential_temperature(T_cube)
    equiv_theta=theta_cube.data * np.exp(((L_v/cp_d) * q_cube.data)/(T_cube.data * (1000-q_cube.data)))
    if np.any(equiv_theta==np.inf):
        pdb.set_trace()
    lat_coord=T_cube.coord('latitude')
    lon_coord=T_cube.coord('longitude')
    t_coord=T_cube.coord('time')
    shape=np.shape(T_cube.data)
    if len(shape)==4:
        nt=shape[0]
    else:
        nt=0
    if nt==0:
        dim_coords_and_dims=[(p_coord, 0), (lat_coord, 1), (lon_coord,2)]
    else:
        dim_coords_and_dims=[(t_coord, 0), (p_coord, 1), (lat_coord, 2), (lon_coord,3)]
    equiv_theta_cube=iris.cube.Cube(equiv_theta, var_name='theta_e', long_name='equivalent potential temperature', units=T_cube.units, dim_coords_and_dims=dim_coords_and_dims)

    return equiv_theta_cube

