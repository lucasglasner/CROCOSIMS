'''
 # @ Author: Your name
 # @ Create Time: 2023-03-22 09:01:45
 # @ Modified by: Your name
 # @ Modified time: 2023-03-22 09:01:49
 # @ Description:
 '''


# ---------------------------------------------------------------------------- #
# ---------------------------------- Imports --------------------------------- #
# ---------------------------------------------------------------------------- #

import xesmf as xe
import xarray as xr
import numpy as np
import pandas as pd
import datetime

# ---------------------------------------------------------------------------- #
# ---------------------------- Utilities functions --------------------------- #
# ---------------------------------------------------------------------------- #


def fill_borders(data, x_name='lon', y_name='lat'):
    """
    Fill all nans with forward and backward filling.

    Args:
        data (xarray): 

    Returns:
        xarray: 
    """
    #fill bays with pixels from the top
    data = data.bfill(y_name, limit=4)
    data = data.ffill(x_name).bfill(x_name)
    data = data.ffill(y_name).bfill(y_name)
    return data

def fix_crocotime(DATA,YORIG):
    """
    Grab simulation time and transform to datetime objects based on given YORIG

    Args:
        DATA (XDataset, XDataArray): CROCO simulation data, with "time" coordinate 
        YORIG (str): Given reference date

    Returns:
        XDataset, XDataArray: Data with fixed time coordinate
    """
    ORIG = pd.to_datetime(YORIG)
    new_time = pd.to_datetime([datetime.timedelta(seconds=t.item())+ORIG
                               for t in DATA.time])
    DATA['time'] = new_time
    return DATA.sortby('time')


def croco_sellonlatbox(data, lonmin, lonmax, latmin, latmax):
    """
    This functions grabs a croco output and slice it to the
    desired latlon bounds. 
    Only works for rho point centered variables (like temp, salt, zeta, etc)
    Other arakawa-C grid variables like zonal and meridional currents must first
    be horizontally interpolated or shifted to the grid center.

    Args:
        data (xarray): loaded dataset (centered in rho points) as xarray object
        lonmin (float): min longitude
        lonmax (float): max longitude
        latmin (float): min latitude
        latmax (float): max latitude

    Returns:
        data: sliced data to the user defined latlon box 
    """
    geoindex = ((data.lon_rho > lonmin) & (data.lon_rho<lonmax) & (data.lat_rho>latmin) & (data.lat_rho < latmax)).load()
    geoindex = np.argwhere(geoindex.values)
    xmin = min(geoindex[:,1])
    xmax = max(geoindex[:,1])
    ymin = min(geoindex[:,0])
    ymax = max(geoindex[:,0])
    # print(xmin,xmax)
    # print(ymin,ymax)
    data = data.sel(eta_rho=slice(ymin,ymax), xi_rho=slice(xmin,xmax))
    return data

def center_crocogrid(data,variables):
    """
    This function grabs a croco model outputs and moves the 
    arakawa-C edges variables to the center of the grid.
    (like water velocities...) 

    Args:
        data (xarray): input dataset
        variables (list): list of variables to transform

    Returns:
        xarray: dataset with variables in the center of the grid 
    """
    for v in variables:
        if 'eta_v' in data[v].dims:
            x = data[v].interp(eta_v=data.eta_rho.values, method='nearest')
            x = x.rename({'eta_v':'eta_rho'})
            data = data.drop(v)
            data[v]=x
        elif 'xi_u' in data[v].dims:
            x = data[v].interp(xi_u=data.xi_rho.values, method='nearest')
            x = x.rename({'xi_u':'xi_rho'})
            data = data.drop(v)
            data[v]=x
        else:
            pass
    data = data.drop(['xi_u','eta_v', 'lon_v', 'lat_v','lon_u','lat_u'])
    return data

def croco_selpoint(data, lon, lat):
    """
    This functions finds the nearest point for the given lat,lon 
    coordinates.

    Args:
        data (_type_): _description_
        lon (_type_): _description_
        lat (_type_): _description_

    Returns:
        _type_: _description_
    """
    eta = abs(data.lat_rho-lat).argmin(axis=0)[0]
    xi  = abs(data.lon_rho-lon).argmin(axis=1)[0]
    return data.sel(eta_rho=eta, xi_rho=xi)


def trim_memory() -> int:
    """
    Following dask docs: 
    "It is possible to forcefully release allocated but unutilized memory as follows:
    This should be only used as a one-off debugging experiment. Watch the dashboard while running the above code.
    If unmanaged worker memory (on the “Bytes stored” plot) decreases significantly after calling client.run(trim_memory),
    then move on to the next section. Otherwise, you likely do have a memory leak.""
    How to use:
        client.run(trim_memory)
    
    """
    import ctypes
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)