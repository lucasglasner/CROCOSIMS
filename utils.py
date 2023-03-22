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
import scipy.signal as signal

# ---------------------------------------------------------------------------- #
# ---------------------------- Utilities functions --------------------------- #
# ---------------------------------------------------------------------------- #

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
    data = data.sortby('eta_rho').sortby('xi_rho')
    geoindex = ((data.lon_rho > lonmin) & (data.lon_rho<lonmax) & (data.lat_rho>latmin) & (data.lat_rho < latmax)).load()
    geoindex = np.argwhere(geoindex.values)
    xmin = min(geoindex[:,1])
    xmax = max(geoindex[:,1])
    ymin = min(geoindex[:,0])
    ymax = max(geoindex[:,0])
    print(xmin,xmax)
    print(ymin,ymax)
    data = data.sel(eta_rho=slice(ymin,ymax), xi_rho=slice(xmin,xmax))
    return data
