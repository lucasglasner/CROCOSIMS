'''
 # @ Author: Your lucas
 # @ Create Time: 2022-07-22 19:35:24
 # @ Modified by: lucas
 # @ Modified time: 2022-07-22 19:35:46
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
# ---------------------------- Numerical data functions ---------------------- #
# ---------------------------------------------------------------------------- #


def rhopoints_depths(h,zeta,s_rho,Cs_r,hc,vtransform=2):
    """ Compute depth of roms sigma coordinates.

    Args:
        h (_type_): _description_
        zeta (_type_): _description_
        s_rho (_type_): _description_
        Cs_r (_type_): _description_
        hc (_type_): _description_
        vtransform (int, optional): _description_. Defaults to 2.

    Returns:
        _type_: _description_
    """
    if vtransform==1:
        Z_rho = hc*(s_rho-Cs_r)+Cs_r*h
        z_rho = Z_rho+zeta*(1+Z_rho/h)
        return z_rho
    else:
        Z_rho = (hc*s_rho+Cs_r*h)/(hc+h)
        z_rho = zeta+(zeta+h)*Z_rho
        return z_rho
    

def haversine(p1,p2):
    """
    Given two points with lat,lon coordinates, compute the distance 
    between those points on the surface of the sphere with the haversine formula
    Args:
        p1 (tuple): first point lat,lon
        p2 (tuple): last point lat,lon

    Returns:
        float: distance
    """
    lat1,lon1 = p1
    lat2,lon2 = p2
    
    lon1,lon2,lat1,lat2 = map(np.deg2rad, [lon1,lon2,lat1,lat2])
    
    dlon = lon2-lon1
    dlat = lat2-lat1
    
    a = np.sin(dlat/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    c = 2*np.arcsin(np.sqrt(a))
    r = 6371
    return c*r


def regrid(data,husk,method='bilinear'):
    """
    Regrid dataset to a new latlon grid.

    Args:
        data (xarray): Data to regrid
        husk (xarray): Data with the new coordinates.
        method (str, optional):
         "bilinear","conservative","nearests2d"
         Defaults to 'bilinear'.

    Returns:
        (xarray): Regridded data using xESMF routines.
    """
    regridder = xe.Regridder(data,husk,method)
    return regridder(data)

def filter_timeseries(ts, order, cutoff, btype='lowpass', fs=1, **kwargs):
    """Given an array, this function apply a butterworth (high/low pass) 
    filter of the given order and cutoff frequency.
    For example:
    If 'ts' is a timeseries of daily samples, filter_timeseries(ts,3,1/20)
    will return the series without the 20 days or less variability using an
    order 3 butterworth filter. 
    In the same way, filter_timeseries(ts,3,1/20, btype='highpass') will
    return the series with only the 20 days or less variability.

    Args:
        ts (array_like): timeseries or 1D array to filter
        order (int): _description_
        cutoff (array_like): Single float for lowpass or highpass filters, 
        arraylike for bandpass filters.
        btype (str, optional): The type of filter. Defaults to 'lowpass'.
        fs (int): Sampling frequency. Defaults to 1.s
        **kwargs are passed to scipy.signal.butter

    Returns:
        output (array): Filtered array
    """
    mask = np.isnan(ts)
    nans = np.ones(len(ts))*np.nan
    if mask.sum()==len(ts):
        return nans
    else:
        b, a = signal.butter(order,cutoff, btype=btype, fs=fs, **kwargs)
        filt=signal.filtfilt(b, a, ts[~mask])
        output=np.ones(len(ts))*np.nan
        output[np.where(~mask)] = filt
        return output
    
def filter_xarray(data, dim, order, cutoff, btype='lowpass', parallel=True, fs=1):
    """Given a 3d DataArray, with time and spatial coordinates, this function apply
    the 1D function filter_timeseries along the time dimension, filter the complete
    xarray data.

    Args:
        data (XDataArray): data
        dim (str): name of the time dimension
        order (int): butterworth filter order
        cutoff (array_like): if float, the cutoff frequency, if array must be the
                            [min,max] frequencys for the bandpass filter.
        btype (str, optional): {lowpass,highpass,bandpass}. Defaults to 'lowpass'.
        parallel (bool, optional): If parallelize with dask. Defaults to True.
        fs (int, optional): Sampling frequency. Defaults to 1.

    Returns:
        XDataArray: filtered data
    """
    if parallel:
        dask='parallelized'
    else:
        dask='forbidden'
    filt = xr.apply_ufunc(filter_timeseries, data, order, cutoff, btype, fs,
                          input_core_dims=[[dim],[],[],[],[]],
                          output_core_dims=[[dim]],
                          exclude_dims=set((dim,)),
                          keep_attrs=True,
                          vectorize=True, dask=dask)
    filt[dim] = data[dim]
    return filt


def utc_to_local(series, gap=4):
    """
    Transform pandas series with datetime index from UTC time to local

    Args:
        series (pd.Series): pandas timeseries with utc data 

    Returns:
        pd.Series: data in local time
    """
    series.index = series.index-pd.Timedelta(hours=gap)
    return series

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
            x = data[v].interp(eta_v=data.eta_rho.values)
            x = x.rename({'eta_v':'eta_rho'})
            data = data.drop(v)
            data[v]=x
        elif 'xi_u' in data[v].dims:
            x = data[v].interp(xi_u=data.xi_rho.values)
            x = x.rename({'xi_u':'xi_rho'})
            data = data.drop(v)
            data[v]=x
        else:
            pass
    data = data.drop(['xi_u','eta_v', 'lon_v', 'lat_v','lon_u','lat_u'])
    return data
