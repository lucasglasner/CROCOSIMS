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
    


def fill_borders(data):
    """
    Fill all nans with forward and backward filling.

    Args:
        data (xarray): 

    Returns:
        xarray: 
    """
    #fill bays with pixels from the top
    data = data.bfill('lat', limit=4)
    data = data.ffill('lon').bfill('lon')
    data = data.ffill('lat').bfill('lat')
    return data

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

def compute_metrics(model,reference, dim='time'):
    """
    Compute model metrics (MBIAS, RMSE, pearsonr) against reference.

    Args:
        model xarray: _description_
        reference xarray: _description_
        dim (str, optional): dimension to reduce along. Defaults to 'time'.

    Returns:
        xarray: Dataset with stored metrics
    """
    BIAS  = model-reference
    MBIAS = BIAS.mean(dim)
    RMSE  = (BIAS**2).mean(dim)**0.5
    CORR  = xr.corr(model,reference,dim)
    METRICS = xr.merge([MBIAS.to_dataset(name='MBIAS'),
                        RMSE.to_dataset(name='RMSE'),
                        CORR.to_dataset(name='CORR')])
    return METRICS


def bias_correct_SST(data,method='linregress'):
    """
    Bias correct SST with delta method or linregress
    Require the files with the coefficient or correction data.

    Args:
        data (xarray): 
         Mercator forecast sea surface temperature
        method (str, optional): 
         linregress or delta. Defaults to 'linregress'.

    Returns:
        xarray: bias corrected mercator sst as xarray object
    """
    if method == 'delta':
        p = '~/storage/FORECAST/MERCATOR/'
        mbias = xr.open_dataset(p+'MBIAS_MERCATORANALYSIS-OSTIA_MONTHLY.nc').mbias
        mbias = mbias.rename({'longitude':'lon','latitude':'lat'})
        mbias = mbias.reindex({'lon':data.lon,'lat':data.lat},
                            method='nearest')
        data = data.groupby('time.dayofyear')-mbias
    if method == 'linregress':
        p = '~/storage/FORECAST/MERCATOR/'
        lr = xr.open_dataset(p+'LINEAR_REGRESSION_COEFFICIENTS.nc')
        lr = lr.reindex({'lon':data.lon,'lat':data.lat},
                            method='nearest')
        data = data*lr.SLOPE+lr.INTERCEPT
    return data


def deg2compass(angle):
    """
    Transform an angle in degrees to the
    windrose string.
    The angle increases clockwise with 0° pointing North.

    Args:
        angle (float): Angle in degrees

    Returns:
        str: vector direction string
    """
    if angle<0:
        angle = angle+360
    if np.isnan(angle):
        return np.nan
    val = int((angle/45)+0.5)
    arr = np.array(['N','NE','E','SE','S','SW','W','NW'])
    return arr[(val%8)]

def coastwinddir(u,v,land_sea_mask):
    """
    Compute wind angle against coastline

    Args:
        u (xarray): zonal wind
        v (xarray): meridional wind
        land_sea_mask (xarray): Land sea mask

    Returns:
        angle: raster with the angle of the wind respect to coastline
    """
    LSM = land_sea_mask
    LSM = (LSM.where(LSM==1).ffill('lon')==1).astype(float)
    LSM_GRADIENT=LSM.differentiate('lon'),LSM.differentiate('lat')
    LSM_GRADIENT_MODULE=(LSM_GRADIENT[0]**2+LSM_GRADIENT[1]**2)**0.5    
    ws = (u**2+v**2)**0.5
    angle = LSM_GRADIENT[0]*u+LSM_GRADIENT[1]*v
    angle = angle/(LSM_GRADIENT_MODULE*ws)*180/np.pi
    angle = angle.where(~LSM.astype(bool))
    return angle

def beaufort_scale(wind):
    """
    Given the wind (knots) return the beaufort number
    
    Args:
        wind (float): wind in knots

    Returns:
        int: beaufort number
    """
    scale=np.array([0,1,4,7,11,17,22,28,34,41,48,56,64,np.inf])
    for n in range(len(scale)-1):
        if scale[n]<=wind and scale[n+1]>wind:
            return int(n)
    

def seasonal_decompose(ts, period, nharmonics=3, bandwidth=2):
    """
    Parameters
    ----------
    ts : Time series data in a pandas series format, with timestamps
         in the index.
    period : period of the season
    nharmonics : Number of harmonics to remove, default is 3.

    Returns
    -------
    season : Seasonal component of the time series.
    anomaly : The time series anomaly without the seasonal cycle.
    """
    if len(ts)%2==0:
        n = len(ts)
    else:
        n = len(ts)+1
    ft = np.fft.fft(ts)
    ft[0] = 0  # Remove mean#
    for i in range(nharmonics):  # Filter cycle#
        pos = n//(period//(i+1))
        ft[pos-bandwidth:pos+bandwidth] = 0
        ft[n-pos-bandwidth:n-pos+bandwidth] = 0
        # ft[pos]=0
        # ft[n-pos]=0
    anomaly = np.fft.ifft(ft).real
    # anomaly = pd.Series(anomaly, index=ts.index)
    season = ts-anomaly
    return season

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

def compute_anomalies(forecast, climatology, timename='leadtime'):
    """Compute anomaly from climatology

    Args:
        forecast (XDataArray): forecast data
        climatology (XDataArray): climatology data

    Returns:
        XDataArray: anomaly
    """
    climatology = climatology.reindex({'lat':forecast.lat,'lon':forecast.lon})
    anomaly = forecast.groupby(timename+'.dayofyear')-climatology
    return anomaly

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

def grabpoint(data,lat,lon,method='nearest'):
    """
    Given a coordinate in lat,lon an a netcdf file in data (as
    an xarray object), this function grabs the timeseries
    in the specified coordinate.

    Args:
        data (xarray): _description_
        lat (float): _description_
        lon (float): _description_
        method (str, optional): Interpolation method. Defaults to 'nearest'.

    Returns:
        pd.Dataframe: timeseries of the pixel as a dataframe
    """
    point = fill_borders(data).interp(lat=lat,lon=lon, method=method)
    point = point.to_dataframe()
    if isinstance(point.index,pd.DatetimeIndex):
        point = utc_to_local(point)
    return point