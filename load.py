#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 12:52:26 2022

@author: lucasg
"""


import xarray as xr
from glob import glob
import numpy as np
import pandas as pd
import datetime

def fix_crocotime(DATA,YORIG):
    ORIG = pd.to_datetime(YORIG)
    new_time = pd.to_datetime([datetime.timedelta(seconds=t.item())+ORIG
                               for t in DATA.time])
    DATA['time'] = new_time
    return DATA.sortby('time')

def center_grid(data,variables):
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
    lat,lon = data.lat_rho.values[:,0],data.lon_rho.values[0,:]-360
    data = data.assign_coords({'lat':('eta_rho',lat),'lon':('xi_rho',lon)})
    data = data.swap_dims({'xi_rho':'lon','eta_rho':'lat'})
    data = data.drop(['xi_rho','eta_rho', 'lon_rho', 'lat_rho'])   
    return data

def load_croco(paths, YORIG='1958-01-01T12:00:00', **kwargs):
    if type(paths) == str:
        data = xr.open_mfdataset(paths,
                                 concat_dim='time',
                                 combine='nested',
                                 **kwargs)
        data = fix_crocotime(data,YORIG)
    elif type(paths) == list:
        data=[]
        for p in paths:
            d = xr.open_dataset(p, **kwargs)
            d = fix_crocotime(d,YORIG)
            data.append(d)
        data = xr.concat(data,'time')
    else:
        raise ValueError('Path to croco output should only be a string or a list of strings.')
    data = center_grid(data, data.keys())
    return data


def rhopoints_depths(h,zeta,s_rho,Cs_r,hc,vtransform=2):
    if vtransform==1:
        Z_rho = hc*(s_rho-Cs_r)+Cs_r*h
        z_rho = Z_rho+zeta*(1+Z_rho/h)
        return z_rho
    if vtransform==2:
        Z_rho = (hc*s_rho+Cs_r*h)/(hc+h)
        z_rho = zeta+(zeta+h)*Z_rho
        return z_rho
    return
    

def haversine(p1,p2):
    lat1,lon1 = p1
    lat2,lon2 = p2
    
    lon1,lon2,lat1,lat2 = map(np.deg2rad, [lon1,lon2,lat1,lat2])
    
    dlon = lon2-lon1
    dlat = lat2-lat1
    
    a = np.sin(dlat/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    c = 2*np.arcsin(np.sqrt(a))
    r = 6371
    return c*r
    





