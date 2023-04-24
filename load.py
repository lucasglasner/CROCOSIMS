#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 12:52:26 2022

@author: lucasg
"""



# ---------------------------------------------------------------------------- #
# ---------------------------------- Imports --------------------------------- #
# ---------------------------------------------------------------------------- #

import xarray as xr
from glob import glob
import numpy as np
import pandas as pd
import datetime
from utils import fix_crocotime,center_crocogrid


# ---------------------------------------------------------------------------- #
# ---------------------------- loading data functions ------------------------ #
# ---------------------------------------------------------------------------- #


def load_croco(paths, YORIG, variables=None, **kwargs):
    """
    Loading function for reading raw CROCO/ROMS outputs into
    xarray objects. **kwargs are passed to xarray.open_dataset() function

    Args:
        paths (str, list): _description_
        YORIG (str): _description_
        variables (list, optional): _description_. Defaults to None.

    Raises:
        ValueError: If input path is not a string or list of strings

    Returns:
        xarray: loaded croco data
    """
    if type(paths) == str:
        if "*" in paths:
            data = xr.open_mfdataset(paths,
                                    concat_dim='time',
                                    combine='nested',
                                    **kwargs)
            data = fix_crocotime(data, YORIG)
        else:
            data = xr.open_dataset(paths, **kwargs)
            data = fix_crocotime(data, YORIG)
    elif type(paths) == list:
        data=[]
        for p in paths:
            d = xr.open_dataset(p, **kwargs)
            d = fix_crocotime(d,YORIG)
            data.append(d)
        data = xr.concat(data,'time')
    else:
        raise ValueError('Path to croco output should only be a string or a list of strings.')
    if variables==None:
        data = center_crocogrid(data, data.keys())
    else:
        data = center_crocogrid(data, variables)[variables]
    return data.sortby('time')




