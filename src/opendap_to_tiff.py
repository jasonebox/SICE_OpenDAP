# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:15:26 2023

@author: rabni

pip install opendap-protocol
conda update --all

"""

import xarray as xr
from pyproj import CRS,Transformer
import os
import numpy as np
from rasterio.transform import Affine
import rasterio as rio

if os.getlogin() == 'jason':
    base_path = '/Users/jason/0_dat/S3/opendap/'

os.chdir(base_path)

WGSProj = CRS.from_string("+init=EPSG:4326")
PolarProj = CRS.from_string("+init=EPSG:3413")
wgs_data = Transformer.from_proj(WGSProj, PolarProj) 

def ExportGeoTiff(x,y,z,crs,path):
    
    "Input: xgrid,ygrid, data paramater, the data projection, export path, name of tif file"
    
    resx = (x[1] - x[0])
    resy = (y[1] - y[0])
    transform = Affine.translation((x[0]),(y[0])) * Affine.scale(resx, resy)
    
    with rio.open(
    path,
    'w',
    driver='GTiff',
    height=z.shape[0],
    width=z.shape[1],
    count=1,
    dtype=z.dtype,
    crs=crs,
    transform=transform,
    ) as dst:
        dst.write(z, 1)
    
    return None 


# a list of dates, here it is set to 2022-09-21
dates = ['2022-09-' + str(dd+1).zfill(2) for dd in np.arange(20,21)]
out_f = os.getcwd() # output folder for tif files


# This is the variables plotted, it has to be in a Dict format like: 
#  plotting_dict[variable] = {'minval' : n, 'maxval' : n,\
#                                        'cm' : colormap}      
#  minval is the min value of the colorscale, maxval is the max value of the colorscale. 
# values outside the scale will not be plotted, cm is the colormap.

# plotting_dict = {}

# plotting_dict['albedo_bb_planar_sw'] = {'minval' : 0, 'maxval' : 1}

# plotting_dict['isnow'] = {'minval' : 1, 'maxval' : 3}
# plotting_dict['factor'] = {'minval' : 0, 'maxval' : 1}

vars_i_want=['albedo_bb_planar_sw','factor']

# area of interest in espg:4326
lat_n = 62
lat_s = 61
lon_w = -48.5
lon_e = -46.5
lon_e= -46.515023
lon_w= -48.354069
lat_n= 62.047626
lat_s= 60.863512


############### This is the Code ##################

west_x,north_y = wgs_data.transform(lon_w, lat_n)
east_x,south_y = wgs_data.transform(lon_e, lat_s)

y_slice = slice(int(north_y),int(south_y))
x_slice = slice(int(west_x),int(east_x))


## how to list all variables
for d in dates:
    DATASET_ID = 'sice_500_' + d.replace('-', '_') + '.nc'
    
    for v in vars_i_want:
        # try:
            
        print(DATASET_ID,v)
        ds = xr.open_dataset(f'https://thredds.geus.dk/thredds/dodsC/SICE_Greenland_500m/{DATASET_ID}')
        # a=np.array(list(ds.keys()))
        # a.sorted()
        yshape,xshape = np.shape(ds[v])
        
        if yshape != 5424: 
            ds = ds.rename({'x2':'xcoor'})
            ds = ds.rename({'y2':'ycoor'})        
        else:
            ds = ds.rename({'x':'xcoor'})
            ds = ds.rename({'y':'ycoor'})
            
        data = ds[v].sel(ycoor=y_slice,xcoor=x_slice)
        x = ds['xcoor'].sel(xcoor=x_slice)
        y = ds['ycoor'].sel(ycoor=y_slice)
        
   
        data = data.where(data <= plotting_dict[v]['maxval'])
        data = data.where(data >= plotting_dict[v]['minval'])
        
        z = data.to_numpy()
        x = x.to_numpy()
        y = y.to_numpy()
        
        path = out_f + os.sep + DATASET_ID[:-3] + '_' + v + '.tif'
        
        ExportGeoTiff(x, y, z, PolarProj, path)
        
        ds.close()
    
        # except:
            
        #     print('data does not exist')
