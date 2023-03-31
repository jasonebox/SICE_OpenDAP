# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:15:26 2023

@author: rabni
"""

import xarray as xr
from pyproj import CRS,Transformer
WGSProj = CRS.from_string("+init=EPSG:4326")
PolarProj = CRS.from_string("+init=EPSG:3413")
wgs_data = Transformer.from_proj(WGSProj, PolarProj) 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.colors import ListedColormap

#%%
plotting_dict = {}

# This is a colormap i made for the isnow diagnostic
cm_isnow = {1:"blue",
          2:"green",
          3:"red"
          }
cm_isnow_rgb = ListedColormap([cm_isnow[x] for x in cm_isnow.keys()])
labels_isnow = np.array(["Clean Snow","Dirty Snow","Bare Ice"])

# This is a list of dates (right now it is set to 2022-09-21):
dates = ['2022-09-' + str(dd+1).zfill(2) for dd in np.arange(20,21)]

# This is the variables plotted, it has to be in a Dict format like: 
#  plotting_dict['albedo_bb_planar_sw'] = {'minval' : n, 'maxval' : n,\
#                                        'cm' : colormap}      
# the minval is min value of the colorscale, maxval is the max value of the colorscale. 
# values outside the scale will not be plotted, cm is the colormap.

 
plotting_dict['albedo_bb_planar_sw'] = {'minval' : 0, 'maxval' : 1,\
                                        'cm' : plt.get_cmap('Blues_r')}

plotting_dict['isnow'] = {'minval' : 1, 'maxval' : 3,\
                          'cm' : cm_isnow_rgb,'cm_dict' : cm_isnow, 'labels' : labels_isnow}

plotting_dict['factor'] = {'minval' : 0, 'maxval' : 1,\
                           'cm' : plt.get_cmap('Greys_r') }

# This is the Box you want data from in espg:4326
lat_N = 62
lat_S = 61
lon_W = -48.5
lon_E = -46.5


############### This is the Code ##################

# This is the extent of map plotting in PlateCarree Coordinates which xarray uses

#Gr_W = -60
#Gr_E = -30.5
#Gr_S = 58
#Gr_N = 85

Gr_W = -48.5
Gr_E = -46.5
Gr_S = 62
Gr_N = 62.8


west_x,north_y = wgs_data.transform(lon_W, lat_N)
east_x,south_y = wgs_data.transform(lon_E, lat_S)


Gr_W_x,Gr_N_y = wgs_data.transform(Gr_W, Gr_N)
Gr_E_x,Gr_S_y = wgs_data.transform(Gr_E, Gr_S)


y_slice = slice(int(north_y),int(south_y))
x_slice = slice(int(west_x),int(east_x))


for d in dates:
    
    DATASET_ID = 'sice_500_' + d.replace('-', '_') + '.nc'
    
    for v in plotting_dict:
        try:
            
            print(DATASET_ID)
            filex=f'https://thredds.geus.dk/thredds/dodsC/SICE_Greenland_500m/{DATASET_ID}'
            ds = xr.open_dataset(filex)
                
            yshape,xshape = np.shape(ds[v])
            
            if yshape != 5424: 
                ds = ds.rename({'x2':'xcoor'})
                ds = ds.rename({'y2':'ycoor'})        
            else:
                ds = ds.rename({'x':'xcoor'})
                ds = ds.rename({'y':'ycoor'})
                
                
            data = ds[v].sel(ycoor=y_slice,xcoor=x_slice)
       
# plot it
            plt.figure(figsize=[80,32])
            ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
            ax.set_global()
            
            ax.set_extent([Gr_W, Gr_E, Gr_S, Gr_N], ccrs.PlateCarree())
            ax.tick_params(labelsize=18)
            #ax.coastlines()
            
            ax.gridlines(color='black', alpha=0.5, linestyle='--')
            
            data = data.where(data <= plotting_dict[v]['maxval'])

            np_data=np.array(data.to_numpy())
            plt.imshow(np_data)
        except:
            
            print('data does not exist')

            #%%
            # data.variables
            np_data=np.array(data.to_numpy())
            plt.imshow(np_data)
            
            import rasterio
            
            profile=np_data.profile
            islands_data=islands.read(1)
            #islands_data+=260
            #sectors_data=rasterio.open(ofile).read(1)
            #v=(islands_data!=260)
            #sectors_data[v]=islands_data[v]
            
            wo=0
            if wo:
                with rasterio.open(ofile, 'w', **profile) as dst:
                    dst.write(dst, 1)
            
            #%%
            p = data.plot(x='xcoor', y='ycoor',
                      vmin = plotting_dict[v]['minval'],
                      vmax = plotting_dict[v]['maxval'],
                      cmap=plotting_dict[v]['cm'],
                      ax=ax,
                      add_labels=False,
                      add_colorbar=False
                      )
            
            plt.title(v, size=60)
            #add separate colorbar
            if v =='isnow':
                cm_isnow = plotting_dict[v]['cm_dict']
                labels = plotting_dict[v]['labels'] 
                cm_isnow_rgb = ListedColormap([cm_isnow[x] for x in cm_isnow.keys()])
                len_lab = len(labels)

                norm_bins = np.sort([*cm_isnow.keys()]) + 0.5
                norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
                print(norm_bins)
                ## Make normalizer and formatter
                norm = matplotlib.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
                fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
                diff = norm_bins[1:] - norm_bins[:-1]
                tickz = norm_bins[:-1] + diff / 2
                cb = plt.colorbar(p, format=fmt, ticks=tickz)
                cb.ax.tick_params(labelsize=40)
            else:
                     
                cb = plt.colorbar(p, shrink=0.99)
                cb.ax.tick_params(labelsize=40)                    

            
            
            gl = p.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=4, color='gray', alpha=0.5, linestyle='--')
            
            gl.xlabels_bottom = False
            gl.ylabels_right = False
            gl.ylabels_top = False
            gl.ylocator = mticker.FixedLocator([60,65,70,75,80,85])
            gl.xlocator = mticker.FixedLocator([-80, -40, 0])
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 40, 'color': 'black'}
            gl.ylabel_style = {'size': 40, 'color': 'black'}
            
            
            #ax.set_extent([Gr_W_x,Gr_E_x, Gr_S_y, Gr_N_y], ccrs.NorthPolarStereo())
            
            
            plt.show()
            
            #ds.close()
    
        except:
            
            print('data does not exist')
