#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 16:01:37 2022

@author: stem
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 00:53:52 2021

@author: stem
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
#import matplotlib.pylab as pl
#from matplotlib.colors import ListedColormap
#from netCDF4 import Dataset
from matplotlib.collections import PolyCollection
#import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from shapely.geometry import Polygon
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm
#%%
simname = 'quickplot'
plotdir =  '/project/ivme/MCH-1/icon-art-BRM/icon_output%s'%(simname)
#%%
"""Create path"""
if not os.path.exists(plotdir):
  os.makedirs(plotdir)
#%%
datapath = '/project/ivme/MCH-1/icon-art-BRM/icon_output/'
filename = 'ICON-ART-OEM_1504NewICUnstr_DOM01_00000000.nc'
gridpath = '/project/ivme/MCH-1/icon-art-BRM/icon_output/'
gridfile = 'ICON-1E_DOM01.nc'

fh_d = Dataset(os.path.join(datapath,filename))
fh_g = Dataset(os.path.join(gridpath,gridfile))
    
datavar = np.ma.array(fh_d.variables['clct'][:])
    
clat = np.array(fh_g.variables["clat"][:]) * 180 / np.pi
vlon = np.array(fh_g.variables["vlon"][:]) * 180 / np.pi
vlat = np.array(fh_g.variables["vlat"][:]) * 180 / np.pi
voc = np.array(fh_g.variables["vertex_of_cell"][:])
#%%
corners = np.zeros((len(clat),3,2))
for ncell in np.arange(corners.shape[0]):
    vert_ind = voc[:,ncell] - 1
    for ic in np.arange(3):
        corners[ncell,ic,:] = vlon[vert_ind[ic]] , vlat[vert_ind[ic]]
#%%
if len(datavar.shape)==2:
    plotvar = datavar[0,...]
elif len(datavar.shape)==4:
    plotvar = datavar[0,-1,...]

minval = np.nanmin( [ np.nanquantile(plotvar, 0.005) , np.nanquantile(plotvar, 0.055) ] )
maxval = np.nanmax( [ np.nanquantile(plotvar, 0.995) , np.nanquantile(plotvar, 0.995) ] )
#%%
cmap = cm.get_cmap('viridis')


colors_s = np.ma.array(plotvar, mask=plotvar.mask)


coll = PolyCollection(corners,array=colors_s,cmap=cmap,
                      edgecolors='none',linewidth=0.)
coll.cmap.set_bad(color='white',alpha=0)
coll.set_array(colors_s)
coll.set_clim(vmin=minval,vmax=maxval)
    
    
#%%
"""countries"""
tsize=37


fig = plt.figure(figsize=(50,29))


ax1 = fig.add_subplot(1,1,1)
ax1.add_collection(coll)
ax1.autoscale_view()
ax1.set_title('ICON-OUTPUT of clct at 0',y=1.01,fontsize=tsize+6)
ax1.set_xticks(np.arange(-2,19,5))
ax1.set_xlim(-2,19)
ax1.set_yticks(np.arange(32,51,5))
ax1.set_ylim(42,51)
ax1.tick_params(axis="both",labelsize=tsize)
plt.grid()
plt.xlabel('longitude [°E]',fontsize=tsize)
plt.ylabel('latitude [°N]',fontsize=tsize)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2i'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2i'))


fig.subplots_adjust(top=0.90)
fig.subplots_adjust(left=0.05)
fig.subplots_adjust(right=0.95)


figname = os.path.join(plotdir,'clct_NewICBC_0')
fig.savefig(figname+'.png')
plt.close(fig)

