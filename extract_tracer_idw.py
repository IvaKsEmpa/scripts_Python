#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 19:38:13 2021

@author: stem
"""


#%%
"""Import"""
from shapely.geometry import Point, Polygon
from netCDF4 import Dataset
import os
import numpy as np
import datetime
from multiprocessing import Pool
from itertools import repeat
from math import sin, cos, sqrt, atan2, radians

#%%

"""Interpolation function"""
    
def intp_icon_data(args):
    
    iloc = args[0]
    gridinfo = args[1]
    datainfo = args[2]
    latitudes = args[3]
    longitudes = args[4]
    asl = args[5]

    nn_sel = np.zeros(gridinfo.nn)
    u=np.zeros(gridinfo.nn)
    
    R = 6373.0 # approximate radius of earth in km
    
    
    if (radians(longitudes[iloc])<np.nanmin(gridinfo.clon)) or (radians(longitudes[iloc])>np.nanmax(gridinfo.clon)):
        u[:] = np.nan
        return np.zeros((gridinfo.nn)), np.zeros((gridinfo.nn)).astype(int), np.zeros((gridinfo.nn)).astype(int), nn_sel[:], u[:]
    
    
    if (radians(latitudes[iloc])<np.nanmin(gridinfo.clat)) or (radians(latitudes[iloc])>np.nanmax(gridinfo.clat)):
        u[:] = np.nan
        return np.zeros((gridinfo.nn)), np.zeros((gridinfo.nn)).astype(int), np.zeros((gridinfo.nn)).astype(int), nn_sel[:], u[:]

    #%
    lat1 = radians(latitudes[iloc])
    lon1 = radians(longitudes[iloc])
    
    #%
    """FIND 4 CLOSEST CENTERS"""
    distances = np.zeros((len(gridinfo.clon)))
    for icell in np.arange(len(gridinfo.clon)):
        lat2 = gridinfo.clat[icell]
        lon2 = gridinfo.clon[icell]
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        distances[icell] = R * c
    nn_sel[:] = [x for _, x in sorted(zip(distances,np.arange(len(gridinfo.clon))))][0:gridinfo.nn]
    nn_sel=nn_sel.astype(int)
    #print('---nn_sel:',nn_sel)
    #print('---distances[0:gridinfo.nn]:',distances[0:gridinfo.nn])
    u[:] = [1./distances[y] for y in nn_sel]
    print('---distances:',[distances[y] for y in nn_sel])
    print('---weights:',u)
    
    #%
    """Calculate vertical interpolation factor"""
    idx_above = -1*np.ones((len(nn_sel))).astype(int)
    idx_below = -1*np.ones((len(nn_sel))).astype(int)
    for nnidx in np.arange(len(nn_sel)):
        for i_mc,mc in enumerate(datainfo.z_mc[:,nn_sel[nnidx]]):
            if mc>=asl[iloc]:
                idx_above[nnidx] = i_mc
            else:
                idx_below[nnidx] = i_mc
                break
    
    
    #in case of data point below lowest midlevel:
    for nnidx in np.arange(len(nn_sel)):
        if idx_below[nnidx]==-1:
            idx_below[nnidx] = idx_above[nnidx]
    if any(np.ravel([idx_above,idx_below])<0):
        print("At least one nearest neigbor has no valid height levels")
        u[:] = np.nan #-999.
        return np.zeros((gridinfo.nn)), np.zeros((gridinfo.nn)).astype(int), np.zeros((gridinfo.nn)).astype(int), nn_sel[:], u
    
    vert_scaling_fact = np.zeros((len(nn_sel)))
    for nnidx in np.arange(len(nn_sel)):
        if idx_below[nnidx] != idx_above[nnidx]:
            vert_scaling_fact[nnidx] = (asl[iloc]-datainfo.z_mc[idx_below[nnidx],nn_sel[nnidx]])/(datainfo.z_mc[idx_above[nnidx],nn_sel[nnidx]]-datainfo.z_mc[idx_below[nnidx],nn_sel[nnidx]])
        else:
            vert_scaling_fact[nnidx] = 0.
                                                                                              
    print('---idx above:',idx_above)
    print('---idx below:',idx_below)
    print('---vert_scaling_fact:',vert_scaling_fact)
    
    #%

    return vert_scaling_fact, idx_below, idx_above, nn_sel[:], u


#%%
    
"""Paths, Variables and Constants"""

ICON_grid = '/scratch/snx3000/msteiner/synthetic_data/input/grid/dyn_grid.nc'
stationdir = '/store/empa/em05/msteiner/CTDAS_data/observations/input_2021'
fname_base = 'ICON-ART-OEM-UNSTR'
DATA_path = '/scratch/snx3000/msteiner/ICON_simulations_synthetic/output_synth_mc_2018010100'
starttime ='2018010200'
enddtime = '2018011200'
nlev = 60 #number of vertical levels
nneighb = 4 #number of nearest neighbors to consider

#%%

"""Variables to extract"""

meta = { 
        'CH4_A': {'offset': 1.2e-06},
        'CH4_BG': {'offset': 0., 'ensemble': 8},
        'u': {'offset': 0.},
        'v': {'offset': 0.},
        'temp': {'offset': 0.},
        'CH4_A-ENS': {'offset': 1.2e-06, 'ensemble': 33},
        }

#%%

"""Get grid data"""

fh_grid = Dataset(ICON_grid,'r')
class gridinfo:
    clon_vertices = np.array(fh_grid.variables['clon_vertices'])
    clat_vertices = np.array(fh_grid.variables['clat_vertices'])
    cells_of_vertex = np.array(fh_grid.variables['cells_of_vertex'])
    vertex_of_cell = np.array(fh_grid.variables['vertex_of_cell'])
    neighbor_cell_index = np.array(fh_grid.variables['neighbor_cell_index'])
    vlon = np.array(fh_grid.variables['vlon'])
    vlat = np.array(fh_grid.variables['vlat'])
    clon = np.array(fh_grid.variables['clon'])
    clat = np.array(fh_grid.variables['clat'])
    ncells = len(fh_grid.dimensions['cell'])
    nn=nneighb

#%%
    
"""Times"""

firstfile=True
startdate = datetime.datetime.strptime(starttime,'%Y%m%d%H')
enddate = datetime.datetime.strptime(enddtime,'%Y%m%d%H')
delta = datetime.timedelta(hours=1)
looptime = startdate

#%%

"""Get locations of measurement stations"""

longitudes = []
latitudes = []
stationnames = []
asl = []
#for station in stationlist:
for filename in os.listdir(stationdir):
    if not filename.endswith('verify.csv'):
        continue
#    print('Processing %s ...'%(filename))

    infile = os.path.join(stationdir,filename)
    f = open(infile, 'r')
    lines = f.readlines()
    f.close()
    ns = filename.split('_')
    newname = ns[1]+'_'+ns[2]
    for iline,line in enumerate(lines):
      if iline<11: continue
      items = line.split(',')
      if (items[0]!='IZO' and items[0]!='ZEP'):
          latitudes.append(float(items[10]))
          longitudes.append(float(items[11]))
#          stationnames.append(items[0])
          stationnames.append(newname)
          asl.append(float(items[-5]))
      if iline==11:
          break
print("Found %i locations."%(len(latitudes)))

"""Add 5 missing stations"""
missing_longitudes = [-1.15, 4.93, 0.23, 8.4, 8.18]
missing_latitudes = [54.36, 51.97, 50.98, 47.48, 47.19]
missing_stationnames = ['bsd', 'cbw', 'hea', 'lae', 'beo']
missing_asl = [628.,200.,250.,872.,1009.]
for imiss in np.arange(len(missing_longitudes)):
    latitudes.append(missing_latitudes[imiss])
    longitudes.append(missing_longitudes[imiss])
    stationnames.append(missing_stationnames[imiss])
    asl.append(missing_asl[imiss])
print("Added %i missing locations."%(len(missing_longitudes)))

#%%

"""Initialize output variables"""
n_det = int(np.nansum([1 for var in meta.keys() if 'ensemble' not in meta[var]]))
n_ens = int(np.nansum([1 for var in meta.keys() if 'ensemble' in meta[var]]))
intp_ICON_data_det = np.zeros((n_det,len(latitudes),0))
maxmem=0
for var in meta.keys():
    if 'ensemble' in meta[var]:
        if meta[var]['ensemble']>maxmem:
            maxmem=meta[var]['ensemble']
maxmem=int(maxmem)
intp_ICON_data_ens = np.zeros((n_ens,maxmem,len(latitudes),0))
#%%

"""Loop over Data Files (=timesteps)"""

datetime_list = []
print('======================================')
date_idx = 0
while looptime <= enddate:

    intp_ICON_data_det = np.concatenate(( intp_ICON_data_det,np.zeros((n_det,len(latitudes),1)) ),axis=2)
    intp_ICON_data_ens = np.concatenate(( intp_ICON_data_ens,np.zeros((n_ens,maxmem,len(latitudes),1)) ),axis=3)
            
    timestring = datetime.datetime.strftime(looptime,'%Y-%m-%dT%H')
    datetime_list.append(timestring)
    DATA_file = os.path.join(DATA_path,'%s_%s:00:00.000.nc' %(fname_base,timestring))

    print('extracting from %s'%(DATA_file), flush=True)

    
    fh_data = Dataset(DATA_file,'r')
    
    class datainfo:
        z_mc = np.array(fh_data.variables['z_mc'])
        z_ifc = np.array(fh_data.variables['z_ifc'])
        
        
    ICON_data_det = np.zeros(( n_det, nlev, gridinfo.ncells ))
    ICON_data_ens = np.zeros(( n_ens, maxmem, nlev, gridinfo.ncells ))
    ivar = 0
    for var in meta.keys():
        if not 'ensemble' in meta[var]:
            ICON_data_det[ivar,...] = np.array(fh_data.variables[var]) - meta[var]['offset']
            ivar+=1
    ivar = 0
    for var in meta.keys():
        if 'ensemble' in meta[var]:
#            for iens in np.arange(n_member):
            for iens in np.arange(meta[var]['ensemble']):
                varnc = var.split('-')[0]+'-%.3i'%(iens+1)
                ICON_data_ens[ivar,iens,...] = np.array(fh_data.variables[varnc]) - meta[var]['offset']
            ivar+=1
#%%
        
    """Since the stations don't walk around, I only call the function at the first timestep"""  
    
    if looptime==startdate:
        args = zip( np.arange(len(latitudes)),repeat(gridinfo),repeat(datainfo),repeat(latitudes),repeat(longitudes),repeat(asl) )
        with Pool(48) as pool:
            vsf, idxb, idxa, neighbours, u_ret = list(zip(*pool.map(intp_icon_data, args)))

        
    vsf = np.array(vsf)
    idxb = np.array(idxb, dtype=int)
    idxa = np.array(idxa, dtype=int)
    neighbours = np.array(neighbours, dtype=int)
    u_ret = np.array(u_ret)

    #Do the interpolation 
    for iloc in np.arange(len(latitudes)):

        ##First, the vertical interpolation:
        vert_intp_data = np.zeros(( n_det, len(idxb[iloc]) ))
        for nn in np.arange(len(idxb[iloc])):
            vert_intp_data[:,nn] = ICON_data_det[:,idxb[iloc,nn],neighbours[iloc,nn]] + vsf[iloc,nn]*(ICON_data_det[:,idxa[iloc,nn],neighbours[iloc,nn]] \
                                      -ICON_data_det[:,idxb[iloc,nn],neighbours[iloc,nn]])
        ##Now the horizontal interpolation:
        ###First, for the det values:
        intp_ICON_data_det[:,iloc,date_idx] = np.nansum([w*vert_intp_data[:,i] for i,w in enumerate(u_ret[iloc,:])],axis=0)/np.nansum(u_ret[iloc,:]) 
        ###Second, the ensemble values:
        vert_intp_data = np.zeros(( n_ens, maxmem, len(idxb[iloc]) ))
        for nn in np.arange(len(idxb[iloc])):
            vert_intp_data[:,:,nn] = ICON_data_ens[:,:,idxb[iloc,nn],neighbours[iloc,nn]] + vsf[iloc,nn]*(ICON_data_ens[:,:,idxa[iloc,nn],neighbours[iloc,nn]] \
                                      -ICON_data_ens[:,:,idxb[iloc,nn],neighbours[iloc,nn]])
        intp_ICON_data_ens[:,:,iloc,date_idx] = np.nansum([w*vert_intp_data[:,:,i] for i,w in enumerate(u_ret[iloc,:])],axis=0)/np.nansum(u_ret[iloc,:])
#%%
    """Update time"""
    looptime += delta
    date_idx += 1
#%%
"""Save as netcdf"""
with Dataset('output/extracted_%s.nc'%(DATA_path.split('/')[-1],), mode='w') as ofile:

    osites = ofile.createDimension('sites', len(latitudes))
    otime = ofile.createDimension('time', (date_idx))

    oname = ofile.createVariable('site_name', str, ('sites'))
    otimes = ofile.createVariable('time', np.unicode_, ('time'))
   
    ivar = 0
    for var in meta.keys():
        if 'ensemble' not in meta[var]:
            ovar = ofile.createVariable(var, np.float32, ('sites','time'))
            ovar[:,:] = intp_ICON_data_det[ivar,:,:]
            ivar+=1
    ivar=0
    for var in meta.keys():
        if 'ensemble' in meta[var]:
            oens = ofile.createDimension('ens_%.2i'%(ivar+1), meta[var]['ensemble'])
            varnc = var.split('-')[0]+'_ENS'
            ovar = ofile.createVariable(varnc, np.float32, ('ens_%.2i'%(ivar+1),'sites','time'))
            ovar[:,:,:] = intp_ICON_data_ens[ivar,0:meta[var]['ensemble'],:,:]
            ivar+=1

    
    oname[:] = np.array(stationnames[:])
    otimes[:] = np.array(datetime_list)
