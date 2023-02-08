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


#%%

"""Interpolation function"""
    
def intp_icon_data(args):
    
    iloc = args[0]
    gridinfo = args[1]
    datainfo = args[2]
    latitudes = args[3]
    longitudes = args[4]
    asl = args[5]

    nn_sel = np.zeros(3)
    u=np.zeros(3)
    
    
    
    if (longitudes[iloc]<np.nanmin(gridinfo.clon)) or (longitudes[iloc]>np.nanmax(gridinfo.clon)):
        u[:] = np.nan
        return 0., 0, 0, nn_sel[:], u
    
    
    if (latitudes[iloc]<np.nanmin(gridinfo.clat)) or (latitudes[iloc]>np.nanmax(gridinfo.clat)):
        u[:] = np.nan
        return 0., 0, 0, nn_sel[:], u

    #%
    """Create shapely point"""
    Ploc = Point(longitudes[iloc],latitudes[iloc])
    #%
    """Determin containing triangle"""
    # make sure we iterate only over valid gridpoint indices
    for nc in np.arange(gridinfo.ncells):
        if Ploc.x < min(gridinfo.clon_vertices[nc,:]):
            continue
        if Ploc.x > max(gridinfo.clon_vertices[nc,:]):
            continue
        if Ploc.y < min(gridinfo.clat_vertices[nc,:]):
            continue
        if Ploc.y > max(gridinfo.clat_vertices[nc,:]):
            continue
        icon_cell = Polygon(list(zip(gridinfo.clon_vertices[nc,:],gridinfo.clat_vertices[nc,:])))
        if icon_cell.contains(Ploc):
            break
        else:
            continue
    #%
    """Calculate vertical interpolation factor"""
    idx_above = -1
    idx_below = -1
    for i_mc,mc in enumerate(datainfo.z_mc[:,nc]):
        if mc>=asl[iloc]:
            idx_above = i_mc
        else:
            idx_below = i_mc
            break
    
    
    #in case of data point below lowest midlevel:
    if idx_below==-1:
        idx_below = idx_above
    if any(np.array([idx_above,idx_below])<0):
        print("Found one location outide domain (or in boundary region)")
        u[:] = np.nan #-999.
        return 0., 0, 0, nn_sel[:], u
    
    
    if idx_below != idx_above:
        vert_scaling_fact = (asl[iloc]-datainfo.z_mc[idx_below,nc])/(datainfo.z_mc[idx_above,nc]-datainfo.z_mc[idx_below,nc])
    else:
        vert_scaling_fact = 0.
        
    #%
    """find the correct corner
    (unfortunately it's not simply the closest corner, but the one corner where the distance between
    data point and corner is shorter than the distance between center point and corner)"""
    nv = 0
    for iv in np.arange(3):
        nvloop = gridinfo.vertex_of_cell[iv,nc]
        distance_dc = Ploc.distance(Point(gridinfo.vlon[nvloop-1],gridinfo.vlat[nvloop-1])) #python-idx
        distance_cc = Point(gridinfo.clon[nc],gridinfo.clat[nc]).distance(Point(gridinfo.vlon[nvloop-1],gridinfo.vlat[nvloop-1])) #python-idx
        if distance_dc<distance_cc:
            nv = nvloop
            break
        
    
    nv = nv - 1 #python-idx
    #%
    """Find the 3 cells surrounding this closest edge"""
    nn = np.zeros(6)
    for icell in np.arange(6):
        nn[icell] = gridinfo.cells_of_vertex[icell,nv]-1
    nn=nn.astype(int)
    
    #find neighbours of nc in the 6 cells surrounding the closest vertex:
    neighb_c = np.zeros(2)
    ie=0
    for icell in nn:
        icell = int(icell)
        if icell in (gridinfo.neighbor_cell_index[:,nc]-1):
            neighb_c[ie] = icell
            ie+=1
            
    #the neighbours of those 2 neighbours in the surrounding cells are the one we look for
    neighb_c=neighb_c.astype(int)
    ie=0
    for icell in nn:
        icell = int(icell)
        if icell in (gridinfo.neighbor_cell_index[:,neighb_c[0]]-1) or icell in (gridinfo.neighbor_cell_index[:,neighb_c[1]]-1):
            nn_sel[ie] = icell
            ie+=1
    nn_sel=nn_sel.astype(int)
    #%
    """Calculate the 4 triangle-areas and from it the 3 weighting factors"""
    ABC = Polygon(list(zip(gridinfo.clon[nn_sel[:]],gridinfo.clat[nn_sel[:]])))
    
    #ABP
    pollist = list(zip(gridinfo.clon[nn_sel[1:]],gridinfo.clat[nn_sel[1:]]))
    pollist.append([Ploc.x,Ploc.y])
    W1POL = Polygon(pollist)
    #BCP
    pollist = list(zip(gridinfo.clon[nn_sel[:][::2]],gridinfo.clat[nn_sel[:][::2]]))
    pollist.append([Ploc.x,Ploc.y])
    W2POL = Polygon(pollist)
    #CAP
    pollist = list(zip(gridinfo.clon[nn_sel[:2]],gridinfo.clat[nn_sel[:2]]))
    pollist.append([Ploc.x,Ploc.y])
    W3POL = Polygon(pollist)
    #%
    """Calculate weighting factors"""
    u[0] = W1POL.area/ABC.area
    u[1] = W2POL.area/ABC.area
    u[2] = W3POL.area/ABC.area
    #%

    return vert_scaling_fact, idx_below, idx_above, nn_sel[:], u


#%%
    
"""Paths, Variables and Constants"""

ICON_grid = '/scratch/snx3000/msteiner/synthetic_data/input/grid/dyn_grid.nc'
stationdir = '/store/empa/em05/msteiner/CTDAS_data/observations/input_2021'
fname_base = 'ICON-ART-OEM-UNSTR'
DATA_path = '/scratch/snx3000/msteiner/synthetic_data/output_synth_2018010100'
starttime ='2018010200'
enddtime = '2018011200'
nlev = 60 #number of vertical levels

#%%

"""Variables to extract"""

meta = { 
        'CH4_A': {'offset': 1.2e-06},
        'CH4_BG': {'offset': 0.},
        'u': {'offset': 0.},
        'v': {'offset': 0.},
        'temp': {'offset': 0.},
        'CH4_A-ENS': {'offset': 1.2e-06, 'ensemble': True},
        }

n_member = 11

#%%

"""Get grid data"""

fh_grid = Dataset(ICON_grid,'r')
class gridinfo:
    clon_vertices = np.array(fh_grid.variables['clon_vertices'])*180/np.pi
    clat_vertices = np.array(fh_grid.variables['clat_vertices'])*180/np.pi
    cells_of_vertex = np.array(fh_grid.variables['cells_of_vertex'])
    vertex_of_cell = np.array(fh_grid.variables['vertex_of_cell'])
    neighbor_cell_index = np.array(fh_grid.variables['neighbor_cell_index'])
    vlon = np.array(fh_grid.variables['vlon'])*180/np.pi
    vlat = np.array(fh_grid.variables['vlat'])*180/np.pi
    clon = np.array(fh_grid.variables['clon'])*180/np.pi
    clat = np.array(fh_grid.variables['clat'])*180/np.pi
    ncells = len(fh_grid.dimensions['cell'])

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
intp_ICON_data_ens = np.zeros((n_ens,n_member,len(latitudes),0))

#%%

"""Loop over Data Files (=timesteps)"""

datetime_list = []
print('======================================')
date_idx = 0
while looptime <= enddate:

    intp_ICON_data_det = np.concatenate(( intp_ICON_data_det,np.zeros((n_det,len(latitudes),1)) ),axis=2)
    intp_ICON_data_ens = np.concatenate(( intp_ICON_data_ens,np.zeros((n_ens,n_member,len(latitudes),1)) ),axis=3)
            
    timestring = datetime.datetime.strftime(looptime,'%Y-%m-%dT%H')
    datetime_list.append(timestring)
    DATA_file = os.path.join(DATA_path,'%s_%s:00:00.000.nc' %(fname_base,timestring))

    print('extracting from %s'%(DATA_file), flush=True)

    
    fh_data = Dataset(DATA_file,'r')
    
    class datainfo:
        z_mc = np.array(fh_data.variables['z_mc'])
        z_ifc = np.array(fh_data.variables['z_ifc'])
        
        
    ICON_data_det = np.zeros(( n_det, nlev, gridinfo.ncells ))
    ICON_data_ens = np.zeros(( n_ens, n_member, nlev, gridinfo.ncells ))
    ivar = 0
    for var in meta.keys():
        if not 'ensemble' in meta[var]:
            ICON_data_det[ivar,...] = np.array(fh_data.variables[var]) - meta[var]['offset']
            ivar+=1
    ivar = 0
    for var in meta.keys():
        if 'ensemble' in meta[var]:
            for iens in np.arange(n_member):
                varnc = var.split('-')[0]+'-%.2i'%(iens+1)
                ICON_data_ens[ivar,iens,...] = np.array(fh_data.variables[varnc]) - meta[var]['offset']
            ivar+=1
#%%
        
    """Since the stations don't walk around, I only call the function at the first timestep"""  
    
    if looptime==startdate:
        args = zip( np.arange(len(latitudes)),repeat(gridinfo),repeat(datainfo),repeat(latitudes),repeat(longitudes),repeat(asl) )
        with Pool(48) as pool:
            vsf, idxb, idxa, neighbours, u_ret = list(zip(*pool.map(intp_icon_data, args)))

        
    vsf = np.array(vsf)
    idxb = np.array(idxb, dtype=np.int)
    idxa = np.array(idxa, dtype=np.int)
    neighbours = np.array(neighbours, dtype=np.int)
    u_ret = np.array(u_ret)

    #Do the interpolation 
    for iloc in np.arange(len(latitudes)):

        ###First, the deterministic values:
        ##First, the vertical interpolation:
        vert_intp_data = ICON_data_det[:,idxb[iloc],neighbours[iloc,:]] + vsf[iloc]*(ICON_data_det[:,idxa[iloc],neighbours[iloc,:]] \
                                  -ICON_data_det[:,idxb[iloc],neighbours[iloc,:]])
        ##Now the horizontal interpolation:
        intp_ICON_data_det[:,iloc,date_idx] = u_ret[iloc,0]*vert_intp_data[:,0]+u_ret[iloc,1]*vert_intp_data[:,1] \
                            +u_ret[iloc,2]*vert_intp_data[:,2]
                           
        ###Second, the ensemble values:
        vert_intp_data = ICON_data_ens[:,:,idxb[iloc],neighbours[iloc,:]] + vsf[iloc]*(ICON_data_ens[:,:,idxa[iloc],neighbours[iloc,:]] \
                                  -ICON_data_ens[:,:,idxb[iloc],neighbours[iloc,:]])
        intp_ICON_data_ens[:,:,iloc,date_idx] = u_ret[iloc,0]*vert_intp_data[:,:,0]+u_ret[iloc,1]*vert_intp_data[:,:,1] \
                            +u_ret[iloc,2]*vert_intp_data[:,:,2]
#%%
    """Update time"""
    looptime += delta
    date_idx += 1
#%%
"""Save as netcdf"""
with Dataset('output/extracted_%s.nc'%(DATA_path.split('/')[-1],), mode='w') as ofile:

    oens = ofile.createDimension('ens', n_member)
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
            varnc = var.split('-')[0]+'_ENS'
            ovar = ofile.createVariable(varnc, np.float32, ('ens','sites','time'))
            ovar[:,:,:] = intp_ICON_data_ens[ivar,:,:,:]
            ivar+=1

    
    oname[:] = np.array(stationnames[:])
    otimes[:] = np.array(datetime_list)
