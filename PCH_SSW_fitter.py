'''
This code takes .npy files (output from SSW_finder.py). It then creates a time-height slice of Polar Cap Height (PCH)anomalies for each event and combines these into separate netCDF files for splits, displacements and unclassifiable events. Slices are taken from 60 days before to 60 days after SSW onset.
PCH anomalies are preprocesssed for the Polar Cap (65-90N) from gepotentail height at all model levels, using Climate Data Operators
'''

import numpy as np
import xarray as xr


#load .npy files of SSW dates
split_dates=np.load('model_split_dSSW_dates.npy',allow_pickle=True)
disp_dates=np.load('model_disp_dSSW_dates.npy',allow_pickle=True)
un_dates=np.load('model_unclass_dates_all.npy',allow_pickle=True)

# open  PCH standardised anomaly file using xarray, extract variable
filename = ('PCH_SA.nc')
DS = xr.open_dataset(filename,use_cftime=True)
PCH =DS['zg']
DS.close()

#calculate time-height slices, and mean time-height slice for each event type, save to netCDF
events=[]
for iday in range(len(split_dates)):
    IDX = np.where(PCH.time==split_dates[iday])[0][0]
    event = PCH[IDX-60:IDX+61,:,:,:]
    event = event.assign_coords(time=range(-60,61))
    events.append(event)

SSW_split_slices=xr.concat(events,dim='events')
SSW_split_mean=xr.DataArray.mean(SSW_split_slices,dim='events',keep_attrs=True)

SSW_split_slices.to_netcdf('model_PCH_split_events.nc')
SSW_split_mean.to_netcdf('model_PCH_split_mean.nc')

events=[]
for iday in range(len(disp_dates)):
    IDX = np.where(PCH.time==disp_dates[iday])[0][0]
    event = PCH[IDX-60:IDX+61,:,:,:]
    event = event.assign_coords(time=range(-60,61))
    events.append(event)
    
SSW_disp_slices=xr.concat(events,dim='events')
SSW_disp_mean=xr.DataArray.mean(SSW_disp_slices,dim='events',keep_attrs=True)
SSW_disp_slices.to_netcdf('model_PCH_disp_events.nc')
SSW_disp_mean.to_netcdf('model_PCH_disp_mean.nc')

events=[]
for iday in range(len(un_dates)):
    IDX = np.where(PCH.time==un_dates[iday])[0][0]
    event = PCH[IDX-60:IDX+61,:,:,:]
    event = event.assign_coords(time=range(-60,61))
    events.append(event)
    
SSW_un_slices=xr.concat(events,dim='events')
SSW_un_mean=xr.DataArray.mean(SSW_un_slices,dim='events',keep_attrs=True)
SSW_un_slices.to_netcdf('model_PCH_un_events_all.nc')
SSW_un_mean.to_netcdf('model_PCH_un_mean_all.nc')
