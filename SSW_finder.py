'''
    this code accompnies the paper "Persistent biases in the CMIP6 representation of stratospheric polar vortex variability" Journal of Geophysical Research: Atmospheres.
    
    This code identifies SSW onset dates, and classifies each SSW event. Dates for each event type are saved as .npy files
    
    The code reads in files for each winter that are preprocessed using Climate Data Operators, for example to select level, area, regridding, zonal means etc.
    
'''

import xarray as xr
import numpy as np
import more_itertools as mit
import glob
import datetime as dt 
import vor_fast
import vor_fast_setup
import pandas as pd


# filepaths to obtain 10hPa 60N zonal mean uwind values from separate winter files
filepaths1=[]
for filepath1 in glob.iglob('ua10/NDJFMA/winter*'):
    filepaths1.append(filepath1)
filepaths1.sort()

# filepaths to obtain 10hPa GPH for northern hemisphere
filepaths2=[]
for filepath2 in glob.iglob('zg10/NDJFMA/winter*'):
    filepaths2.append(filepath2)
filepaths2.sort()

# open data files, extract variable, lon and lat
vortex_class=[]
SSW_dates=[]
for i in range(0,len(filepaths1)):
    filename1 = (filepaths1[i])
    DS = xr.open_dataset(filename1,use_cftime=True)
    u10=DS['ua']
    u10 = u10.squeeze()
    DS.close()
    
    filename2 = (filepaths2[i])
    DS = xr.open_dataset(filename2,use_cftime=True)
    gph=DS['zg']
    gph = gph.squeeze()
    lons=DS.lon
    lats=DS.lat
    days=DS.time
    DS.close()
    
# Set up cartesian mapping xypoints and restrict to NH
    gph_nh, lats_nh, xypoints = vor_fast_setup.setup(gph.values,lats.values,lons.values,'NH')
    
# Set up moment diagnostics
    aspect = np.empty(0)
    latcent = np.empty(0)
    
# find indices of dates of wind reversal
    idx=np.where(u10<0)

#group the indices into events
    events=[list(group) for group in mit.consecutive_groups(idx[0])]
    print(events)
    if not (events):
        print ('no events')
        continue

# check last event is not a final warming. If it is, remove
    if events[-1][-1] >len(u10)-10:
        del events[-1]
    if not (events):
        print ('no events')
        continue
# check that there are at least 10 consecutive days of westerlies at some point after last remaining event
    if len(events)==1:
        idx_end=np.where(u10[events[0][-1]+1:-1]>0)
        events_end=[list(group) for group in mit.consecutive_groups(idx_end[0])]
        if len(events_end[0])<10:
            del events[0]
    else:
        for i in range(len(events)-1,0,-1):
            idx_end=np.where(u10[events[i][-1]+1:-1]>0)
            events_end=[list(group) for group in mit.consecutive_groups(idx_end[0])]
            consec_pos=[]
            for k in range (0,len(events_end)):
                length=len(events_end[k])
                consec_pos.append(length)
            if all (ele <10 for ele in consec_pos):
                del events[i]

# delete any events starting in April or November
    months=[]
    for e in range(0,len(events)):
        mon=int(DS.time[events[e][0]].dt.strftime('%m'))
        months.append(mon)

    remove=[]
    for e in range(0,len(events)):
        if months[e]==4 or months[e]==11:
            remove.append(e)

    for index in sorted(remove,reverse=True):
        del events[index]
    print(events)

    if not (events):
        print ('no events')

#for  events, check if 20 consecutive days of westerlies before events
    remove=[]
    for j in np.arange(len(events)):
        if np.any(u10[events[j][0]-20:events[j][0]]<0):
            remove.append(j)
    print('to remove= ', remove)

    for index in sorted(remove,reverse=True):
        del events[index]
    print(events)
        
    if not (events):
        print ('no events')

#final check for last event, does it have 10 consecutive days of westerlies afterwards
    if len(events)>0:
        idx_end=np.where(u10[events[-1][-1]+1:-1]>0)
        events_end=[list(group) for group in mit.consecutive_groups(idx_end[0])]
        consec_pos=[]
        for k in range (0,len(events_end)):
            length=len(events_end[k])
            consec_pos.append(length)
        if all (ele <10 for ele in consec_pos):
            del events[-1]

    elif not (events):
        print ('no events')

#extract the onset dates from the events
    event_length=[]
    start_idx=[]
    for i in range(0,len(events)):
        start=events[i][0]
        length=len(events[i])
        start_idx.append(start)
        event_length.append(length)

    print('indices of onset dates=',start_idx)
    print('event lengths=',event_length)

    onset_dates=DS.time[start_idx].values
    print('onset dates =',onset_dates)

    #using vortex moments classify vortex as split, displacement or unclassifiable
    #NB calculate vorTex edge for each model separately and add below
    vortex_ID=[]
    for idx in range(0,len(start_idx)):
        date_ind=start_idx[idx]
    
        for iday in range(date_ind-10,date_ind+11):
            print('Calculating moments for day '+str(iday))
            moments = vor_fast.calc_moments(gph_nh[iday,:,:],lats_nh,lons.values,xypoints,
                                        hemisphere='NH',field_type='GPH',
                                        edge=3.019e4,resolution='low')
            aspect = np.append(aspect, moments['aspect_ratio'])
            latcent = np.append(latcent, moments['centroid_latitude'])
                                        
        print(aspect)
        print(latcent)
        aspect_days=np.count_nonzero(aspect>2.4)
        latcent_days=np.count_nonzero(latcent<66)
        print('days above aspect threshold=',aspect_days)
        print('days below latcent threshold=',latcent_days)
        vortex_type=[]
        if aspect_days > latcent_days:
            vortex_type ="split"
        elif aspect_days < latcent_days:
            vortex_type = "displaced"
        else:
            vortex_type = "unclassified"
        vortex_ID.append(vortex_type)

    SSW_dates.extend(onset_dates)
    vortex_class.extend(vortex_ID)

print('SSW dates are:',SSW_dates)
print('vortex types:', vortex_class)
print(len(SSW_dates))

# make dataframe of results for model
model_SSWs=pd.DataFrame()
model_SSWs['onset dates']=SSW_dates
model_SSWs['vortex type']=vortex_class

split_model=model_SSWs[(model_SSWs['vortex type']=='split')]
disp_model=model_SSWs[(model_SSWs['vortex type']=='displaced')]
unclass_model=model_SSWs[(model_SSWs['vortex type']=='unclassified')]
split_dates=np.array(split_model['onset dates'])
disp_dates=np.array(disp_model['onset dates'])
unclass_dates=np.array(unclass_model['onset dates'])

#save SSW dates as .npy files for use in subsequent code.
np.save('model_split_dates_all.npy',split_dates)
np.save('model_disp_dates_all.npy',disp_dates)
np.save('model_unclass_dates_all.npy',unclass_dates)
np.save('model_all_SSW_dates.npy',SSW_dates)

print('number of displacements: ',len(disp_dates))
print('number of splits: ',len(split_dates))
print('number of unclassified: ',len(unclass_dates))
