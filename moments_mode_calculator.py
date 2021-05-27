"""
This example uses vor_fast to calculate daily moment diagnostics
for DJFM for a model. It also calculates the modal centroid latitude and aspect ratio for DJFM for the model.

"""

from netCDF4 import Dataset
import numpy as np
import vor_fast
import vor_fast_setup
from scipy import stats

# Read in NetCDF file with geopotential height values 10hPa northern hemisphere
ncin = Dataset('DJFM_data.nc', 'r')
gph = ncin.variables['zg'][:].squeeze()
lons = ncin.variables['lon'][:]
lats = ncin.variables['lat'][:]
days = ncin.variables['time'][:]
ncin.close()

# Set up cartesian mapping xypoints and restrict to NH
gph_nh, lats_nh, xypoints = vor_fast_setup.setup(gph,lats,lons,'NH')

# Set up moment diagnostics
aspect = np.empty(0)
latcent = np.empty(0)


# Calculate diagnostics for each day
# vortex edge should be calculated for each individual model
for iday in range(len(days)):
    print('Calculating moments for day '+str(iday))
    moments = vor_fast.calc_moments(gph_nh[iday,:,:],lats_nh,lons,xypoints,
                                    hemisphere='NH',field_type='GPH',
                                    edge=3.016e4,resolution='low')
    aspect = np.append(aspect, moments['aspect_ratio'])
    latcent = np.append(latcent, moments['centroid_latitude'])

np.save('model_aspect.npy',aspect)
np.save('model_centlat.npy',latcent)

#remove nans and inf values from arrays and limit aspect values to <100
aspect=aspect[np.logical_not(np.isnan(aspect))]
aspect=aspect[aspect < float('+inf')]
aspect=aspect[aspect < 100.0]
latcent=latcent[np.logical_not(np.isnan(latcent))]
latcent=latcent[latcent < float('+inf')]

#for latcent, cube the values to transform the distribution
latcent3=np.power(latcent,3)
#fit Gaussian distribution and use KS test for fit
dist=getattr(stats,'norm')
parameters=dist.fit(latcent3)
(mean,SD)=parameters
ks_stat,ks_pval=stats.kstest(latcent3,"norm",parameters)
# convert latitude back (cube root)
mode=np.power(mean,1/3)

print('latcent mode:',mode)
print('norm KS statistic and p value:',ks_stat,ks_pval)

#fit GEV to aspect data
dist2=getattr(stats,'genextreme')
parameters=dist2.fit(aspect)
(shape,location,scale)=parameters
ks_stat_GEV,ks_pval_GEV = stats.kstest(aspect,"genextreme",parameters)

print('location:', location)
print(' GEV KS statistic and p value:',ks_stat_GEV,ks_pval_GEV)

#aspect mode calculated according to Seviour et al. 2016 (S16)
asp_params=stats.genextreme.fit(aspect)
asp_pdf=stats.genextreme.pdf(aspect,asp_params[0],loc=asp_params[1],scale=asp_params[2])
index=np.argmax(asp_pdf)
max_asp_S16=aspect[index]

print ('S16 aspect:',max_asp_S16)

#kde alternative
# aspect
min=1
max=np.round(np.max(aspect))+1
divs=(max-min)*100

asp_xs=np.linspace(min,max,divs)
asp_density=stats.gaussian_kde(aspect)
asp_ys=asp_density(asp_xs)
index=np.argmax(asp_ys)
max_asp_kde=asp_xs[index]

#latcent
min1=np.round(np.min(latcent))-1
max1=90
divs1=(max1-min1)*100

lat_xs=np.linspace(min1,max1,divs1)
lat_density=stats.gaussian_kde(latcent)
lat_ys=lat_density(lat_xs)
index1=np.argmax(lat_ys)
max_lat_kde=lat_xs[index1]

print ('kde aspect:',max_asp_kde)
print ('kde latcent:',max_lat_kde)






