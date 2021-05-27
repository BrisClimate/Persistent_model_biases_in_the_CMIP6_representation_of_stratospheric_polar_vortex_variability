# Persistent_model_biases_in_the_CMIP6_representation_of_stratospheric_polar_vortex_variability
Code to support the paper

SSW_finder.py. This code identifies the onset dates of SSWs.It also classifies each event as split, displacement or unclassifiable using the vortex moments approach

vor_fast.py and vor_fast_Setup.py are provided from William Seviour's Github page (wseviour) and are installed by SSW_finder.py to calculate the vortex moments. If you use the vortex moments code please cite:
Seviour, W. J. M., D. M. Mitchell, and L. J. Gray (2013), A practical method to identify displaced and split stratospheric polar vortex events, Geophys. Res. Lett., 40, 5268-5273 doi: 10.1002/grl.50927

PCH_SSW_fitter.py. This code calculates time-height slices of SSWs, using .npy files from SSW_finder.py

moments_mode_calculator.py. This codecalculates the aspect ratio and centroid latitude modes for a model.
