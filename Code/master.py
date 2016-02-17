#------------------------------------------------------------------------------
# Name:        master.py
# Purpose:     Run MCMC calibration of WBM for multiple sites.
#
# Author:      James Sample
#
# Created:     27/12/2014
# Copyright:   (c) James Sample and JHI, 2014
#------------------------------------------------------------------------------
""" This script uses the subprocess module to run MCMC-based autocalibration of
    the water balance model for multiple sites. The script runs one site at a 
	time, and the MCMC runs for each site are distributed over the specified 
	number of processors. The use of subprocess ensures that all threads are 
	closed properly before processing for the next site begins.
	
	Note that this script will take a long time to run for many sites - expect
	days to weeks on virtual machine 'arcgis1'.
"""
import subprocess as sp, sys

# #############################################################################
# User input
in_h5_path = r'Y:\Water_Balance_Modelling_2014\WBM_2014_Monthly_Input_File.h5'
obs_flow_path = r'Y:\Water_Balance_Modelling_2014\Observed_Flows\Obs_Monthly_Flows_1971_2010.csv'
stn_path = r'Y:\Water_Balance_Modelling_2014\Observed_Flows\Observed_Flows_Sites.csv'
mask_path = r'Y:\Water_Balance_Modelling_2014\ASCII_Grids\all_obs_catchments_1km.txt'
slave_path = r'Y:\Water_Balance_Modelling_2014\Python\Monthly_WBM_2014\mcmc_slave.py'
out_fold = r'Y:\Water_Balance_Modelling_2014\Calibration_Plots\Test_Delete'

# Site of interest
sites_list = ['7002', '11001', '14001', '15006', '19001', '21009', '77002',
              '81002', '84005', '94001']

# Period of interest for modelling
st_yr, end_yr = 1961, 1990

# Period of interest for calibration. This is different as you need to run the
# model for a few years first to allow it to equilibrate
cal_st_yr, cal_end_yr = 1971, 1990

# Models of interest
models = 'base'

# Select land use grid for PET to AET factors. Choose from:
# lcms_1988, lcm_2007 or lu_2050
lu_grid_path = r'lcms_1988'

# MCMC params
ndim, nwalkers = 4, 50
nsteps, nburn = 600, 300
n_proc = 4 # Number of processors available

# Define limits for uniform priors. Walkers will start at random points
# between min and max values entered here
k_s_min, k_s_max = 0.1, 2.			# Soil water residence time
k_g_min, k_g_max = 0.005, 0.04      # Groundwater residence time
et_min, et_max = 0.25, 1.           # ET correction factor 
sigma_min, sigma_max = 0, 10.       # Likelihood sigma for basic iid Gaussian errors

# Apply square root transformation to mod and obs values before evaluating
# likelihood
sqrt_trans = True
# #############################################################################

# Loop over sites
for site_code in sites_list:
    print 'Processing: %s.' % site_code
	
	# Call main MCMC script
    arg_list = [sys.executable, slave_path, in_h5_path, obs_flow_path, stn_path, 
                mask_path, out_fold, site_code, str(st_yr), str(end_yr), str(cal_st_yr), 
                str(cal_end_yr), models, lu_grid_path, str(ndim), str(nwalkers),
                str(nsteps), str(nburn), str(n_proc), str(k_s_min),
                str(k_s_max), str(k_g_min), str(k_g_max), str(et_min), 
                str(et_max), str(sigma_min), str(sigma_max), str(sqrt_trans)]
    sp.call(arg_list, shell=True)