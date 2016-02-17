#------------------------------------------------------------------------------
# Name:        mcmc_slave.py
# Purpose:     MCMC calibration script for the new water balance model.
#
# Author:      James Sample
#
# Created:     09/12/2014
# Copyright:   (c) James Sample and JHI, 2014
#------------------------------------------------------------------------------
""" Uses the Emcee module to calibrate k_s and k_g. This version assumes the 
    parameters are are spatially invariant and sets priors that are uniform
    between specified max and min values.
    
    Also includes a correction factor for pet as an additional parameter. It
    is clear from manual experimentation that the total volume of water entering
    the system according to my weather data is lower than the total volume
    observed in rivers. Either my rainfall is too low or ET is too high. This
    code introduces a correction factor for ET.
"""
import h5py, numpy as np, calendar, pandas as pd, datetime as dt, emcee
import input_output as io, drainage as dr, corner, bottleneck as bn
import matplotlib.pyplot as plt, sys, os

np.seterr(invalid='ignore')
np.set_printoptions(precision=2)

def log_prior(p):
    """ Calculate log prior probability.
	
	Args:
		p: Vector of parameters [k_s, k_g, et, sigma].
	
	Returns:
		0 if parameter set is within prior range; -inf otherwise.
    """
    # Unpack parameters
    k_s, k_g, et, sigma = p
    
    # Assume uniform priors 
    if (k_s<=k_s_min) or (k_s>=k_s_max):
        return -1.*np.inf 
    elif (k_g<=k_g_min) or (k_g>=k_g_max):
        return -1.*np.inf
    elif (et<=et_min) or (et>=et_max):
        return -1.*np.inf
    elif (sigma<=sigma_min) or (sigma>=sigma_max):
        return -1.*np.inf 
    elif (k_s==k_g):
        return -1.*np.inf # k_s cannot equal k_g        
    else:    
        return 0 # The zero gives us a flat prior

def log_likelihood(p, obs):
    """ Calculate log likelihood.

	Args:
		p: Vector of parameters [k_s, k_g, et, sigma].
	
	Returns:
		Log probability of the data given the parameters assuming
		iid Gaussian error structure.
    """
    # Unpack parameters
    k_s, k_g, et, sigma = p

    # Assume uniform priors 
    if (k_s<=k_s_min) or (k_s>=k_s_max):
        return -1.*np.inf 
    elif (k_g<=k_g_min) or (k_g>=k_g_max):
        return -1.*np.inf
    elif (et<=et_min) or (et>=et_max):
        return -1.*np.inf
    elif (sigma<=sigma_min) or (sigma>=sigma_max):
        return -1.*np.inf 
    elif (k_s==k_g):
        return -1.*np.inf # k_s cannot equal k_g         
    else:           
        # Estimate model results
        mod = run_wbm(k_s, k_g, et)
        
        # Truncate to calibration period
        mod = mod.truncate(before='%s-01-01' % cal_st_yr, 
                           after='%s-12-31' % cal_end_yr)
        
        # Convert to array
        mod = mod['Runoff_m3/s'].values
        
        assert obs.shape==mod.shape, 'Observed and modelled arrays are not same shape.'

        # Apply data transformation if specified
        if (sqrt_trans==True):
            mod = mod**0.5
            obs = obs**0.5
            
        # Calculate likelihood    
        lp = -0.5*np.sum(np.log(2.*np.pi*sigma**2) + (((mod - obs)**2)/sigma**2))
        
        return lp
  
def log_posterior(p, obs):
    """ Calculate log posterior (ignore normalising constant).

	Args:
		p:   Vector of parameters [k_s, k_g, et, sigma].
		obs: Vector of runoff observations.
		
	Returns:
		Log posterior probability.
    """
    return log_prior(p) + log_likelihood(p, obs)
    
def check_params():
    """ Validates user input.
    """
    # Check that co-ords lie in Scotland
    if (xmin < 0) or (xmax > 485000) or (ymin < 520000) or (ymax > 1235000):
        raise ValueError('The bounding box supplied does not lie '
                         'within Scotland.')

    # Check that the co-ordinates are an even multiple of 1000
    if ((xmin % 1000 != 0) or (xmax % 1000 != 0) or
        (ymin % 1000 != 0) or (ymax % 1000 != 0)):
        raise ValueError('Please enter co-ordinates that are an integer '
                         'multiple of 1000m')
               
    # Check years
    assert st_yr >= 1961, 'Start year cannot be before 1961.'
    assert end_yr <= 2090, 'End year cannot be after 2090.'
    
    # Check models are OK
    all_models = ['base', 'snow_corr_base', 'afixa', 'afixc', 'afixl', 'afixm',
                  'afixo', 'afixh', 'afixi', 'afixj', 'afixk', 'afgcx', 'afixq']
    for model in models:
        assert model in all_models, '%s is not a valid climate model.' % model

    # If the run includes the baseline, check that period of interest is 
    # compatible
    if ('base' in models) or ('snow_corr_base' in models):
        assert (st_yr>=1961) and (end_yr<=2010), ('Period selected is not '\
                                                  'compatible with baseline '\
                                                  'model runs.')
        
def get_grid_indices(xmin, xmax, ymin, ymax):
    """ Take the OSGB-1936 co-ordinates of a bounding rectangle and
        determine the array indices required to define the desired area.
	
	Args:
		xmin:  OSGB-1936 metres for minimum Easting.
		xmax:  OSGB-1936 metres for maximum Easting.
		ymin:  OSGB-1936 metres for minimum Northing.
		ymax:  OSGB-1936 metres for maximum Northing.
	
	Returns:
		Tuple of array indices (xmin_idx, xmax_idx, ymin_idx, ymax_idx).
    """
    # Calculate indices
    xmin_idx = xmin / 1000
    xmax_idx = xmax / 1000
    ymin_idx = (1235000-ymax) / 1000
    ymax_idx = (1235000-ymin) / 1000

    return (xmin_idx, xmax_idx, ymin_idx, ymax_idx)
   
def run_wbm(k_s, k_g, et):
    """ Run the WBM with the specified parameters.

	Args:
		k_s: Soil residence time (days).
		k_g: Groundwater residence time (days).
		et:  Correction factor applied to the ET data (dimensionless).
	
	Returns:
		Pandas dataframe of monthly river flows.
    """   
	
    # Define dicts storing number of days in each period. One dict for leap years
    # the other for non-leap years
    days_in_month_dict = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30,
                          10:31, 11:30, 12:31}
    days_in_month_lpyr_dict = {1:31, 2:29, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 
                               9:30, 10:31, 11:30, 12:31}
                                   
    # Validate user input
    check_params()
    
    # Get array indices for bounding box
    xmin_idx, xmax_idx, ymin_idx, ymax_idx = get_grid_indices(xmin, xmax, 
                                                              ymin, ymax)
                                                                 
    # Open H5 file
    h5 = h5py.File(in_h5_path, 'r')
    
    # Get soil properties
    fc, sat, bfi = io.read_soil_properties(h5, xmin_idx, xmax_idx, ymin_idx, 
                                           ymax_idx)

    # Get LU grid
    lu_grid = io.read_land_use(h5, lu_grid_path, xmin_idx, xmax_idx, ymin_idx, 
                               ymax_idx)
                           
    # Get soil and groundwater time constant grids
    k_s = np.ones(shape=bfi.shape)*k_s
    k_g = np.ones(shape=bfi.shape)*k_g
    
    # Water between fc and sat
    phi = sat - fc
    
    # Loop over models
    for model in models:  
        # Initial water level mid-way between fc and sat
        surf_lev = (fc + sat)/2.
        gw_lev = np.zeros(shape=phi.shape)
        
        # Dict to store data
        data_dict = {'Date':[],
                     'of':[],
                     'ssf':[],
                     'gwf':[],
                     'Runoff_m3/s':[]}
                     
        # Loop over years
        for year in range(st_yr, end_yr+1):
            for month in range(1, 13):                
                # Get date
                date = dt.date(year, month, 1)
                
                # Get number of days in month, allowing for leap years
                if calendar.isleap(year):
                    days_in_month = days_in_month_lpyr_dict[month]
                else:
                    days_in_month = days_in_month_dict[month]
                           
                # Get met data
                pptn, pet = io.read_met_data(h5, model, xmin_idx, xmax_idx,
                                             ymin_idx, ymax_idx, year, month)

                # Correct grass reference PET for land use and apply ET
                # correction factor
                pet = pet*lu_grid*et
                
                # Convert from mm/month to mm/day
                pptn = pptn/days_in_month
                pet = pet/days_in_month
                
                # Calculate HER
                her = pptn - pet # mm/day
    
                # Get drainage params for this time step
                drainage_params = dr.calculate_drainage(her, surf_lev, gw_lev, 
                                                        fc, days_in_month, phi,
                                                        k_s, k_g, bfi)
                
                # NB: the net_her value returned here is for the whole month
                # i.e. mm/month NOT mm/day like her, above
                net_her, of, ssf, gwf, surf_lev, gw_lev = drainage_params
                
                # Calculate monthly runoff
                ro = of + ssf + gwf
                
                # Apply mask
                ro[mask!=site_code] = np.nan
                of[mask!=site_code] = np.nan
                ssf[mask!=site_code] = np.nan
                gwf[mask!=site_code] = np.nan
                
                # Calculate monthly total in m3/s
                ro_mon = 1000.*bn.nansum(ro)/(days_in_month*24*60*60)
                of_mon = 1000.*bn.nansum(of)/(days_in_month*24*60*60)
                ssf_mon = 1000.*bn.nansum(ssf)/(days_in_month*24*60*60)
                gwf_mon = 1000.*bn.nansum(gwf)/(days_in_month*24*60*60)
                
                # Append to data dict
                data_dict['Date'].append(date)
                data_dict['of'].append(of_mon)
                data_dict['ssf'].append(ssf_mon)
                data_dict['gwf'].append(gwf_mon)
                data_dict['Runoff_m3/s'].append(ro_mon)
                
    # Close file
    h5.close()
    
    # Build df
    df = pd.DataFrame(data_dict)
    df.index = pd.to_datetime(df['Date'])
    del df['Date']
    
    return df
    
def main():
    """ The main function for the MCMC analysis.
    """         
    # Generate starting guesses
    k_s_guess = np.random.uniform(low=k_s_min, high=k_s_max, size=nwalkers)
    k_g_guess = np.random.uniform(low=k_g_min, high=k_g_max, size=nwalkers)
    et_guess = np.random.uniform(low=et_min, high=et_max, size=nwalkers)
    sigma_guess = np.random.uniform(low=sigma_min, high=sigma_max, size=nwalkers)
    starting_guesses = np.column_stack((k_s_guess, k_g_guess, 
                                        et_guess, sigma_guess))
    
    # Prepare to sample. The params, p, are automatically passed to log_posterior
    # as part of ndim. "args" lists the other params that are also necessary
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, 
                                    threads=n_proc, args=[obs,])
    
    # Run sampler
    sampler.run_mcmc(starting_guesses, nsteps)

    # Print some stats. based on run properties
    print '\n'
    print 'Acceptance fraction: ', sampler.acceptance_fraction
    print 'Autocorrelation time: ', sampler.acor
    
    # Get results
    # Plot traces, including burn-in
    fig, axes = plt.subplots(nrows=4, ncols=1)
    axes[0].plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
    axes[0].set_title('k_s')
    
    axes[1].plot(sampler.chain[:,:,1].T, '-', color='k', alpha=0.3)
    axes[1].set_title('k_g')

    axes[2].plot(sampler.chain[:,:,2].T, '-', color='k', alpha=0.3)
    axes[2].set_title('et')
    
    axes[3].plot(sampler.chain[:,:,3].T, '-', color='k', alpha=0.3)
    axes[3].set_title('sigma')
    
    plt.subplots_adjust(hspace=0.5)    
    
    # Save fig 
    out_path = os.path.join(out_fold, '%s_Traces.png' % site_code)
    plt.savefig(out_path)
    plt.close()
    
    # Discard burn-in
    sample = sampler.chain[:, nburn:, :].reshape((-1, ndim))
    
    # Print results
    print '\n'
    print 'Mean k_s: %.2f.' % sample[:, 0].mean()
    print 'Mean k_g %.2f.' % sample[:, 1].mean()
    print 'Mean et %.2f.' % sample[:, 2].mean()
    print 'Mean sigma %.2f.' % sample[:, 3].mean()
    
    # Plot using corner
    tri = corner.corner(sample,
                          labels=['k_s', 'k_g', 'et', 'sigma'], 
                          quantiles=[0.10, 0.50, 0.90])

    # Save fig 
    out_path = os.path.join(out_fold, '%s_Triangle_Plot.png' % site_code)
    tri.savefig(out_path)
    plt.close()
    
    # Save the ENTIRE chain (inc. burn-in period) for later analysis
    out_path = os.path.join(out_fold, '%s_Samples.npy' % site_code)
    np.save(out_path, sampler.chain.reshape((-1, ndim)))

# #############################################################################
# Handle input from master.py
in_h5_path = sys.argv[1]
obs_flow_path = sys.argv[2]
stn_path = sys.argv[3]
mask_path = sys.argv[4]
out_fold = sys.argv[5]

# Site of interest
site_code = int(sys.argv[6])

# Period of interest for modelling
st_yr, end_yr = int(sys.argv[7]), int(sys.argv[8])

# Period of interest for calibration. This is different as you need to run the
# model for a few years first to allow it to equilibrate
cal_st_yr, cal_end_yr = int(sys.argv[9]), int(sys.argv[10])

# Models of interest
models = [sys.argv[11], ]

# Select land use grid for PET to AET factors. Choose from:
# lcms_1988, lcm_2007 or lu_2050
lu_grid_path = sys.argv[12]

# MCMC params
ndim, nwalkers = int(sys.argv[13]), int(sys.argv[14])
nsteps, nburn = int(sys.argv[15]), int(sys.argv[16])
n_proc = int(sys.argv[17]) # Number of processors available

# Parameter guesses for starting walkers. Walkers will start at random points
# between min and max values entered here
k_s_min, k_s_max = float(sys.argv[18]), float(sys.argv[19])
k_g_min, k_g_max = float(sys.argv[20]), float(sys.argv[21])
et_min, et_max = float(sys.argv[22]), float(sys.argv[23])
sigma_min, sigma_max = float(sys.argv[24]), float(sys.argv[25])

# Apply square root transformation to mod and obs values before evaluating
# likelihood
sqrt_trans = bool(sys.argv[26])
# #############################################################################

# Read the station data
stn_df = pd.read_csv(stn_path, index_col=1)

# Read the observed data
obs_df = pd.read_csv(obs_flow_path, parse_dates=True, index_col=0, 
                     dayfirst=True)
obs_df = obs_df.truncate(before='%s-01-01' % cal_st_yr, 
                         after='%s-12-31' % cal_end_yr)
obs = obs_df[str(site_code)].values

# Get spatial data for this site
xmin = stn_df.ix[int(site_code)]['xmin']
xmax = stn_df.ix[int(site_code)]['xmax']
ymin = stn_df.ix[int(site_code)]['ymin']
ymax = stn_df.ix[int(site_code)]['ymax']
rows = (ymax - ymin)/1000
cols = (xmax - xmin)/1000        

# Read mask grid
mask = io.read_ascii(mask_path, xmin=xmin, xmax=xmax, 
                     ymin=ymin, ymax=ymax, exptd_rows=rows, 
                     exptd_cols=cols)
                      
if __name__ == "__main__":
    import time
    st_time = time.time()
    main()
    end_time = time.time()
    print '\n'
    print 'Finished.'
    print 'Processing time: %.2f hours.' % ((end_time - st_time)/3600.)