#------------------------------------------------------------------------------
# Name:        drainage.py
# Purpose:     Functions for calculating drainage in the new WBM.
#
# Author:      James Sample
#
# Created:     31/07/2014
# Copyright:   (c) James Sample and JHI, 2014
#------------------------------------------------------------------------------
""" Functions for estimating drainage.
"""

def check_params(S_0, k_s, k_g, phi, b):
    """ Check drainage parameters.
    """   
    assert (S_0 <= k_s*phi).all(), 'S_0 cannot be greater than k_s*phi.'
    assert ((b>=0).all() and (b<=1).all()), 'BFI must be between 0 and 1.'
    assert (k_s != k_g).all(), 'k_s cannot equal k_g.'

def t_critical(R, S_0, phi, k_s):
    """ Calculates the time when water level equals saturation capacity.
	
	Args:
		R:   Grid of inputs (essentially HER).
		S_0: Grid of initial soil outflows.
		phi: Grid of saturation capacity.
		k_s: Soil residence time.
		
	Returns:
		Array of critical times.
    """
    import numpy as np
    
    t_c = np.log((R - S_0)/(R - k_s*phi))/k_s

    # Handle edge cases.
    # 1. If the above equn becomes log(0 / 0) it evaluates to nan, but for 
    # consistency it should give zero, as R = S_0 = k_s*phi
    bool_mask = (R == k_s*phi) & (S_0 == k_s*phi)
    t_c[bool_mask] = 0
    
    # 2. If R < k_s*phi or R = S_0, water never reaches sc, so should give +inf
    # rather than Nan or -inf. S_0 can never exceed k_s*phi, so the following
    # is sufficient
    t_c[R<k_s*phi] = np.inf
    
    return t_c

def soil_outflow_equn(t, R, S_0, k_s):
    """ Evaluates the soil outflow rate at time t.

	Args:
		t:	 Time of interest.	
		R:   Grid of inputs (essentially HER).
		S_0: Grid of initial soil outflows.
		k_s: Soil residence time.
		
	Returns:
		Array.
    """
    import numpy as np
    
    return R - (R - S_0)*np.exp(-1.*k_s*t)

def gw_outflow_equn(t, R, S_0, D_0, k_s, k_g, b):
    """ Evaluates the gw outflow rate at time t. 

	Args:
		t:	 Time of interest.	
		R:   Grid of inputs (essentially HER).
		S_0: Grid of initial soil outflows.
		D_0: Grid of initial groundwater outflows.
		k_s: Soil residence time.
		k_g: Groundwater residence time.
		b:   Base Flow Index.
		
	Returns:
		Array.
    """
    import numpy as np
    
    return (R*b*(1 - np.exp(-1.*k_g*t)) + 
            k_g*b*(R - S_0)*(np.exp(-1.*k_s*t) - np.exp(-1.*k_g*t))/(k_s - k_g) +
            D_0*np.exp(-1.*k_g*t))
            
def unsaturated_vol(R, S_0, t_1, t_c, k_s):
    """ Calculates the drainage taking place between t_1 and t_c when the water
        level is below saturation capacity.

	Args:
		R:   Grid of inputs (essentially HER).
		S_0: Grid of initial soil outflows.
		t_1: Initial time of interest (t_1 < t_c).
		t_c: Time at which water level reaches saturation capacity.
		k_s: Soil residence time.
		
	Returns:
		Array.
    """
    import numpy as np
    
    return (R*(t_c - t_1) + 
            (R - S_0)*(np.exp(-1.*k_s*t_c) - np.exp(-1.*k_s*t_1))/k_s)

def gw_vol(R, S_0, D_0, t_1, t_2, k_s, k_g, b):
    """ Evaluates the GW drainage between t_1 and t_2. The parameters R, S_0
        and D_0 can be specified explicitly here as it useful to be able to 
        change them from the global R and S_0 values used by the other 
        functions.

	Args:
		R:   Grid of inputs (essentially HER).
		S_0: Grid of initial soil outflows.
		D_0: Grid of initial groundwater outflows.
		t_1: Initial time of interest (t_1 < t_2).
		t_2: Final time of interest (t_1 < t_2).		
		k_s: Soil residence time.
		k_g: Groundwater residence time.
		b:   Base Flow Index.
		
	Returns:
		Array.
    """
    import numpy as np
    
    return (R*b*(t_2 - t_1) +
            R*b*(np.exp(-1.*k_g*t_2) - np.exp(-1.*k_g*t_1))/k_g +
            b*(R - S_0)*(np.exp(-1.*k_g*t_2) - np.exp(-1.*k_g*t_1))/(k_s - k_g) +
            k_g*b*(R - S_0)*(np.exp(-1.*k_s*t_1) - np.exp(-1.*k_s*t_2))/(k_s*(k_s - k_g)) +
            D_0*(np.exp(-1.*k_g*t_1) - np.exp(-1.*k_g*t_2))/k_g)

def saturated_vol(t_sat, phi, k_s):
    """ Calculates the drainage taking place between t_c and t_2 when the water
        level is above saturation capacity.

	Args:
		t_sat: Time for which soil is saturated (t_2 - t_c).
		phi:   Grid of saturation capacities.
		k_s:   Soil residence time.
		
	Returns:
		Array.
    """
    sat_vol = k_s*phi*t_sat

    return sat_vol

def overland_flow(R, t_sat, phi, k_s):
    """ Calculates the overland flow when the water level is above the 
        saturation capacity.
	
	Args:
		R:     Grid of inputs (essentially HER).
		t_sat: Time for which soil is saturated (t_2 - t_c).
		phi:   Grid of saturation capacities.
		k_s:   Soil residence time.
		
	Returns:
		Array.
    """
    of = (R - k_s*phi)*t_sat

    return of

def t_zero(R, S_0, k_s):
    """ Calculates the time when water level equals zero. Only relevant if
        R < 0.
		
	Args:
		R:   Grid of inputs (essentially HER).
		S_0: Grid of initial soil outflows.
		k_s: Soil residence time.
		
	Returns:
		Array.
    """
    import numpy as np
    
    t_0 = np.log(1 - (S_0/R))/k_s
   
    return t_0

def reduce_et_if_dry(wat_lev, fc, excess_et):
    """ Function to reduce ET when the soil is dry. Based on exponential decay
        of the form y = A*exp(-kx). If assume that at FC, 1 mm of ET reduces 
        the water level by 1 mm (i.e. ET takes place normally), this provides
        an additional constraint and the equation for the decay curve becomes:
        
                wat_lev = FC*exp(-ET/FC)

	Args:
		wat_lev:   Grid of current water levels.
		fc:        Grid of field capacities.
		excess_et: Grid of 'raw' PET values.
		
	Returns:
		Arrays [actual_et, corrected_water_level].
    """
    import numpy as np

    # Calculate how much drying has already taken place to reach this point
    et_0 = -1.*fc*np.log(wat_lev/fc)

    # Add this to the additional ET for this step
    et_1 = et_0 + excess_et

    # Calculate new water level
    new_wat_lev = fc*np.exp(-1.*et_1/fc)

    # Actual amount of ET after reduced
    act_et = wat_lev - new_wat_lev    
       
    return act_et, new_wat_lev
            
def calculate_drainage(R, surf_lev, gw_lev, fc, days_in_period, phi, k_s, 
                       k_g, b):
    """ Calculates drainage via OF, SSF and GWF pathways. Also calculates S_1
        and D_1, the soil and groundwater flows at the end of the time step.
		
	Args:
		R:              Grid of inputs (essentially HER). 
		surf_lev:       Soil wtaer level.   
		gw_lev:         Groundwater level.
		fc:             Field capacity grid.
		days_in_period: Length of time step (i.e. of month) in days.
		phi:            Saturation capacity grid.
		k_s:            Soil residence time.
		k_g:            Groundwater residence time.
		b:              Base Flow Index.
		
	Returns:
		Tuple of arrays: (net_her, of, ssf, gwf, new_surf_lev, new_gw_lev).
    """
    import numpy as np

    # Make a copies of R so that modification doesn't affect HER gris in main 
    # script
    R1 = R.copy()
    R2 = R.copy()
        
    # 1. First consider only cells where R >= 0 i.e. all cells where R < 0 
    # should be NaN in these calculations

    # 1.1. Cacl time to reach FC
    t_fc = (fc - surf_lev)/R1 # days

    # If t_fc is less than zero, we're already at or above fc
    t_fc[t_fc<0] = 0
    t_fc[R1<0] = np.nan

    # 1.1a. Added 09/12/2014 as initially forgot groundwater drainage while
    # soil level blow FC. Calc groundwater drainage when soil level below fc
    # Get time in step where soil level below FC
    t_fc2 = t_fc.copy()
    t_fc2[t_fc2>days_in_period] = days_in_period

    # Get gw outflow rate from previous step    
    D_0 = gw_lev*k_g
    D_0[R1<0] = np.nan

    # Calc outflow after t_fc2    
    gw_out_0 = gw_outflow_equn(t_fc2, 0, 0, D_0, k_s, k_g, b)
    gw_out_0[R1<0] = np.nan

    # Cacl volume after t_fc2
    gw_vol_0 = gw_vol(0, 0, D_0, 0, t_fc2, k_s, k_g, b)
    gw_vol_0[R1<0] = np.nan
    
    # 1.2. If t_fc is greater than days_in_period, the water doesn't reach fc 
    # in this step. We can calculate the level simply by adding water at rate
    # R for length of period. 
    temp_wat_lev_1 = surf_lev.copy()
    temp_wat_lev_1[t_fc>days_in_period] = (surf_lev[t_fc>days_in_period] +
                                           (R1*days_in_period)[t_fc>days_in_period])
    temp_wat_lev_1[t_fc<=days_in_period] = np.nan
    temp_wat_lev_1[R1<0] = np.nan

    # For the cells, where the water level doesn't reach fc, the processing is
    # finished for this time step. However, none zero values for R in these
    # cells will affect subsequent processing. Set R to zero for these cells,
    # then patch water levels in from temp_wat_lev_1 at end.
    R1[t_fc>days_in_period] = 0
    
    # 1.3. For those cells initially below fc, but which reach fc in this time 
    # step, set the levl to fc.
    temp_wat_lev_2 = surf_lev.copy()
    temp_wat_lev_2[t_fc>0] = fc[t_fc>0]
    temp_wat_lev_2[R1<0] = np.nan
    
    # 1.4. Calc time in period for which levels are at or above fc
    t_rem = days_in_period - t_fc
    t_rem[t_rem<0] = 0
    t_rem[R1<0] = np.nan

    # 1.5. Calc time to sc
    S_0 = k_s*(temp_wat_lev_2 - fc)
    S_0[S_0<0] = 0
    S_0[R1<0] = np.nan
    t_c = t_critical(R1, S_0, phi, k_s)
    t_c[R1<0] = np.nan

    # 1.6. Calc time spent in unsaturated zone
    t_fc_sc = np.minimum(t_rem, t_c)

    # 1.7. Calc soil outflow rate at end of unsaturated portion (same as 
    # outflow rate at end of time step)
    soil_out_1 = soil_outflow_equn(t_fc_sc, R1, S_0, k_s)
    soil_out_1[R1<0] = np.nan

    # 1.8. Calc groundwater flow at end of unsaturated period   
    gw_out_1 = gw_outflow_equn(t_fc_sc, R1, S_0, gw_out_0, k_s, k_g, b)
    gw_out_1[R1<0] = np.nan
    
    # 1.9. Calc soil outflow volume while unsaturated
    soil_vol_1 = unsaturated_vol(R1, S_0, 0, t_fc_sc, k_s)
    soil_vol_1[R1<0] = np.nan
    
    # 1.10. Calc groundwater volume while unsaturated
    gw_vol_1 = gw_vol(R1, S_0, gw_out_0, 0, t_fc_sc, k_s, k_g, b)
    gw_vol_1[R1<0] = np.nan

    # 1.11. Calc the time for which the soil is saturated
    t_sat = t_rem - t_c
    t_sat[t_sat<0] = 0
    t_sat[R1<0] = np.nan

    # 1.12. Calc soil drainage while saturated
    soil_vol_2 = saturated_vol(t_sat, phi, k_s)
    soil_vol_2[R1<0] = np.nan
    
    # 1.13. Calc overland flow
    of1 = overland_flow(R1, t_sat, phi, k_s)
    of1[R1<0] = np.nan
    
    # 1.14. Calc groundwater drainage volume while saturated
    gw_vol_2 = gw_vol(k_s*phi, k_s*phi, gw_out_1, 0, t_sat, k_s, k_g, b)
    gw_vol_2[R1<0] = np.nan

    # 1.15. Calc groundwater outflow at end of time step
    gw_out_2 = gw_outflow_equn(t_sat, k_s*phi, k_s*phi, gw_out_1, k_s, k_g, b)
    gw_out_2[R1<0] = np.nan
    
    # 1.16. Sum volumes
    ssf1 = (1 - b)*(soil_vol_1 + soil_vol_2)
    gwf1 = gw_vol_0 + gw_vol_1 + gw_vol_2

    # 1.17. Water levels at end of time step
    new_gw_lev1 = gw_out_2/k_g
    new_surf_lev1 = soil_out_1/k_s + fc
    
    # Path levels from cells still below field capacit (step 1.2)
    new_surf_lev1[t_fc>days_in_period] = temp_wat_lev_1[t_fc>days_in_period]

    # 1.18. Calc actual HER for this step
    her1 = R*days_in_period
    her1[R<0] = np.nan
    
    # 2. Now consider cells where R < 0.
    
    # 2.1. Calculate S_0 where water level is above FC
    S_0 = (surf_lev - fc)*k_s
    S_0[S_0<0] = 0
    S_0[R2>=0] = np.nan
    
    # 2.2. Calculate the time taken to drain to FC. For cells with water above
    # FC, this must be a positive real number (because R<0). Values for t_0
    # less than zero will be handled separately, so set to zero here
    t_0 = t_zero(R2, S_0, k_s)
    t_0[R2>=0] = np.nan

    # 2.3. Calc time between FC and SC
    t_fc_sc = np.minimum(days_in_period, t_0)
    t_fc_sc[R2>=0] = np.nan

    # 2.4. Calc soil outflow at end of step
    soil_out_1 = soil_outflow_equn(t_fc_sc, R2, S_0, k_s)
    soil_out_1[t_0<days_in_period] = 0 # Avoids floating point errors later    
    soil_out_1[R2>=0] = np.nan

    # 2.5. Calc soil drainage volume
    soil_vol_1 = unsaturated_vol(R2, S_0, 0, t_fc_sc, k_s)
    soil_vol_1[soil_vol_1<0] = 0
    soil_vol_1[R2>=0] = np.nan

    # 2.6. Calc groundwater outflow between SC and FC
    D_0 = gw_lev*k_g
    D_0[R2>=0] = np.nan
    gw_out_1 = gw_outflow_equn(t_fc_sc, R2, S_0, D_0, k_s, k_g, b)
    gw_out_1[gw_out_1<0] = 0
    gw_out_1[R2>=0] = np.nan

    # 2.7. Calc groundwater volume between SC and FC
    gw_vol_1 = gw_vol(R2, S_0, D_0, 0, t_fc_sc, k_s, k_g, b)
    gw_vol_1[gw_vol_1<0] = 0
    gw_vol_1[R2>=0] = np.nan

    # 2.8. Calc time at or below FC
    t_dry = days_in_period - t_fc_sc

    # 2.9. For cells where water level is initially above FC, but which reach
    # FC in this time step, set the level to FC. These cells have 
    # days_in_period > t_dry > 0.
    temp_wat_lev_1 = surf_lev.copy()
    temp_wat_lev_1[(t_dry>0)&(t_dry<days_in_period)] = fc[(t_dry>0) & 
                                                          (t_dry<days_in_period)]
    temp_wat_lev_1[R2>=0] = np.nan

    # 2.10. Calc excess ET
    ex_et = -1.*R2*t_dry
    ex_et[R2>=0] = np.nan

    # 2.11. Reduce ET where dry
    act_et, new_surf_lev2 = reduce_et_if_dry(temp_wat_lev_1, fc, ex_et)
    
    # 2.12. Patch in values where level is still above FC at end of step
    temp_wat_lev_2 = (soil_out_1/k_s) + fc
    new_surf_lev2[soil_out_1>0] = temp_wat_lev_2[soil_out_1>0]
    new_surf_lev2[R2>=0] = np.nan

    # 2.13. Calc groundwater outflow at end of step
    gw_out_2 = gw_outflow_equn(t_dry, 0, 0, gw_out_1, k_s, k_g, b)
    gw_out_2[gw_out_2<0] = 0
    gw_out_2[R2>=0] = np.nan

    # 2.14. Calc groundwater volume once soil dry
    gw_vol_2 = gw_vol(0, 0, gw_out_1, 0, t_dry, k_s, k_g, b)
    gw_vol_2[gw_vol_2<0] = 0
    gw_vol_2[R2>=0] = np.nan

    # 2.15. OF is zero for all cells where R < 0
    of2 = R2.copy()
    of2[R2<0] = 0
    of2[R2>=0] = np.nan

    # 2.16. Sum volumes
    ssf2 = (1 - b)*soil_vol_1
    gwf2 = gw_vol_1 + gw_vol_2

    # 2.17. Groundwater level
    new_gw_lev2 = gw_out_2/k_g 

    # 2.18. Calculate the actual HER for this time step
    her2 = R2*t_fc_sc - act_et
    her2[R2>=0] = np.nan
    
    # 3. Combine results from the two stages
    # 3.1. OF
    of = np.dstack((of1, of2))
    of = np.nanmax(of, axis=2)
    
    # 3.2. SSF
    ssf = np.dstack((ssf1, ssf2))
    ssf = np.nanmax(ssf, axis=2)
    
    # 3.3. GWF
    gwf = np.dstack((gwf1, gwf2))
    gwf = np.nanmax(gwf, axis=2) 
    
    # 3.4. Surface level
    new_surf_lev = np.dstack((new_surf_lev1, new_surf_lev2))
    new_surf_lev = np.nanmax(new_surf_lev, axis=2)

    # 3.5. Groundwater level
    new_gw_lev = np.dstack((new_gw_lev1, new_gw_lev2))
    new_gw_lev = np.nanmax(new_gw_lev, axis=2)
    
    # 3.6. HER
    net_her = np.dstack((her1, her2))
    net_her = np.nanmax(net_her, axis=2)

    return (net_her, of, ssf, gwf, new_surf_lev, new_gw_lev)