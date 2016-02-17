#------------------------------------------------------------------------------
# Name:        input_output.py
# Purpose:     Functions to handles I/O for the new version of the water 
#              balance model
#
# Author:      James Sample
#
# Created:     25/07/2014
# Copyright:   (c) James Sample and JHI, 2014
#------------------------------------------------------------------------------
""" Functions for handling I/O operations from the main HDF5 file.
"""

def read_ascii(ascii_path,
               xmin=0,
               xmax=485000,
               ymin=520000,
               ymax=1235000,
               exptd_rows=715,
               exptd_cols=485,
               exptd_px_wd=1000,
               exptd_px_ht=-1000,
               exptd_ndv=-9999):
    """ Reads an ASCII grid file, clips it to the specified bounding box and
        returns a numpy array.

	Args:
		ascii_path:  Path to ASCII grid file
		xmin:        OSGB-1936 metres for minimum Easting of area of interest.
		xmax:        OSGB-1936 metres for maximum Easting of area of interest.
		ymin:        OSGB-1936 metres for minimum Northing of area of interest.
		ymax:        OSGB-1936 metres for maximum Northing of area of interest.
		exptd_rows:  Expected number of rows in array.
		exptd_cols:  Expected number of columns in array.
		exptd_px_wd: Expected pixel width in metres.
		exptd_px_ht: Expected pixel height in -1*metres.
		exptd_ndv:   Expected no data value.
	
	Returns:
		Array.
    """
    from osgeo import gdal, gdalconst

    # Register drivers
    gdal.AllRegister()

    # Process the file with GDAL
    ds = gdal.Open(ascii_path, gdalconst.GA_ReadOnly)
    if ds is None:
    	print 'Could not open ' + ascii_path
    	sys.exit(1)

    # In order to select the first cell correctly, choose a point just within
    # the top left corner of the specified bounding box.
    x = xmin + 10
    y = ymax - 10

    # Dataset properties
    geotransform = ds.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    # Calculate number of rows and cols to return
    rows = abs(int((ymax-ymin)/pixelHeight))
    cols = int((xmax-xmin)/pixelWidth)

    # Select starting pixel
    xOffset = int((x - originX) / pixelWidth)
    yOffset = int((y - originY) / pixelHeight)

    band = ds.GetRasterBand(1)
    no_data_val = band.GetNoDataValue()

    # Simple checking
    assert rows == exptd_rows
    assert cols == exptd_cols
    assert pixelWidth == exptd_px_wd
    assert pixelHeight == exptd_px_ht
    assert no_data_val == exptd_ndv

    # Read the data to an array
    data = band.ReadAsArray(xOffset, yOffset, cols, rows)

    # Close the dataset
    ds = None

    return data
	
def read_array_from_h5(h5, dset_path, slice_idx, ndv, xmin_idx, xmax_idx, 
                       ymin_idx, ymax_idx):
    """ Reads an array from the specified location in an H5 file. Clips it to
        the desired co-ordinates and sets NoData values to NaN.
	
	Args:
	    h5:        Path to the (open) h5 object to read from.
		dset_path: Path to dataset of interest in h5. 
		slice_idx: Time slice. If the array of interest is 2D (not 3D), pass 
		           None.
		ndv:       No data value. 
		xmin_idx:  Array index for minimum Easting.
		xmax_idx:  Array index for maximum Easting.
		ymin_idx:  Array index for minimum Northing.
		ymax_idx:  Array index for maximum Northing.
	
	Returns:
		Array.
    """
    import numpy as np
	
    if slice_idx == None:
        data = h5.get(dset_path)[ymin_idx:ymax_idx, 
                                 xmin_idx:xmax_idx].astype(float)
    else:
        data = h5.get(dset_path)[ymin_idx:ymax_idx, 
                                 xmin_idx:xmax_idx,
                                 slice_idx].astype(float)
    
    # Set NoData to NaN
    data[data==ndv] = np.nan
    
    return data
    
def read_soil_properties(h5, xmin_idx, xmax_idx, ymin_idx, ymax_idx):
    """ Returns the FC, SC and BFI arrays.
	
		Args:
	    h5:        Path to the (open) h5 object to read from.
		xmin_idx:  Array index for minimum Easting.
		xmax_idx:  Array index for maximum Easting.
		ymin_idx:  Array index for minimum Northing.
		ymax_idx:  Array index for maximum Northing.
	
	Returns:
		List of arrays [fc, sat, bfi].
    """
    # Empty list to store arrays
    grid_list = []

    # List of desired grids to extract
    wanted_grids = ['fc', 'sc', 'bfi']

    for grid in wanted_grids:
        data = read_array_from_h5(h5, r'/soil_properties/%s' % grid, None,
                                  -99, xmin_idx, xmax_idx, ymin_idx, 
                                  ymax_idx)
        grid_list.append(data)
    
    return grid_list # [fc, sat, bfi]

def read_land_use(h5, lu_grid, xmin_idx, xmax_idx, ymin_idx, ymax_idx):
    """ Returns the PET factors for the desired land use grid.
	
	Args:
	    h5:       Path to the (open) h5 object to read from.
		lu_grid:  Path to land use dataset of interest. 
		xmin_idx: Array index for minimum Easting.
		xmax_idx: Array index for maximum Easting.
		ymin_idx: Array index for minimum Northing.
		ymax_idx: Array index for maximum Northing.
	
	Returns:
		Array.
    """
    data = read_array_from_h5(h5, r'/land_use/%s' % lu_grid, None,
                              -99, xmin_idx, xmax_idx, ymin_idx, 
                              ymax_idx)
    
    return data
    
def read_met_data(h5, model, xmin_idx, xmax_idx, ymin_idx, ymax_idx,
                  year, month):
    """ Reads precipitation and PET grids for the specified month.
	
	Args:
	    h5:       Path to the (open) h5 object to read from.
		model:    'base', 'snow_corr_base' or a valid future model
		          code for the FF data (e.g. 'afixa' etc.). 
		xmin_idx: Array index for minimum Easting.
		xmax_idx: Array index for maximum Easting.
		ymin_idx: Array index for minimum Northing.
		ymax_idx: Array index for maximum Northing.
		year:     Year of interest.
		month:    Months of interest.
	
	Returns:
		Arrays [pptn, pet].
    """
    if model=='base':
        # Pptn
        pptn_dset = r'/mo_data/rainfall/rainfall_%s' % year
        pptn = read_array_from_h5(h5, pptn_dset, (month - 1), -99,
                                  xmin_idx, xmax_idx, ymin_idx, ymax_idx)                                
        
        # Convert units
        pptn = pptn/100.
        
        # PET                          
        pet_dset = r'/mo_data/pm_pet/pm_pet_%s' % year
        pet = read_array_from_h5(h5, pet_dset, (month - 1), -99,
                                 xmin_idx, xmax_idx, ymin_idx, ymax_idx)
        
        # Convert units
        pet = pet/100.
    
    elif model=='snow_corr_base':
        # Pptn
        pptn_dset = r'/mo_data/snow_corr_rainfall/snow_corr_rainfall_%s' % year
        pptn = read_array_from_h5(h5, pptn_dset, (month - 1), -99,
                                  xmin_idx, xmax_idx, ymin_idx, ymax_idx)                                
        
        # Convert units
        pptn = pptn/100.
        
        # PET                          
        pet_dset = r'/mo_data/snow_corr_pm_pet/snow_corr_pm_pet_%s' % year
        pet = read_array_from_h5(h5, pet_dset, (month - 1), -99,
                                 xmin_idx, xmax_idx, ymin_idx, ymax_idx)
        
        # Convert units
        pet = pet/100.
    
    else:
        # Pptn
        pptn_dset = r'/ff_data/%s/rainfall/rainfall_%s' % (model, year)
        pptn = read_array_from_h5(h5, pptn_dset, (month - 1), -99,
                                  xmin_idx, xmax_idx, ymin_idx, ymax_idx)                                
        
        # Convert units
        pptn = pptn/100.
        
        # PET                          
        pet_dset = r'/ff_data/%s/pet/pet_%s' % (model, year)
        pet = read_array_from_h5(h5, pet_dset, (month - 1), -99,
                                 xmin_idx, xmax_idx, ymin_idx, ymax_idx)
        
        # Convert units
        pet = pet/100.
    
    return pptn, pet