# A spatially distributed water balance model for Scotland

This repository contains code for a simple, spatially distributed water balance model. Calibration of the model parameters is achieved using parallel-tempered MCMC. The model was developed by the [James Hutton Institute](http://www.hutton.ac.uk/) and has been used to inform hydropower research. Associated work can be found [here](http://www.sciencedirect.com/science/article/pii/S1364032115007182) and [here](http://meetingorganizer.copernicus.org/EGU2015/EGU2015-3351.pdf).

## Requirements

  * [Python 2.7](https://www.python.org/download/releases/2.7/)
  * [Numpy](http://www.numpy.org/)
  * [Pandas](http://pandas.pydata.org/)
  * [H5py](http://www.h5py.org/)
  * [GDAL](http://www.gdal.org/)
  * [emcee](http://dan.iel.fm/emcee/current/)
  * [corner](https://github.com/dfm/corner.py)

If you're installing from scratch on Windows then I thoroughly recommend the [WinPython](http://winpython.github.io/) installation.

## Documentation

The model is still under development, but a brief introduction to the basic structure and calibration process can be found [here]().
