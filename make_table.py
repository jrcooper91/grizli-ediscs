import glob
import time
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import astropy.io.fits as pyfits
import drizzlepac

import grizli

directory = '.'

for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith("full.fits"): 
         fit_hdu = pyfits.open(filename) 
         ids = fit_hdu[0].header['ID']
         ra = fit_hdu[0].header['RA']
         dec = fit_hdu[0].header['DEC']
         z = fit_hdu[0].header['REDSHIFT']
         Ha = fit_hdu[2].header['FLUX_004']
         Ha_err = fit_hdu[2].header['ERR_004']
         sys.stdout = open("t.txt", "a")
         print(ids,ra,dec,z,Ha,Ha_err) 
         continue
     else:
         continue
