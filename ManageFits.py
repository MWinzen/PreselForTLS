# -*- coding: utf-8 -*-
"""
Functions to manage Kepler Fits Files from MAST observational data:

- DownloadFits()
- GetListofFitsAvail()

DATE        AUTHOR          Modification
07/06/2020  Marcel Winzen   Initial Version
05/07/2020  Marcel Winzen   Added reserve Downlink Funktion
"""
import pathlib
import time
from astroquery.mast import Mast
from astroquery.mast import Observations


__all__ = ['DownloadFits', 'GetListofFitsAvail']

def DownloadFits(KeplerId,FileFormat):
    """
    This function is to download all fits files available of a given KeplerID
    """    

    FilesAvailable = True
    KeplerIdstr = 'kplr' + '{0:09}'.format(KeplerId)
    keplerObs = Observations.query_criteria(target_name=KeplerIdstr, obs_collection='Kepler')
    keplerProds = Observations.get_product_list(keplerObs)
    if FileFormat == "SLC":
        yourProd = Observations.filter_products(keplerProds, extension='slc.fits',  mrp_only=False)
    if FileFormat == "LLC":
        yourProd = Observations.filter_products(keplerProds, extension='llc.fits',  mrp_only=False)
    
    if len(yourProd) == 0:
        FilesAvailable = False
    Observations.download_products(yourProd, mrp_only = False, cache = False)
 
    
    return FilesAvailable

def GetListofFitsAvail(KeplerId):
    """
    Get a list of all Fits-Files available for a Kepler-ID
    """    
    KeplerIdstr = 'kplr' + '{0:09}'.format(KeplerId)
    FitsFiles = list(pathlib.Path('./mastDownload/Kepler').glob('*/*'+KeplerIdstr +'*.fits'))
    FitsFiles.sort(reverse=False)
    return FitsFiles
  
