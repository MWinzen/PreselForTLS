# -*- coding: utf-8 -*-
"""
Functions to prepare data for different use:

- FilterSAPQual()
- FilterPDCErr()
- Average()
- getClosestNeighbour()
- pairTransits()
- removesimilarSol()
- processSol()

DATE        AUTHOR          Modification
01/05/2020  Marcel Winzen   Derived from FirstTransitFcn
16/05/2020  Marcel Winzen   Added removesimilarSol
21/05/2020  Marcel Winzen   Added processSol
06/11/2020  Marcel Winzen   Added new Solution Delivery
11/11/2020  Marcel Winzen   3 waves need to be found for recognition
12/11/2020  Marcel Winzen   removesimalarSol filters adapted
18/11/2020  Marcel Winzen   Change for sampling
29/05/2021  Marcel Winzen   Added pairTransits
"""

import numpy as np
import math
import astropy.constants as co
from astropy.table import Column, Table, vstack
import progressbar


__all__ = ['FilterTimeErr','FilterPDCErr','FilterSAPQual', 'Average', 'getClosestNeighbour', 'pairTransits', 'removesimilarSol', 'processSol']

def FilterPDCErr(df):
    """
    This function is to check that the PDC Err has a value
    """    
    df = df[np.isfinite(df['PDCSAP_FLUX_ERR'])]
    return df

def FilterTimeErr(df):
    """
    This function is to check that the Time has a value
    """    
    df = df[np.isfinite(df['TIME'])]
    return df

def FilterSAPQual(df):
    """
    This function is to check that the SAP Quality Flag for the whole table df
    If the SAP_QUALITY is other than 0 it is deleted from df
    """    
    mask = df['SAP_QUALITY'] == 0
    df = df[mask]
    return df
    
def Average(lst): 
    return sum(lst) / len(lst) 

def getClosestNeighbour(inputvalue, column, df ):
    """
    This function is to find the closest value to an input
    within a column and returns the index
    """
    index = (np.abs(inputvalue-df[column])).argmin() 
    return index #, df[index][column]

def pairTransits(df, t_data):
    """
    This function is to find similar solutions and 
    keep only the first one in time
    """
    #
    df.sort('t_0')
    Solutions = Table([[],[],[],[],[],[],[]], names=('T', 'tran_duration', 'delta_F', 't_0', 'Amount Transits','SDE','SDE_Dis'))
    print('Pair transits to find possible orbit periods:')
    bar = progressbar.ProgressBar(maxval=(len(df)), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
                                                
    # pair transit
    for i in range(0, len(df) ):
        for i2 in range (i+1, len(df)):
            if ((df[i]['tran_duration'] > 0.8 * df[i2]['tran_duration']) and (df[i]['tran_duration'] < 1.25 * df[i2]['tran_duration']) ):
                if ((df[i]['delta_F'] > 0.6 * df[i2]['delta_F']) and (df[i]['delta_F'] < 1.67 * df[i2]['delta_F']) ):
                    if (( abs(df[i2]['t_0']-df[i]['t_0']) > 50.0) and (abs(df[i2]['t_0']-df[i]['t_0'])<t_data/3)):
                                T_new = abs(df[i2]['t_0']-df[i]['t_0'])
                                tran_duration_new= Average([df[i]['tran_duration'] ,df[i2]['tran_duration']])
                                delta_F_new= Average([df[i]['delta_F'] ,df[i2]['delta_F']])            
                                Amounttransits= df[i]['Amount Transits'] +df[i2]['Amount Transits']
                                Solutions.add_row([T_new, tran_duration_new, delta_F_new, df[i]['t_0'], Amounttransits, 0.0, 0.0])
        bar.update(i+1)
    bar.finish()       
    return Solutions                 


def removesimilarSol(df):
    """
    This function is to find similar solutions and 
    keep only the first one in time
    """
    #
    df.sort('T')        

    # Collect transits of same exoplanet
    print('Collect transits of same exoplanet:')
    bar = progressbar.ProgressBar(maxval=(len(df)), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
        
    for i in range(0, len(df) ):
        if df[i]['T'] != 0:
            for i2 in range (i+1, len(df)):
                if df[i2]['T'] != 0:
                    if abs(df[i2]['T']-df[i]['T']) < (df[i]['tran_duration']/4):
                        res = abs(df[i2]['t_0']-df[i]['t_0'])% df[i]['T']
                        mult =  1.1+round(abs(df[i2]['t_0']-df[i]['t_0'])/ df[i]['T'])
                        if ( res< (df[i]['tran_duration']*(2-2/mult))) or (res> (df[i]['T']-df[i]['tran_duration']*(2-2/mult))):
                            res = abs(df[i2]['t_0']-df[i]['t_0'])% df[i2]['T']
                            mult = 1.1+round(abs(df[i2]['t_0']-df[i]['t_0'])/ df[i2]['T'])
                            if ( res< (df[i2]['tran_duration']*(2-1/mult))) or (res> (df[i2]['T']-df[i2]['tran_duration']*(2-1/mult))):
                                df[i2]['T'] = 0.0
                                df[i]['Amount Transits']= df[i]['Amount Transits'] +df[i2]['Amount Transits']
                                df[i]['tran_duration']= Average([df[i]['tran_duration'] ,df[i2]['tran_duration']])
                                df[i]['delta_F']= Average([df[i]['delta_F'] ,df[i2]['delta_F']])
                                
        bar.update(i+1)
    bar.finish()  
    mask = ( df['T'] != 0.0)
    df = df[mask]                    
                     


    return df, False

def processSol(df, R_star, M_star):
    """
    This function is to process the data:
    Time should be in days instead of seconds
    Also calculate Planet data r_p and R_p
    """
    df['T'] = df['T']  *86400
    df['tran_duration'] = df['tran_duration']*86400  
    df['t_0'] = df['t_0']*86400
    r_p = Column(np.arange(len(df)), name='r_p')
    df.add_column(r_p, index=2)
    R_p = Column(np.arange(len(df)), name='R_p')
    df.add_column(R_p, index=2)
    df['r_p'] = ((df['T'])**2 * ((co.G.value*co.M_sun.value) * M_star)/(4 * (math.pi)**2))**(1/3) / co.au
    df['R_p'] = ((R_star*co.R_sun.value)**2*df['delta_F'])**(1/2)/co.R_earth.value
    df['T'] = df['T']  /86400
    df['tran_duration'] = df['tran_duration']/86400  
    df['t_0'] = df['t_0'] /86400
    # df.remove_column('delta_F')
    return df 
