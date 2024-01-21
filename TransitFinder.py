# -*- coding: utf-8 -*-
"""
Functions to find a possible transit:

- getTransit()
- isWing()
- getTransitTime()
- getCircOrbitData()
- enoughData()
- getTransitdFlux()

DATE        AUTHOR          Modification
01/05/2020  Marcel Winzen   Derived from CheckFirstSol
21/05/2020  Marcel Winzen   Added getTransitdFlux
25/06/2020  Marcel Winzen   Added getTransitneg and isWingneg
11/11/2020  Marcel Winzen   Two different approaches for twing 0.02 and twing 0.04
16/11/2020  Marcel Winzen   PCD Error sigma considered
26/11/2020  Marcel Winzen   Added approach for twing = 0.06
09/12/2020  Marcel Winzen   Added new calculation of delta_F
11/12/2020  Marcel Winzen   Added approach for twing = 0.08
03/10/2021  Marcel Winzen   f_1 is dependent on percentile instead of error
03/10/2021  Marcel Winzen   Both limb darking durations must be of same size
"""

import numpy as np
import math
import astropy.constants as co
from DataUseTools import * 


__all__ = ['getTransit', 'isWing', 'getTransitTime',
           'getCircOrbitData','enoughData','getTransitdFlux']

def getTransit(i_start, i_last, t_SMA, f_2, f_1,  df):
    """
    This function searches for a possible planet transit
    and returns the breakpoints of it
    """ 
    wingpoint = 0 
    t_start = df[i_start]['TIME']
    t_end = df[i_last]['TIME']
    index_start = getClosestNeighbour((t_start), 'TIME', df) 
    index_end   = getClosestNeighbour((t_end-t_SMA), 'TIME', df) - 2
    c_0 = i_start
    c_1 = i_start
    wingflag = False
    for i in range(index_start, index_end):
        wingflag = isWing(i, t_SMA, f_2, f_1, wingpoint, df)
        if wingflag == True:
            wingpoint = 1
            c_0 = i
            c_1= c_0
            break
    if wingflag == False:
        c_0 = index_end
        c_1= index_end
    
    if wingpoint == 1:
        c1max_end = getClosestNeighbour((df[c_0]['TIME']+3*t_SMA), 'TIME', df)
        for i in range(c_0, c1max_end):
            wingflag = isWing(i, t_SMA, f_2, f_1, wingpoint, df)
            if wingflag == True:
                wingpoint = 2
                c_1 = i
                break

            
    if wingpoint == 2:
        c2max_end = getClosestNeighbour((df[c_1]['TIME']+10*t_SMA), 'TIME', df)
        for i in range(c_1+1, c2max_end):
            wingflag = isWing(i, t_SMA, f_2, f_1, wingpoint, df)
            if wingflag == True:
                wingpoint = 3
                c_2 = i
                break

            
    if wingpoint == 3:
        c3min_start = getClosestNeighbour((df[c_2]['TIME']+(df[c_1]['TIME']-df[c_0]['TIME'])/2), 'TIME', df) 
        c3max_end = getClosestNeighbour((df[c_2]['TIME']+(df[c_1]['TIME']-df[c_0]['TIME'])*2), 'TIME', df) 
        if c3min_start < c_2+1:
            c3min_start = c_2+1
        if c3max_end < c3min_start+1:
            c3max_end = c3min_start+1
        for i in range(c3min_start, c3max_end):
            wingflag = isWing(i, t_SMA, f_2, f_1, wingpoint, df)
            if wingflag == True:
                wingpoint = 4
                c_3 = i
                break

            
    if wingpoint != 4:
        #Dummy initialisation to fail continue
        c_2= 0
        c_3 = 0
                    

    return c_0, c_1, c_2, c_3

def isWing(index, t_SMA, f_2, f_1, wingpoint, df):
    """
    This function searches for the characteristic breakpoints 
    of the transit brightness profile
    """    
    wingflag = False
    t_i = df[index]['TIME']   
    h = getClosestNeighbour((t_i-t_SMA), 'TIME', df)
    k = getClosestNeighbour((t_i+t_SMA), 'TIME', df)
    if (h == k) or not df[k]['PDCSAP_FLUX']or not df[h]['PDCSAP_FLUX']:
        return wingflag  
    
    if (t_SMA > 0.09): #meaning t_SMA = 0.1 days
        f_1fac= 1-0.2*f_1 #c_SMA =0.2
        F_i_Av = Average(df[(index-2):(index+2)]['PDCSAP_FLUX'])
        F_h_Av = Average(df[(h-2):(h+2)]['PDCSAP_FLUX'])
        F_k_Av = Average(df[(k-2):(k+2)]['PDCSAP_FLUX'])
    elif (t_SMA > 0.07) and (t_SMA < 0.09): #meaning t_SMA = 0.08 days
        f_1fac= 1-0.25*f_1 #c_SMA =0.25
        F_i_Av = Average(df[(index-1):(index+2)]['PDCSAP_FLUX'])
        F_h_Av = Average(df[(h-1):(h+2)]['PDCSAP_FLUX'])
        F_k_Av = Average(df[(k-1):(k+2)]['PDCSAP_FLUX'])
    elif (t_SMA > 0.05) and (t_SMA < 0.07): #meaning t_SMA = 0.06 days
        f_1fac= 1-0.33*f_1  #c_SMA =0.33
        F_i_Av = Average(df[(index-1):(index+1)]['PDCSAP_FLUX'])
        F_h_Av = Average(df[(h-1):(h+1)]['PDCSAP_FLUX'])
        F_k_Av = Average(df[(k-1):(k+1)]['PDCSAP_FLUX'])

  
    if wingpoint == 0:
        if ((F_k_Av - F_i_Av) < (f_2 * (F_i_Av - F_h_Av))) and (F_k_Av < (f_1fac* F_i_Av)):
            wingflag = True
    elif wingpoint == 1:
        if ((F_i_Av - F_h_Av) < (f_2 * (F_k_Av - F_i_Av))) and ((f_1fac*F_k_Av) < F_i_Av):
            wingflag = True
    elif wingpoint == 2:
        if ((F_k_Av - F_i_Av) > (f_2 * (F_i_Av - F_h_Av))) and ((f_1fac*F_k_Av) > F_i_Av):
            wingflag = True
    elif wingpoint == 3:
        if ((F_i_Av - F_h_Av) > (f_2 * (F_k_Av - F_i_Av))) and (F_k_Av > (f_1fac* F_i_Av)):
            wingflag = True
            
    return wingflag
             
 
        

def getTransitTime(c0, c1, c2, c3, df):
    """
    This function is to refine the start and end of the transit found
    and to receive the transit time tran_duration [sec]
    """ 

    # Calculate the transit time tran_duration
    t_0 = df[c0]['TIME'] * 86400  # in sec
    t_1 = df[c3]['TIME'] * 86400  # in sec
    tran_duration = t_1 - t_0 # in sec
    return tran_duration


def getCircOrbitData(tran_duration, R_star, M_star):
    """
    This function is to iterate the best solution for planet orbit 
    period T (sec). In this step a circular orbit 
    is assumed and that the mass of the star is much bigger 
    than the mass of the planet. 
    """ 
    # With geometry, Kepler's third rule (M_star >> m_p)
    # and r_p >> R_Star
    # the period can be calculated from:
    T = (co.G.value/co.R_sun.value**3*co.M_sun.value) * M_star * (tran_duration**3) * math.pi / (4 * R_star**3)
    return T
        
def enoughData(T, df):
    """
    This function is to confirm that the time range of data
    is more than 3 times the planet orbit period T to ensure
    enough data to confirm a solution
    """    
    last_index = len(df) - 1 
    t_start = df[0]['TIME']
    t_end = df[last_index]['TIME']
    t_data = t_end - t_start
    if t_data > (3 * T/86400):
        return True
    else:
        return False 
    
def getTransitdFlux(c0, c1, c2, c3, df):
    """
    Get the Delta Light Flux in order to determine the planet radius
    """ 
    dPDCSAP =  (Average(df[c3:(c3+1)]['PDCSAP_FLUX']) + Average(df[(c0-1):c0]['PDCSAP_FLUX']))/2-Average(df[c1:c2]['PDCSAP_FLUX'])
    delta_F = dPDCSAP / ((Average(df[c3:(c3+1)]['PDCSAP_FLUX']) + Average(df[(c0-1):c0]['PDCSAP_FLUX']))/2)
    return delta_F

