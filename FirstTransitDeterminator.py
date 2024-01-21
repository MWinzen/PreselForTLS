# -*- coding: utf-8 -*-
"""
Functions to find a first possible transit:

- getfirstTransitData()

DATE        AUTHOR          Modification
01/05/2020  Marcel Winzen   Initial version
21/05/2020  Marcel Winzen   Also retrieve Transit PCD_Flux difference
01/11/2020  Marcel Winzen   Scale tran_duration to 25 days for efficiency reason
09/12/2020  Marcel Winzen   Added new filter of transit by delta_F
03/10/2021  Marcel Winzen   f_1 is dependent on percentile instead of error
"""

from TransitFinder import * 
from DataUseTools import * 

__all__ = ['getfirstTransitData']



def getfirstTransitData(index, t_SMA, f_2, f_1, R_star, M_star, df):
    """
    This function finds breakpoints of a first transit,
    checks if it is realistic and if there is enough data to verify it
    """ 
    # define delta_Fmin criteria 
    if (t_SMA > 0.09):
        delta_Fmin = 0.2*f_1
    elif (t_SMA > 0.07) and (t_SMA < 0.09):
        delta_Fmin = 0.25*f_1
    elif (t_SMA > 0.05) and (t_SMA < 0.07):      
        delta_Fmin = 0.33*f_1
    c0, c1, c2, c3 = getTransit(index,(len(df) - 1), t_SMA, f_2, f_1, df)


    if (c3 != 0):
        delta_F = getTransitdFlux(c0, c1, c2, c3, df)
        if delta_F > delta_Fmin:
            tran_duration = getTransitTime(c0, c1, c2, c3, df)
            if ((tran_duration/86400 < 25.0) and (tran_duration/86400 > 3*t_SMA)):
                #T = getCircOrbitData(tran_duration, R_star, M_star)
                T = 0
                #Continue = enoughData(T, df)
                return True, T, tran_duration, delta_F,c0, c3,c1
    return False, 0, 0, 0, c0, 0, c1

