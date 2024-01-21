# -*- coding: utf-8 -*-
"""
Functions to find a first possible transit:

- doTLS()

DATE        AUTHOR          Modification
05/01/2021  Marcel Winzen   Initial version
02/05/2021  Marcel Winzen   Smoothen transitions between datasets
20/05/2021  Marcel Winzen   Create tls plots
16/06/2021  Marcel Winzen   Added group validation to accelare process
01/04/2022  Marcel Winzen   Only save SDE starting i_iter >=3
"""

from PlanetScanTLS import transitleastsquares
import os
import numpy as np
import math
import astropy.constants as co
from astropy.table import Table, vstack
import astropy.io.fits as iof
from astropy.io.votable import parse_single_table
from ManageFits import GetListofFitsAvail, DownloadFits
from DataUseTools import * 
from wotan import flatten
#from PlotsProducer import *

__all__ = ['doTLS']



def doTLS(Keplerid, path_to_csv, Tremove, t0remove,  M_star, R_star, confirmedplanets):
    """
    This function performes TLS for a specified 
    """ 



    samples = np.genfromtxt(path_to_csv+'Solutions_'+str(Keplerid)+ '.csv', dtype=float, delimiter=';', names=True)
    


    print('Download data with KeplerID: ' + str(Keplerid))
    # FilesAvailable = DownloadFits(Keplerid,'LLC')
    # if FilesAvailable == 0:
    #     print('No Planet Scan performed, since no data available for KeplerID: ' +str(Keplerid) )
    #     return
    
    FitsFiles = GetListofFitsAvail(Keplerid)
    
    print('Read in the light curve data and the star properties (radius, mass) from the Fits file given:')
    Fitsfile = FitsFiles[0]
    df = Table.read(Fitsfile)   
    df = FilterTimeErr(df)
    df = FilterSAPQual(df)
    df = FilterPDCErr(df)
    df['PDCSAP_FLUX']= df['PDCSAP_FLUX']/np.median(df['PDCSAP_FLUX'])
    
    for i_Files in range(1,(len(FitsFiles))):
        df_add = Table.read(FitsFiles[i_Files])
        df_add = FilterTimeErr(df_add)
        df_add = FilterSAPQual(df_add)
        df_add = FilterPDCErr(df_add)
        if len(df_add) > 0:
            df_add['PDCSAP_FLUX']= df_add['PDCSAP_FLUX']/np.median(df_add['PDCSAP_FLUX'])#/Average(df_add['PDCSAP_FLUX'])
            print('Read in the light curve data:')
            if (df_add[1]['TIME']>df[len(df) - 1]['TIME']):
                print(FitsFiles[i_Files])             
                df_trans = len(df)
                df = vstack([df, df_add], join_type='inner')
                # make dataset transtion smooth
                for idx in range (df_trans-5,df_trans+5):
                    df[idx]['PDCSAP_FLUX'] = Average([df[idx-5]['PDCSAP_FLUX'],df[idx+5]['PDCSAP_FLUX']])
            
    # Get and filter confirmed planets
    KIDplanets = confirmedplanets[np.isclose(confirmedplanets['kepid'], float(Keplerid), rtol = 0.0)]
    for i_remove in range (0,len(KIDplanets)):
        t0_low = KIDplanets[i_remove]['t_0'] - KIDplanets[i_remove]['T'] * int((KIDplanets[i_remove]['t_0']-df[1]['TIME'])/KIDplanets[i_remove]['T'])
        print('Remove the transits from the lightcurve with T = ' +str(KIDplanets[i_remove]['T'])+ ' and t0 = ' +str(t0_low))
        for i_transit in range (0, (int((df[len(df) - 1]['TIME']-t0_low)/KIDplanets[i_remove]['T'])+1)):
            idx_0 = getClosestNeighbour((t0_low+KIDplanets[i_remove]['T']*i_transit-1.5*KIDplanets[i_remove]['delta_t']/24.0), 'TIME', df)
            idx_1 = getClosestNeighbour((t0_low+KIDplanets[i_remove]['T']*i_transit+1.5*KIDplanets[i_remove]['delta_t']/24.0), 'TIME', df)                
            try:
                PDC_patch=Average([Average(df[(max(idx_0-5,0)):(idx_0+1)]['PDCSAP_FLUX']),Average(df[(idx_1):(min(idx_1+6,len(df)))]['PDCSAP_FLUX'])])
                mask = np.logical_not(np.logical_or((df['TIME'] < (t0_low+KIDplanets[i_remove]['T']*i_transit-1.5*KIDplanets[i_remove]['delta_t']/24.0 )),(df['TIME'] > (t0_low+KIDplanets[i_remove]['T']*i_transit+1.5*KIDplanets[i_remove]['delta_t']/24.0 ))))
                np.place(df['PDCSAP_FLUX'],mask,PDC_patch)     
            except ZeroDivisionError:
                mask = np.logical_or((df['TIME'] < (t0_low+KIDplanets[i_remove]['T']*i_transit-1.5*KIDplanets[i_remove]['delta_t']/24.0 )),(df['TIME'] > (t0_low+KIDplanets[i_remove]['T']*i_transit+1.5*KIDplanets[i_remove]['delta_t']/24.0 )))
                df = df[mask]
        df['PDCSAP_FLUX']= df['PDCSAP_FLUX']/np.median(df['PDCSAP_FLUX'])#/Average(df['PDCSAP_FLUX'])          
                
    if (Tremove != 0):
        print('Remove the transits from the lightcurve with T = ' +str(Tremove)+ ' and t0 = ' +str(t0remove))
        # Calculate first possible transit to be removed
        t0_low = t0remove - Tremove * int((t0remove-df[1]['TIME'])/Tremove)
        for i_transit in range (0, (int((df[len(df) - 1]['TIME']-t0_low)/Tremove)+1)):
            idx_0 = getClosestNeighbour((t0_low+Tremove*(i_transit-0.015)), 'TIME', df)
            idx_1 = getClosestNeighbour((t0_low+Tremove*(i_transit+0.015)), 'TIME', df)
            try:
                mask = np.logical_not(np.logical_or((df['TIME'] < (t0_low+Tremove*(i_transit-0.015) )),(df['TIME'] > (t0_low+Tremove*(i_transit+0.015) ))))
                PDC_patch=Average([Average(df[(max(idx_0-5,0)):(idx_0+1)]['PDCSAP_FLUX']),Average(df[(idx_1-1):(min(idx_1+5,len(df)))]['PDCSAP_FLUX'])])
                np.place(df['PDCSAP_FLUX'],mask,PDC_patch)  
            except ZeroDivisionError:
                mask = np.logical_or((df['TIME'] < (t0_low+Tremove*(i_transit-0.015) )),(df['TIME'] > (t0_low+Tremove*(i_transit+0.015) )))
                df = df[mask]			
        df['PDCSAP_FLUX']= df['PDCSAP_FLUX']/np.median(df['PDCSAP_FLUX'])#/Average(df['PDCSAP_FLUX'])
        
    
    # sort into groups regarding transitduration
    stepsize_sg = max(2, math.ceil(len(samples)/4))
    (samples).sort(order='tran_duration')
    samplegroup = [[],[],[],[]]
    group = 0
    for i_sg in range(0, len(samples), stepsize_sg ):
        samplegroup[group]=samples[i_sg:i_sg +stepsize_sg]
        (samplegroup[group]).sort(order='T')
        group = group +1
    
    for group in range (0, len(samplegroup)):   
        stepsize = max(2,math.ceil(len(samplegroup[group])/20))
        print('Split samples into groups of: ' + str(stepsize))
        for i in range(0, len(samplegroup[group]), stepsize ):
            if (samplegroup[group][i]['SDE'] < 7.0):# and (samplegroup[group][i]['T'] > 350.0 ):       
                transitfound = False
                
                flatten_lc, trend_lc = flatten(df['TIME'], df['PDCSAP_FLUX'], window_length=np.median(samplegroup[group][i:i +stepsize]['tran_duration'])*3, method='biweight', return_trend=True, robust=True)    
                print('The mean transit duration of the samples considered is: ' + str(np.median(samplegroup[group][i:i +stepsize]['tran_duration'])))
                    
                model = transitleastsquares(df['TIME'], flatten_lc)
                            
                
                results = model.power(period_max=math.ceil(max(samplegroup[group][i:i +stepsize]['T'])),
                                      period_min=int(min(samplegroup[group][i:i +stepsize]['T'])),
                                      R_star=R_star, 
                                      M_star=M_star, 
                                      duration=np.median(samplegroup[group][i:i +stepsize]['tran_duration']), 
                                      duration_range = 0.05, 
                                      duration_grid_step=0.005, 
                                      oversampling_factor=5,
                                      transit_depth_min=0.1*np.median(samplegroup[group][i:i +stepsize]['delta_F']) )        
                try:
                    
                    print('Period', format(results.period, '.5f'), 'd')
                    print(len(results.transit_times), 'transit times in time series:', \
                            ['{0:0.5f}'.format(i) for i in results.transit_times])
                    print('Transit depth', format(results.depth, '.5f'))
                    print('Best duration (days)', format(results.duration, '.5f'))
                    print('Signal detection efficiency (SDE):', results.SDE)
                        
                    for i2 in range(i,min(i +stepsize,len(samplegroup[group]))):    
                        # Check if TLS found the same transit
                        for i_transittime in results.transit_times:
                            if ((abs(i_transittime-samplegroup[group][i2]['t_0'])< (0.005* samplegroup[group][i2]['T'])) and
                                (abs(results.period-samplegroup[group][i2]['T'])< (0.001* samplegroup[group][i2]['T']))):
                                # Update with more accurate data from TLS and record SDE
                                samplegroup[group][i2]['T'] = results.period
                                samplegroup[group][i2]['t_0'] = i_transittime
                                # samplegroup[group][i2]['SDE'] = results.SDE
                                # samplegroup[group][i2]['SDE_Dis'] = abs(results.SDE/np.percentile(results.power, 99.99))
                                # samplegroup[group][i2]['tran_duration'] = max(results.duration,2/3*samplegroup[group][i2]['tran_duration'])
                                transitfound = True
                        
                                
                            
                        if (transitfound == True):
                            flatten_lc, trend_lc = flatten(df['TIME'], df['PDCSAP_FLUX'], window_length=samplegroup[group][i2]['tran_duration']*3, method='biweight', return_trend=True, robust=True)    
                                
                            for i_iter in range(3,10): 
                                transitfound = False
                                SDE_old = results.SDE
                            
                                model = transitleastsquares(df['TIME'], flatten_lc)
                                
                               
                                results = model.power(period_max=((1+i_iter**2/100)*max(samplegroup[group][i:i +stepsize]['T'])),period_min=((1-i_iter**2/100)*min(samplegroup[group][i:i +stepsize]['T'])),R_star=R_star, M_star=M_star, duration=samplegroup[group][i2]['tran_duration'], transit_depth_min=0.1*samplegroup[group][i2]['delta_F'], duration_grid_step=0.005, duration_range = 0.05,oversampling_factor=5 )        
                                try:
                                    
                                    print('Period', format(results.period, '.5f'), 'd')
                                    print(len(results.transit_times), 'transit times in time series:', \
                                            ['{0:0.5f}'.format(i) for i in results.transit_times])
                                    print('Transit depth', format(results.depth, '.5f'))
                                    print('Best duration (days)', format(results.duration, '.5f'))
                                    print('Signal detection efficiency (SDE):', results.SDE)
                                        
                                        
                                    # Check if TLS found the same transit
                                    for i_transittime in results.transit_times:
                                        if abs(i_transittime-samplegroup[group][i2]['t_0'])< (0.001* samplegroup[group][i2]['T']):
                                            # Update with more accurate data from TLS and record SDE
                                            samplegroup[group][i2]['T'] = results.period
                                            samplegroup[group][i2]['t_0'] = i_transittime
                                            samplegroup[group][i2]['SDE'] = results.SDE
                                            samplegroup[group][i2]['SDE_Dis'] = abs(results.SDE/np.percentile(results.power, 99.99))
                                            transitfound = True
                                        

                                        
                                    
                                except TypeError:
                                    print('No result in biweight filtered run') 
                                if (transitfound == False) or (results.SDE*1.5 < SDE_old):                        
                                    break
             
                except TypeError:
                    print('No result') 
    
    samples = samplegroup[0]
    for group in range (1, len(samplegroup)):   
        samples =np.concatenate([samples,samplegroup[group]])
    
    
    print('Update csv file')
    np.savetxt(path_to_csv+'Solutions_'+str(Keplerid)+ '.csv', samples,header='T;tran_duration;R_p;r_p;delta_F;t_0;Amount Transits;SDE;SDE_Dis', delimiter=";")


    print('Remove downloaded files')
    for i_Files in range (0,(len(FitsFiles))):
        os.remove(FitsFiles[i_Files])
        
    return 
                         


if __name__ == "__main__":
    main()

