# -*- coding: utf-8 -*-
"""
This program scans light curve files in FITS format to identify 
possible planet transits and calculate the planet data such as planet 
orbit radius and period by assuming a circular orbit.

DATE        AUTHOR          Modification
23/02/2021  Marcel Winzen   Initial version from PlanetScan
20/05/2021  Marcel Winzen   Save Results into specified folder
03/10/2021  Marcel Winzen   f_1 is dependent on percentile instead of error
"""

import os
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits as iof
from astropy.io.votable import parse_single_table
import astropy.constants as co
import FirstTransitDeterminator as FTD
from DataUseTools import * 
from ManageFits import GetListofFitsAvail, DownloadFits
import progressbar
from TLS_executer import doTLS

def main(FileFormat, min_kid, max_f_1):
    """ Main program """
    
    kepler_cat = parse_single_table("nph-nstedAPI.xml").to_table()
    kepler_cat['kepid'] = [int(i) for i in kepler_cat['kepid']]
    
    kepler_cat.sort('kepid')
    
    # import 
    confirmedplanets = np.genfromtxt('confirmedplanetsAPI.csv', delimiter=';',names=True,dtype=float)
    confirmedplanets = np.sort(confirmedplanets, order='kepid')
    KIDStart = int(confirmedplanets['kepid'][0])
    KIDEnd = int(confirmedplanets['kepid'][len(confirmedplanets)-1])
    
    
    index_start = np.searchsorted(kepler_cat['kepid'], KIDStart, side="left")
    index_end = np.searchsorted(kepler_cat['kepid'], KIDEnd, side="left")
    
    
    
    for index_data in range(index_start, index_end, 5):
        nextindex = False
        Keplerid = kepler_cat['kepid'][index_data]
        if (Keplerid < min_kid):
            nextindex = True
        if not np.any(np.isclose(confirmedplanets['kepid'], float(Keplerid), 0, 0.5)):
            print('No Planet Scan performed, since no confirmed planets available for KeplerID: ' +str(Keplerid) )
            nextindex = True         
        if (nextindex == False):
            FilesAvailable = DownloadFits(Keplerid,FileFormat)
            print('Download data with KeplerID: ' + str(Keplerid))
            FitsFiles = GetListofFitsAvail(Keplerid)
                        
            print('Read in the light curve data and the star radius from the Fits file given:')
            Fitsfile = FitsFiles[0]
            df = Table.read(Fitsfile)   
            df = FilterTimeErr(df)
            df = FilterSAPQual(df)
            df = FilterPDCErr(df)
            df['PDCSAP_FLUX']= df['PDCSAP_FLUX']/np.median(df['PDCSAP_FLUX'])
            
            
            R_star =  iof.getval(Fitsfile, 'RADIUS')
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
                
            
            R_star =  iof.getval(Fitsfile, 'RADIUS')
            Mass_lst = kepler_cat['mass'][(index_data):(index_data+5)]
            M_star =  Average(Mass_lst[np.isfinite(Mass_lst)])
            
   
        
        
            print('Kepler-ID: ' + str(Keplerid))
            print('Stellar radius: ' + str(iof.getval(Fitsfile, 'RADIUS')))
            
            # Get and filter confirmed planets
            KIDplanets = confirmedplanets[np.isclose(confirmedplanets['kepid'], float(Keplerid), rtol = 0.0)]
            for i_remove in range (0,len(KIDplanets)):
                t0_low = KIDplanets[i_remove]['t_0'] - KIDplanets[i_remove]['T'] * int((KIDplanets[i_remove]['t_0']-df[1]['TIME'])/KIDplanets[i_remove]['T'])
                print('Remove the transits from the lightcurve with T = ' +str(KIDplanets[i_remove]['T'])+ ' and t0 = ' +str(t0_low))
                for i_transit in range (0, (int((df[len(df) - 1]['TIME']-t0_low)/KIDplanets[i_remove]['T'])+1)):
                    idx_0 = getClosestNeighbour((t0_low+KIDplanets[i_remove]['T']*i_transit-1.0*KIDplanets[i_remove]['delta_t']/24.0), 'TIME', df)
                    idx_1 = getClosestNeighbour((t0_low+KIDplanets[i_remove]['T']*i_transit+1.0*KIDplanets[i_remove]['delta_t']/24.0), 'TIME', df)   
                                 
                    mask = np.logical_or((df['TIME'] < df[max(idx_0,0)]['TIME']),(df['TIME'] > df[min(idx_1,len(df)-1)]['TIME']))
                    df = df[mask]
                    
            
            print('Prepare filter parameters for dataset')
            last_index = len(df) - 1 
            t_SMA_min = 0.06 #2 * abs(df[last_index]['TIME'] - df[0]['TIME'])/ last_index
            df = FilterSAPQual(df)
            df = FilterPDCErr(df)
            y = np.array(df['PDCSAP_FLUX'])
            #if min_f_1 > 0.995:
            #    min_f_1 = 0.995
            # request last index after SAP Quality filter
            Solutions = Table([[],[],[],[],[],[],[],[],[],[]], names=('T', 'tran_duration', 'delta_F', 't_0', 't_SMA', 'f_2','f_1*c_SMA','Amount Transits','SDE','SDE_Dis'))
                                                                                                                                                                            
            for percentilefac in range(0, 7):
                last_index = len(df) - 1 
                if last_index > 1:
                    f_1 = max(1-np.percentile(df['PDCSAP_FLUX'],0.5*2**percentilefac), (np.percentile(df['PDCSAP_FLUX'], (100-0.5*2**percentilefac))-np.percentile(df['PDCSAP_FLUX'], 0.5*2**percentilefac))/2)
                        
                    
                    #return df
                
                    print('Run Test for the following parameters:')
                    print('t_SMA = 0.1; 0.08; 0.06 days')
                    print('f_1 = ' +str(f_1)  )
                    print('f_2 = 2' )
                    print('Number of data indices = ' +str(last_index))
                            
                try:
                    os.mkdir('Res_'+str(Keplerid))
                except OSError as exc:
                    pass
    
                for t_SMA_cnt in range(0, 3):
                    if (np.isfinite(M_star) == False) or (np.isfinite(R_star) == False) or (R_star > 3.5):
                        print('No Planet Scan performed, since no star radius or mass available for KeplerID: ' +str(Keplerid) )
                        print('Remove downloaded files')
                        for i_Files in range (0,(len(FitsFiles))):
                            os.remove(FitsFiles[i_Files])            
                        break                                                                                                     
                    t_SMA = 0.1 - 0.02* t_SMA_cnt  
                    if (t_SMA > 0.09):
                        c_SMA =  0.2
                    elif (t_SMA > 0.07) and (t_SMA < 0.09):
                        c_SMA =  0.25
                    elif (t_SMA > 0.05) and (t_SMA < 0.07):      
                        c_SMA =  0.33
                    
                    if (last_index>1):
                        for f_2 in range(2, 1, -1):
                            if f_1 < max_f_1:  # It can be assumed that only small planets can be discovered
                                c1 = max(getClosestNeighbour((df[0]['TIME']+2*t_SMA), 'TIME', df), 5)
                                print('Find planet transit solution with f_1 = ' + str(f_1)+' , t_SMA = '+ str(t_SMA) )
                                file1 = open("LastKeplerID.txt","w+") 
                                file1.write('KeplerID: '+str(Keplerid)+' f_1 = ' + str(f_1)+' , t_SMA = '+ str(t_SMA) )    
                                file1.close()    
                                bar = progressbar.ProgressBar(maxval=(last_index), \
                                widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
                                bar.start()
                                for index in range (2, last_index):
                                    if index > c1:
                                        Continue, T, tran_duration, delta_F, c0, c3, c1 = FTD.getfirstTransitData(index, t_SMA, f_2, f_1, R_star, M_star, df)
                                        if Continue == True:
                                            t_0 =  (df[c0]['TIME']+df[c3]['TIME'])/2*86400 
                                            print('First Transit time: ' + str(t_0/86400))
                                            print('Planet orbit eclipse time [days]: ' + str(tran_duration/86400))
                                            #print('Planet orbit period [days]: ' + str( T/ 86400)) 
                                            Solutions.add_row([T/86400, tran_duration/86400, delta_F, t_0/86400, t_SMA, f_2, f_1*c_SMA, 1, 0.0, 0.0])
                                            Solutions.write('Res_'+str(Keplerid)+'/Solutions_'+str(Keplerid)+ '.html', format='ascii.html', overwrite = True)    
                                            
                                    bar.update(index+1)
                                bar.finish()
    
                
                
                
                                    
        
        
                if (len(Solutions) > 250) or ((len(Solutions) > 0)and (percentilefac==6)):
                    t_start = df[0]['TIME']
                    t_end = df[last_index]['TIME']
                    t_data = t_end - t_start      
                                
                    Solutions = pairTransits(Solutions, t_data)
                    Solutions.write('Res_'+str(Keplerid)+'/Solutions_'+str(Keplerid)+ '.csv', format='csv', overwrite = True, delimiter = ';')    
                            
                    if (len(Solutions) > 0) :                
                        # Sort out repeating solutions and such that are due to missed intermediate eclipse
                        Solution, DoTLS = removesimilarSol(Solutions)
                        Solution.write('Res_'+str(Keplerid)+'/Solutions_'+str(Keplerid)+ '.csv', format='csv', overwrite = True, delimiter = ';')    
                
                        if DoTLS == True:
            
                            Solution = processSol(Solution, R_star, M_star)
                
                            Solution.write('Res_'+str(Keplerid)+'/Solutions_'+str(Keplerid)+ '.csv', format='csv', overwrite = True, delimiter = ';')    
                
                            doTLS(Keplerid, 'Res_'+str(Keplerid)+'/', 0, 0,  M_star, R_star, confirmedplanets)
                            
                        
                    break                            

            
    return 
                         


if __name__ == "__main__":
    main()

