##########
#Started by Rodrigo Cordova on April 4, 2018
#Usage: Creation of an opacity file for the EHT for a set of given sites, times, etc. 

#packages


import os
import numpy as np
from amread import reader
import astropy.table as tab
import datetime as dt
import argparse
import dateutil.parser as dparser
import matplotlib.pyplot as plt
from math import floor

###function which reads in the data from am ion
def gfs2file(lat, lon, h, date, ut, prod, datafile, newfile, outfile):
    f = open(str(datafile),'w')
    os.system( str('python gfs2am_v3.py ') + str( lat) +' ' + str( lon) +' '+ str( h) +' '+ str( date) +' '+ str( ut) +' ' + str( prod) +' ' + '>' + str(datafile))

    f = open(str(datafile),'r')
    newf = open(str(newfile),'w')
    lines = f.readlines() # read old content

    newf.write('f 221.1 GHz  221.1 GHz  50 MHz \n') # write new content at the beginning
    newf.writelines(['output f GHz  tau neper tx none \n'])
    newf.write('T0 2.7 K')

    for line in lines: # write old content after new
        newf.write(line)
    newf.close()
    f.close()


    os.system('am' + ' ' + str(newfile) + '>' + str(outfile))
    a,b,c =  reader(str(outfile), 'Tau', 'Transmittance')
    os.system('rm ' +str(newfile) )
    os.system('rm ' +str(datafile) )
    os.system('rm ' +str(outfile) )
    return a, b, c




##all telesocopes

scope_data = [['SMT',32.7016, -109.891, 3185],['LMT', 18.9858, -97.3147, 4640],['SMA', 19.8242, -155.478, 4080],['PV', 37.0661, -3.3925, 2850], ['SPT', -90, 0, 2800],['ALMA', -23.0193, -67.7532, 5058.7],['PDB', 44.6339, 5.9081, 2550],['GLT', 76.5 , -68.7 ,84]]

scope_len = len(scope_data)


SMT   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('SMT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
LMT   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('LMT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
SMA   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('SMA')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
PV    = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('PV')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
SPT   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('SPT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
ALMA  = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('ALMA')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
PDB   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('PDB')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
GLT = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('GLT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))


###Design selection tool for which telescopes you want to measure for
for i in np.arange(scope_len):
    print str(i) + "= " + scope_data[i][0]

print "the first argument are the indexes of telescopes we want to run the forecast tool for separated by commas, for the numbers listed above for all the array sites listed. To add more, go into code and add in the antenna name, latitude, longitude, and height."

Instructions = '''

Definitions of all the parameters of the forecasting functions

telescopes      array selection, integer numbers with commas in
                 between them, all within a string

start_date       day on which the forecast begins, YYYYMMDD
                 structure, must be a day within the last 14
                 days

start_hour       UT hour at which the forecast begins

hours_future     how many hours into the future we look
                 into (000h - 384h)
                 
prior            type yes or no directing if this calculation
                 has been run previously and the relevant files
                 are already in your directory, if they are, type
                 yes, if not, type no, selecting yes will trigger
                 the plotting routines, must be a string"
                 
    Example:
    forecasting('1,5,6', 20180411, 6, 3, 'yes')
        calls for 3 hours of forecast for the LMT, ALMA, and PDB
        for April 11, 2018, at UT 6, and indicates there has
        already been a prior run of forecast and would now want
        to plot
        
        
                '''
print Instructions

def forecasting(telescopes, start_date, start_hour,hours_future, prior):
    ##############
    telescope_list = []
    for i in np.arange(len(telescopes)):
        try:
            telescope_list.append(int(telescopes[i]))
        except ValueError:
            continue

    forecast_ants = []
    for i in np.arange(len(telescope_list)):
        forecast_ants.append(scope_data[telescope_list[i]])
    print forecast_ants

    ant_len = len(forecast_ants)
    ##############

    #### Function for the creation of forecast file structure/ calling in the correct time/date/location
    def site_hour(date, ut, prod, hour, length, forecast_list):
        SITES = np.chararray(length, 6)
        freqs = np.zeros((length))
        taus = np.zeros((length))
        trans = np.zeros((length))
    
        for i in np.arange(length):
            f, t, tran    = gfs2file(forecast_list[i][1], forecast_list[i][2], forecast_list[i][3], date, ut, str(prod), str(date)+str(ut)+str(forecast_list[i][0])+'layers.amc', str(date)+str(ut)+str(forecast_list[i][0])+'plugged.amc', str(date)+str(ut)+str(forecast_list[i][0])+'.out')
            freqs[i]= f
            taus[i] = t
            trans[i] = tran
            SITES[i] = forecast_list[i][0]
        if taus[0] != 0.0:
            print 'Atmospheric Model Tasks Complete for ' + str(date) + str(hour)

            dates = np.array(np.zeros(len(freqs)))
            hours = np.array(np.zeros(len(freqs)))
            
            for i in np.arange(len(freqs)):
                dates[i] = date + floor(int(prod[1:])/24.0)
                hours[i] = hour - (24* floor(int(prod[1:])/24.0))

            for i in np.arange(len(hours)):
                if hours[i] >= 24:
                    hours[i] = hours[i] - (24* floor(int(hours[i])/24.0))
                    dates[i] = dates[i] + floor(int(hours[i])/24.0)
                else:
                    continue

            t = tab.Table([dates, hours, SITES , freqs, taus, trans], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str(date)})
            for i in np.arange(length):
                index = np.where(t['Site'] == str(forecast_list[i][0]))
                eval(t['Site'][i]).add_row(t[index][0])
    ############## initiating parameters
    days = [str(start_date)]

    day_len = len(days)
    b = int(hours_future) #Hours of forecasting to explore
    c = int(start_hour)#UT hour from which to start the forecasting
    a = np.arange(1,b)

    if prior == 'no':
        for i in np.arange(len(a)):
            if a[i] < 10:
                site_hour(int((days[0])), c, 'f00' + str(a[i]),  c+int(a[i]), ant_len, forecast_ants)
            if (a[i] < 100) & (a[i] >= 10):
                site_hour(int((days[0])), c, 'f0' + str(a[i]),  c+int(a[i]), ant_len, forecast_ants)
            if a[i] >= 100:
                site_hour(int((days[0])), c, 'f' + str(a[i]),  c+int(a[i]), ant_len, forecast_ants)
            print 'Hour ' + str(a[i]) + ' completed'
        for i in np.arange(len(forecast_ants)):
            eval(forecast_ants[i][0]).write(forecast_ants[i][0]+'_forecast_' + str(days[0])+ '_'+str(b+c)+'hrs'+'.dat', format='ascii')


    if prior == 'yes': ####this continues into the plotting routines, the files have already been created
        forecasts = []
        for i in np.arange(len(forecast_ants)):
            forecast = tab.Table.read(forecast_ants[i][0]+'_forecast_' + str(days[0])+ '_'+str(b+c)+'hrs'+'.dat', format='ascii')

            y = forecast['Opacity']
            x = forecast['Date [YYYYMMDD]'] + (forecast['Hour UT']/24.0)
            plt.figure(i,figsize=(6,4))
            plt.ylim(0,1)
            plt.plot(x , y)
            plt.xlabel('Date [YYYYMMDD.Hour]')
            plt.ylabel(r"Opacity ($\tau$ at 221.1 GHz)")
            plt.title('Forecast Opacity for the '+ str(forecast_ants[i][0]) + ' for ' + str(forecast['Date [YYYYMMDD]'][0]))
        
            
            forecasts.append(forecast)
        
        plt.figure(21, figsize=(8,6))
        for i in np.arange(len(forecasts)):
            y = forecasts[i]['Opacity']
            x = forecasts[i]['Date [YYYYMMDD]'] + (forecasts[i]['Hour UT']/24.0)
            plt.plot(x , y, label = str(forecast_ants[i][0]))
            plt.xlabel('Date [YYYYMMDD.Hour]')
            plt.ylabel(r"Opacity ($\tau$ at 221.1 GHz)")
            plt.title('Forecast Opacity for the EHT array for ' + str(forecasts[i]['Date [YYYYMMDD]'][0]))
        plt.legend()

        plt.show()
        return forecasts


'''
telescope_list = []
for i in np.arange(len(telescopes)):
    try:
        telescope_list.append(int(telescopes[i]))
    except ValueError:
        continue
print telescope_list

forecast_ants = []
for i in np.arange(len(telescope_list)):
    forecast_ants.append(scope_data[telescope_list[i]])
ant_len = len(forecast_ants)
print forecast_ants

for i in np.arange(len(forecast_ants)):
    eval(forecast_ants[i][0]).write(forecast_ants[i][0]+'_forecast_' + str(days[0])+ '_'+str(b+c)+'hrs'+'.dat', format='ascii')

    '''
