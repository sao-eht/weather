import os
import numpy as np
from amread import reader
import astropy.table as tab
import matplotlib.pyplot as plt
import datetime as dt

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
    return a, b, c

scope_data = [['JCMT', 19.8228, -155.477, 4092], ['GLT', 76.5 , -68.7 ,84], ['LMT', 18.9858, -97.3147, 4640]]
#scope_data = [['GLT', 76.5 , -68.7 ,84]]
scope_len = len(scope_data)


hours = np.arange(0,24,1)
hour_len = len(hours)
#######################################################################################################################################


GLT    = tab.Table([[] ,[]   ,[]  ,[]  ,[]  ], names = ('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('GLT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
JCMT   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  ], names = ('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('JCMT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
#PICO = tab.Table([[] ,[]   ,[]  ,[]  ,[]  ], names = ('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('PICO')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
LMT    = tab.Table([[] ,[]   ,[]  ,[]  ,[]  ], names = ('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('LMJCMTT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
#ALMA   = tab.Table([[] ,[]   ,[]  ,[]  ,[]  ], names = ('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('ALMA')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))

def sites(date, ut, prod):
    SITES = np.chararray(scope_len, 6)
    freqs = np.zeros((scope_len))
    taus = np.zeros((scope_len))
    trans = np.zeros((scope_len))

    for i in np.arange(scope_len):
        f, t, tran    = gfs2file(scope_data[i][1], scope_data[i][2], scope_data[i][3], date, ut, str(prod), str(date)+str(ut)+str(scope_data[i][0])+'layers.amc', str(date)+str(ut)+str(scope_data[i][0])+'plugged.amc', str(date)+str(ut)+str(scope_data[i][0])+'.out')
        freqs[i]= f
        taus[i] = t
        trans[i] = tran
        SITES[i] = scope_data[i][0]

    print 'Atmospheric Model Tasks Complete for ' + str(date)

    dates = np.array(np.zeros(scope_len))

    for i in np.arange(scope_len):
        dates[i] = date

    t = tab.Table([dates, SITES , freqs, taus, trans], names = ('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str(date)})

    for i in np.arange(scope_len):
        index = np.where(t['Site'] == str(scope_data[i][0]))
        eval(t['Site'][i]).add_row(t[index][0])


    return t

GLT_hourly = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('GLT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
JCMT_hourly = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('JCMT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
#PICO_hourly = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('PICO')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
LMT_hourly = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('LMT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
#ALMA_hourly = tab.Table([[] ,[]   ,[]  ,[]  ,[]  , []], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('ALMA')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))


def site_hour(date, ut, prod, hour):
    SITES = np.chararray(scope_len, 6)
    freqs = np.zeros((scope_len))
    taus = np.zeros((scope_len))
    trans = np.zeros((scope_len))
    
    for i in np.arange(scope_len):
        f, t, tran    = gfs2file(scope_data[i][1], scope_data[i][2], scope_data[i][3], date, ut, str(prod), str(date)+str(ut)+str(scope_data[i][0])+'layers.amc', str(date)+str(ut)+str(scope_data[i][0])+'plugged.amc', str(date)+str(ut)+str(scope_data[i][0])+'.out')
        freqs[i]= f
        taus[i] = t
        trans[i] = tran
        SITES[i] = scope_data[i][0]
    if taus[0] != 0.0:
        print 'Atmospheric Model Tasks Complete for ' + str(date) + str(hour)

        dates = np.array(np.zeros(len(freqs)))
        hours = np.array(np.zeros(len(freqs)))

        for i in np.arange(len(freqs)):
            dates[i] = date
            hours[i] = hour

        t = tab.Table([dates, hours, SITES , freqs, taus, trans], names = ('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str(date)})

        for i in np.arange(scope_len):
            index = np.where(t['Site'] == str(scope_data[i][0]))
            eval(t['Site'][i]+'_hourly').add_row(t[index][0])

##############
##############
#EDITING AREA#
##############
##############



days = [str(20180329)] ####day in YYYYMMDD from which the forecasting begins, does not go back farther than 14 days prior to the present day at which the reader is editing this script

day_len = len(days)
b = 100 #Hours of forecasting to explore
c = 12#UT hour from which to start the forecasting
a = np.arange(1,b)

#######Un-commenting the following portion of the script initializes the data search for the forecast data on the NOMADS database, if the above search has already been previously done, then it can be skipped if the data tables have been saved locally in the same directory as this script

for i in np.arange(len(a)):
    if a[i] < 10:
        site_hour(int((days[0])), c, 'f00' + str(a[i]),  c+int(a[i]))
    if (a[i] < 100) & (a[i] >= 10):
        site_hour(int((days[0])), c, 'f0' + str(a[i]),  c+int(a[i]))
    if a[i] >= 100:
        site_hour(int((days[0])), c, 'f' + str(a[i]),  c+int(a[i]))
    print 'Hour ' + str(a[i]) + ' completed'

GLT_hourly.write('GLT_forecast_' + str(days[0])+ '_'+str(b+c)+'h_'+'.dat', format='ascii')
JCMT_hourly.write('JCMT_forecast_' + str(days[0]) + '_'+str(b+c)+'h_'+'.dat', format='ascii')
#PICO_hourly.write('PICO_forecast_' + str(days[0]) + '_'+str(b+c)+'h_'+'.dat', format='ascii')
LMT_hourly.write('LMT_forecast_' + str(days[0])  + '_'+str(b+c)+'h_'+'.dat', format='ascii')
#ALMA_hourly.write('ALMA_forecast_' + str(days[0]) + '_'+str(b+c)+'h_'+'.dat', format='ascii')

##Reads in data which has already been calculated
GLT_forecast = tab.Table.read('GLT_forecast_' + str(days[0])+ '_'+str(b+c)+'h_'+'.dat', format='ascii')
JCMT_forecast = tab.Table.read('JCMT_forecast_' + str(days[0])+ '_'+str(b+c)+'h_'+'.dat', format='ascii')
#PICO_forecast= tab.Table.read('PICO_forecast_' + str(days[0])+ '_'+str(b+c)+'h_'+'.dat', format='ascii')
LMT_forecast= tab.Table.read('LMT_forecast_' + str(days[0])+ '_'+str(b+c)+'h_'+'.dat', format='ascii')
#ALMA_forecast= tab.Table.read('ALMA_forecast_' + str(days[0])+ '_'+str(b+c)+'h_'+'.dat', format='ascii')

#Plotting routines
for i in np.arange(scope_len):
    x = eval(str(scope_data[i][0])+'_forecast')['Hour UT']
    y = eval(str(scope_data[i][0])+'_forecast')['Opacity']
    z = eval(str(scope_data[i][0])+'_forecast')['Date [YYYYMMDD]'][0]
    plt.figure(i,figsize=(6,4))
    plt.ylim(0,1)
    plt.plot(x , y)
    plt.xlabel('Hour [UT]')
    plt.ylabel(r"Opacity ($\tau$ at 221.1 GHz)")
    plt.title('Forecast Opacity for the '+ str(scope_data[i][0]) + ' for ' + str(z))

plt.figure(26, figsize=(7,6))
plt.ylim(0,1)
for i in np.arange(scope_len):
    x = eval(str(scope_data[i][0])+'_forecast')['Hour UT']
    y = eval(str(scope_data[i][0])+'_forecast')['Opacity']
    z = eval(str(scope_data[i][0])+'_forecast')['Date [YYYYMMDD]'][0]
    plt.plot(x , y, label = str(scope_data[i][0]))

plt.legend()
plt.xlabel('Hour [UT]')
plt.ylabel(r"Opacity ($\tau$ at 221.1 GHz)")
plt.title('Forecast Opacity for the EHT DR for ' + str(z))
plt.show()



