import os
import numpy as np
from amread import reader
import astropy.table as tab
import matplotlib.pyplot as plt
import datetime as dt

ddirec = 'data/'
adirec = 'data/amc/'
odirec = 'data/out/'

if (os.path.isdir(ddirec) == True):
	print('Output directory exists. Aborting...')
	exit()
os.system('mkdir '+ddirec)
os.system('mkdir '+adirec)
os.system('mkdir '+odirec)

def gfs2file(lat, lon, h, date, ut, prod, datafile, newfile, outfile):
    os.system( str('python local_grb_gfs2am.py ') + str( lat) +' ' + str( lon) +' '+ str( h) +' '+ str( date) +' '+ str( ut) +' ' + str( prod) +' ' + '>' + adirec + str(datafile))

    f = open(adirec + str(datafile),'r')
    newf = open(adirec + str(newfile),'w')
    lines = f.readlines() # read old content

    newf.write('f 221.1 GHz  221.1 GHz  50 MHz \n') # write new content at the beginning
    newf.writelines(['output f GHz  tau neper tx none \n'])
    newf.write('T0 2.7 K')

    for line in lines: # write old content after new
        newf.write(line)
    newf.close()
    f.close()


    os.system('am' + ' ' + adirec + str(newfile) + '>' + odirec + str(outfile))
    a,b,c =  reader(odirec + str(outfile), 'Tau', 'Transmittance')

    return a, b, c

#######################################################################################################################################
#scope_data = [['SMT',32.7016, -109.891, 3185],['LMT', 18.9858, -97.3147, 4640],['JCMT', 19.8228, -155.477, 4092],['SMA', 19.8242, -155.478, 4080],['Iram30', 37.0661, -3.3925, 2850], ['SPT', -90, 0, 2800],['ALMA', -23.0193, -67.7532, 5058.7],['APEX', -23.0058, -67.7592, 5100]]
scope_data = [['SMT',32.7016, -109.891, 3185],['LMT', 18.9858, -97.3147, 4640],['SMA', 19.8242, -155.478, 4080],['PV', 37.0661, -3.3925, 2850], ['SPT', -90, 0, 2800],['ALMA', -23.0193, -67.7532, 5058.7],['PDB', 44.6339, 5.9081, 2550]]

scope_len = len(scope_data)


#days = np.arange(20070201, 20070232, 1)
#days = np.append(days, np.arange(20070301, 20070332, 1))
days = np.arange(20070301, 20070332, 1)
days = np.append(days, np.arange(20070401, 20070432, 1))
days = np.append(days, np.arange(20070501, 20070532, 1))
#days = np.append(days, np.arange(20070601, 20070632, 1))
#days = np.append(days, np.arange(20070701, 20070732, 1))
#days = np.append(days, np.arange(20070801, 20070832, 1))
#days = np.append(days, np.arange(20070901, 20070932, 1))
#days = np.append(days, np.arange(20071001, 20071032, 1))
#days = np.append(days, np.arange(20071101, 20071132, 1))
#days = np.append(days, np.arange(20071201, 20071232, 1))

#days = np.append(days, np.arange(20080101, 20080132, 1))
#days = np.append(days, np.arange(20080201, 20080232, 1))
days = np.append(days, np.arange(20080301, 20080332, 1))
days = np.append(days, np.arange(20080401, 20080432, 1))
days = np.append(days, np.arange(20080501, 20080532, 1))
#days = np.append(days, np.arange(20080601, 20080632, 1))
#days = np.append(days, np.arange(20080701, 20080732, 1))
#days = np.append(days, np.arange(20080801, 20080832, 1))
#days = np.append(days, np.arange(20080901, 20080932, 1))
#days = np.append(days, np.arange(20081001, 20081032, 1))
#days = np.append(days, np.arange(20081101, 20081132, 1))
#days = np.append(days, np.arange(20081201, 20081232, 1))

#days = np.append(days, np.arange(20090101, 20090132, 1))
#days = np.append(days, np.arange(20090201, 20090232, 1))
days = np.append(days, np.arange(20090301, 20090332, 1))
days = np.append(days, np.arange(20090401, 20090432, 1))
days = np.append(days, np.arange(20090501, 20090532, 1))
#days = np.append(days, np.arange(20090601, 20090632, 1))
#days = np.append(days, np.arange(20090701, 20090732, 1))
#days = np.append(days, np.arange(20090801, 20090832, 1))
#days = np.append(days, np.arange(20090901, 20090932, 1))
#days = np.append(days, np.arange(20091001, 20091032, 1))
#days = np.append(days, np.arange(20091101, 20091132, 1))
#days = np.append(days, np.arange(20091201, 20091232, 1))

#days = np.append(days, np.arange(20100101, 20100132, 1))
#days = np.append(days, np.arange(20100201, 20100232, 1))
days = np.append(days, np.arange(20100301, 20100332, 1))
days = np.append(days, np.arange(20100401, 20100432, 1))
days = np.append(days, np.arange(20100501, 20100532, 1))
#days = np.append(days, np.arange(20100601, 20100632, 1))
#days = np.append(days, np.arange(20100701, 20100732, 1))
#days = np.append(days, np.arange(20100801, 20100832, 1))
#days = np.append(days, np.arange(20100901, 20100932, 1))
#days = np.append(days, np.arange(20101001, 20101032, 1))
#days = np.append(days, np.arange(20101101, 20101132, 1))
#days = np.append(days, np.arange(20101201, 20101232, 1))

#days = np.append(days, np.arange(20110101, 20110132, 1))
#days = np.append(days, np.arange(20110201, 20110232, 1))
days = np.append(days, np.arange(20110301, 20110332, 1))
days = np.append(days, np.arange(20110401, 20110432, 1))
days = np.append(days, np.arange(20110501, 20110532, 1))
#days = np.append(days, np.arange(20110601, 20110632, 1))
#days = np.append(days, np.arange(20110701, 20110732, 1))
#days = np.append(days, np.arange(20110801, 20110832, 1))
#days = np.append(days, np.arange(20110901, 20110932, 1))
#days = np.append(days, np.arange(20111001, 20111032, 1))
#days = np.append(days, np.arange(20111101, 20111132, 1))
#days = np.append(days, np.arange(20111201, 20111232, 1))

#days = np.append(days, np.arange(20120101, 20120132, 1))
#days = np.append(days, np.arange(20120201, 20120232, 1))
days = np.append(days, np.arange(20120301, 20120332, 1))
days = np.append(days, np.arange(20120401, 20120432, 1))
days = np.append(days, np.arange(20120501, 20120532, 1))
#days = np.append(days, np.arange(20120601, 20120632, 1))
#days = np.append(days, np.arange(20120701, 20120732, 1))
#days = np.append(days, np.arange(20120801, 20120832, 1))
#days = np.append(days, np.arange(20120901, 20120932, 1))
#days = np.append(days, np.arange(20121001, 20121032, 1))
#days = np.append(days, np.arange(20121101, 20121132, 1))
#days = np.append(days, np.arange(20121201, 20121232, 1))

#days = np.append(days, np.arange(20130101, 20130132, 1))
#days = np.append(days, np.arange(20130201, 20130232, 1))
days = np.append(days, np.arange(20130301, 20130332, 1))
days = np.append(days, np.arange(20130401, 20130432, 1))
days = np.append(days, np.arange(20130501, 20130532, 1))
#days = np.append(days, np.arange(20130601, 20130632, 1))
#days = np.append(days, np.arange(20130701, 20130732, 1))
#days = np.append(days, np.arange(20130801, 20130832, 1))
#days = np.append(days, np.arange(20130901, 20130932, 1))
#days = np.append(days, np.arange(20131001, 20131032, 1))
#days = np.append(days, np.arange(20131101, 20131132, 1))
#days = np.append(days, np.arange(20131201, 20131232, 1))

#days = np.append(days, np.arange(20140101, 20140132, 1))
#days = np.append(days, np.arange(20140201, 20140232, 1))
days = np.append(days, np.arange(20140301, 20140332, 1))
days = np.append(days, np.arange(20140401, 20140432, 1))
days = np.append(days, np.arange(20140501, 20140532, 1))
#days = np.append(days, np.arange(20140601, 20140632, 1))
#days = np.append(days, np.arange(20140701, 20140732, 1))
#days = np.append(days, np.arange(20140801, 20140832, 1))
#days = np.append(days, np.arange(20140901, 20140932, 1))
#days = np.append(days, np.arange(20141001, 20141032, 1))
#days = np.append(days, np.arange(20141101, 20141132, 1))
#days = np.append(days, np.arange(20141201, 20141232, 1))

#days = np.append(days, np.arange(20150101, 20150132, 1))
#days = np.append(days, np.arange(20150201, 20150232, 1))
days = np.append(days, np.arange(20150301, 20150332, 1))
days = np.append(days, np.arange(20150401, 20150432, 1))
days = np.append(days, np.arange(20150501, 20150532, 1))
#days = np.append(days, np.arange(20150601, 20150632, 1))
#days = np.append(days, np.arange(20150701, 20150732, 1))
#days = np.append(days, np.arange(20150801, 20150832, 1))
#days = np.append(days, np.arange(20150901, 20150932, 1))
#days = np.append(days, np.arange(20151001, 20151032, 1))
#days = np.append(days, np.arange(20151101, 20151132, 1))
#days = np.append(days, np.arange(20151201, 20151232, 1))

#days = np.append(days, np.arange(20160101, 20160132, 1))
#days = np.append(days, np.arange(20160201, 20160232, 1))
days = np.append(days, np.arange(20160301, 20160332, 1))
days = np.append(days, np.arange(20160401, 20160432, 1))
days = np.append(days, np.arange(20160501, 20160532, 1))
#days = np.append(days, np.arange(20160601, 20160632, 1))
#days = np.append(days, np.arange(20160701, 20160732, 1))
#days = np.append(days, np.arange(20160801, 20160832, 1))
#days = np.append(days, np.arange(20160901, 20160932, 1))
#days = np.append(days, np.arange(20161001, 20161032, 1))
#days = np.append(days, np.arange(20161101, 20161132, 1))
#days = np.append(days, np.arange(20161201, 20161232, 1))

#days = np.append(days, np.arange(20170101, 20170132, 1))
#days = np.append(days, np.arange(20170201, 20170232, 1))
days = np.append(days, np.arange(20170301, 20170332, 1))
days = np.append(days, np.arange(20170401, 20170432, 1))
days = np.append(days, np.arange(20170501, 20170532, 1))
#days = np.append(days, np.arange(20170601, 20170632, 1))
#days = np.append(days, np.arange(20170701, 20170732, 1))
#days = np.append(days, np.arange(20170801, 20170832, 1))
#days = np.append(days, np.arange(20170901, 20170932, 1))
#days = np.append(days, np.arange(20171001, 20171032, 1))
#days = np.append(days, np.arange(20171101, 20171132, 1))
#days = np.append(days, np.arange(20171201, 20171232, 1))



day_len = len(days)

hours = np.arange(0,24,1)
hour_len = len(hours)
#######################################################################################################################################


SMT   = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('SMT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
LMT   = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('LMT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
#JCMT  = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('JCMT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
SMA   = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('SMA')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
#Iram30= tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('Iram30')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
PV    = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('PV')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
SPT   = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('SPT')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
ALMA  = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('ALMA')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
#APEX  = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('APEX')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))
PDB   = tab.Table([[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('PDB')},dtype=('float64', 'str', 'float64', 'float64', 'float64'))


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

'''
for i in np.arange(len(days)):
    t = sites(days[i], 00, 'anl')

plt.figure(1)
for i in np.arange(scope_len):
    x = eval(str(scope_data[i][0])+'_hourly')['Date [YYYYMMDD]']
    y = eval(str(scope_data[i][0])+'_hourly')['Opacity']
    #ax.set_xticks([i for i in range(len(y))])
    plt.semilogy(y, label = str(scope_data[i][0]))
    
    
    
plt.ylabel('Optical Depth')
plt.title('2016')
plt.legend()

'''
SMT_hourly  = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('SMT') },dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
LMT_hourly  = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('LMT') },dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
#JCMT_hourly = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('JCMT')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
SMA_hourly  = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('SMA') },dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
SPT_hourly  = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('SPT') },dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
ALMA_hourly = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('ALMA')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
#APEX_hourly = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('APEX')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
#Iram30_hourly = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]', 'Hour UT', 'Site', 'Frequency [GHz]', 'Opacity', 'Transmittance'), meta={'name': str('Iram30')},dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
PV_hourly   = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('PV')  },dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))
PDB_hourly  = tab.Table([[],[],[],[],[],[]],names=('Date [YYYYMMDD]','Hour UT','Site','Frequency [GHz]','Opacity','Transmittance'),meta={'name':str('PDB') },dtype=('int', 'float64', 'str', 'float64', 'float64', 'float64'))


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

for i in np.arange(len(days)):
    site_hour(int((days[i])), 00, 'anl',  00)
    #site_hour(int((days[i])), 00, 'f001', 01)
    #site_hour(int((days[i])), 00, 'f002', 02)
    site_hour(int((days[i])), 00, 'f003', 03)
    #site_hour(int((days[i])), 00, 'f004', 04)
    #site_hour(int((days[i])), 00, 'f005', 05)
    site_hour(int((days[i])), 06, 'anl',  06)
    #site_hour(int((days[i])), 06, 'f001', 07)
    #site_hour(int((days[i])), 06, 'f002', 8)
    site_hour(int((days[i])), 06, 'f003', 9)
    #site_hour(int((days[i])), 06, 'f004', 10)
    #site_hour(int((days[i])), 06, 'f005', 11)
    site_hour(int((days[i])), 12, 'anl',  12)
    #site_hour(int((days[i])), 12, 'f001', 13)
    #site_hour(int((days[i])), 12, 'f002', 14)
    site_hour(int((days[i])), 12, 'f003', 15)
    #site_hour(int((days[i])), 12, 'f004', 16)
    #site_hour(int((days[i])), 12, 'f005', 17)
    site_hour(int((days[i])), 18, 'anl',  18)
    #site_hour(int((days[i])), 18, 'f001', 19)
    #site_hour(int((days[i])), 18, 'f002', 20)
    site_hour(int((days[i])), 18, 'f003', 21)
    #site_hour(int((days[i])), 18, 'f004', 22)
    #site_hour(int((days[i])), 18, 'f005', 23)
    print 'One Day Completed'



#dat.write('table.dat', format='ascii')

SMT_hourly.write(ddirec+'SMT_hourly_' + str(days[0]) +'-' + str(days[-1]) +'.dat', format='ascii')
LMT_hourly.write(ddirec+'LMT_hourly_' + str(days[0]) +'-' + str(days[-1]) +'.dat', format='ascii')
#JCMT_hourly.write(ddirec+'JCMT_hourly_' + str(days[0])+ '-' + str(days[-1]) +'.dat', format='ascii')
SMA_hourly.write(ddirec+'SMA_hourly_' + str(days[0]) +'-' + str(days[-1]) +'.dat', format='ascii')
#Iram30_hourly.write(ddirec+'Iram30_hourly_' + str(days[0])+ '-' + str(days[-1]) +'.dat', format='ascii')
PV_hourly.write(ddirec+'PV_hourly_' + str(days[0])+ '-' + str(days[-1]) +'.dat', format='ascii')
SPT_hourly.write(ddirec+'SPT_hourly_' + str(days[0]) +'-' + str(days[-1]) +'.dat', format='ascii')
ALMA_hourly.write(ddirec+'ALMA_hourly_' + str(days[0])+ '-' + str(days[-1]) +'.dat', format='ascii')
#APEX_hourly.write(ddirec+'APEX_hourly_' + str(days[0]) +'-' + str(days[-1]) +'.dat', format='ascii')
PDB_hourly.write(ddirec+'PDB_hourly_' + str(days[0])+ '-' + str(days[-1]) +'.dat', format='ascii')


'''
def indexing(site):
    std = np.zeros(hour_len)
    mean = np.zeros(hour_len)
    for i in np.arange(hour_len):
        index = np.where(eval(str(site)+'_hourly')['Hour UT'] == hours[i])
        indv_hour = eval(str(site)+'_hourly')[index]
        std[i] = np.std(indv_hour['Opacity'])
        mean[i] = np.mean(indv_hour['Opacity'])
        plt.scatter(indv_hour['Hour UT'],(indv_hour['Opacity']), s = 15)
    #plt.title(str(site)+' Hourly Measurements over January 1-31 for 2008 through 2017')
    plt.xlabel('Hour (UT)')
    plt.ylabel('Opacity')
    #plt.savefig(str(site)+'_Hourly_Measurements_over_January_1-31_for_2008_through_2017.png', dpi=300)
    plt.show()
    return std, mean

for i in np.arange(scope_len):
    error, mean = indexing(scope_data[i][0])
    plt.errorbar(hours, mean, yerr = error, fmt = 'o')
    #plt.title(str(scope_data[i][0])+' Hourly Average over January 1-31 for 2008 through 2017')
    plt.xlabel('Hour (UT)')
    plt.ylabel('Opacity')
    #plt.savefig(str(scope_data[i][0])+'_Hourly_Average_over_January_1-31_for_2008_through_2017.png', dpi=300)
    plt.show()


def indexing_stat(site, i_days):
    std = np.zeros(hour_len)
    mean = np.zeros(hour_len)
    for i in np.arange(hour_len):
        index = np.where(eval(str(site)+'_hourly')['Hour UT'] == hours[i])
        indv_hour = eval(str(site)+'_hourly')[index]
        day_ind = np.where(str(indv_hour['Date [YYYYMMDD]'])[4:] == i_days)
        std[i] = np.std(day_ind['Opacity'])
        mean[i] = np.mean(day_ind['Opacity'])
        plt.scatter(indv_hour['Hour UT'],(indv_hour['Opacity']), s = 15)
    plt.title(str(site)+' Hourly Measurements over '+str(day)+'for 2010 through 2017')
    plt.xlabel('Hour (UT)')
    plt.ylabel('Opacity')
    plt.show()
    return std, mean



def plot_stat(i_day, del_day):
    if len(str(i_day)) == 3:
        i_day = '0' + str(i_day)
    np.where(days ==
    for i in np.arange(scope_len):
        for i in np.arange
            error, mean = indexing_stat(scope_data[i][0], i_days)
        plt.errorbar(hours, mean, yerr = error, fmt = 'o')
        plt.title(str(scope_data[i][0])+' Hourly Average over January' +str(i_day) +'-' + str(i_day+ del_day) +'for 2008 through 2017')
        plt.xlabel('Hour (UT)')
        plt.ylabel('Opacity')
        #plt.savefig(str(scope_data[i][0])+'_Hourly_Average_over_June_2-4_for_2010_through_2017.png', dpi=300)
        plt.show()





for i in np.arange(scope_len):
    
    x = eval(str(scope_data[i][0])+'_hourly')['Date [YYYYMMDD]']
    y = eval(str(scope_data[i][0])+'_hourly')['Opacity']
    #ax.set_xticks([i for i in range(len(y))])
    plt.plot(y, label = str(scope_data[i][0]))



plt.ylabel('Opacity')
plt.title('April 2007-2017')
plt.legend()


years = [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017]

for i in np.arange(len(years)):
    np.where(int(str(eval(str(scope_data[i][0])+'_hourly')['Date [YYYYMMDD]'][i])[:4]) == years[i]
    x = eval(str(scope_data[i][0])+'_hourly')['Date [YYYYMMDD]']
    y = eval(str(scope_data[i][0])+'_hourly')['Opacity']

'''


######Figures of Merit
'''def alpha(date, hour, scopes):
    SITES = np.chararray(scope_len, 6)
    std = np.zeros(hour_len)
    mean = np.zeros(hour_len)
    
    for i in np.arange(hour_len):
        index = np.where(eval(str(site)+'_hourly')['Hour UT'] == hours[i])
        indv_hour = eval(str(site)+'_hourly')[index]
        std[i] = np.std(indv_hour['Opacity'])
        mean[i] = np.mean(indv_hour['Opacity'])
    
    for i in np.arange(scope_len):
        SITES[i] = scope_data[i][0]
        index = np.where(eval(str(site)+'_hourly')['Hour UT'] == hours[i])
        error, mean = indexing(scope_data[i][0])

    alpha  = np.sqrt(np.prod())

alpha(20170620, 00, scope_data[0][0:3])
'''
'''


means = (np.mean(ALMA0['Opacity']), np.mean(ALMA1['Opacity']), np.mean(ALMA2['Opacity']), np.mean(ALMA3['Opacity']), np.mean(ALMA4['Opacity']), np.mean(ALMA5['Opacity']), np.mean(ALMA6['Opacity']), np.mean(ALMA7['Opacity']), np.mean(ALMA8['Opacity']), np.mean(ALMA9['Opacity']), np.mean(ALMA10['Opacity']), np.mean(ALMA11['Opacity']), np.mean(ALMA12['Opacity']), np.mean(ALMA13['Opacity']), np.mean(ALMA14['Opacity']), np.mean(ALMA15['Opacity']), np.mean(ALMA16['Opacity']), np.mean(ALMA17['Opacity']), np.mean(ALMA18['Opacity']), np.mean(ALMA19['Opacity']), np.mean(ALMA20['Opacity']), np.mean(ALMA21['Opacity']), np.mean(ALMA22['Opacity']), np.mean(ALMA23['Opacity']))
errors = (std0, std1, std2,std3,std4,std5,std6,std7,std8,std9,std10,std11,std12,std13,std14,std15,std16,std17,std18,std19,std20,std21,std22, std23)
plt.errorbar(day, means, yerr = errors, fmt = 'o')'''



'''

SMT_hourly = tab.Table.read('SMT_hourly_20070401-20170431.dat', format='ascii')
APEX_hourly = tab.Table.read('APEX_hourly_20070401-20170431.dat', format='ascii')
ALMA_hourly = tab.Table.read('ALMA_hourly_20070401-20170431.dat', format='ascii')
Iram30_hourly = tab.Table.read('Iram30_hourly_20070401-20170431.dat', format='ascii')
JCMT_hourly = tab.Table.read('JCMT_hourly_20070401-20170431.dat', format='ascii')
LMT_hourly = tab.Table.read('LMT_hourly_20070401-20170431.dat', format='ascii')
SMA_hourly = tab.Table.read('SMA_hourly_20070401-20170431.dat', format='ascii')
SPT_hourly = tab.Table.read('SPT_hourly_20070401-20170431.dat', format='ascii')





'''
