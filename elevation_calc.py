import os
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import datetime as dt
from scipy.interpolate import interp1d
import argparse
from datetime import date
from dateutil.rrule import rrule, DAILY
from jd_calc import jd_calc



scope_data = [['SMT',32.7016, -109.891, 3185],['LMT', 18.9858, -97.3147, 4640],['JCMT', 19.8228, -155.477, 4092],['SMA', 19.8242, -155.478, 4080],['Iram30', 37.0661, -3.3925, 2850], ['SPT', -90, 0, 2800],['ALMA', -23.0193, -67.7532, 5058.7],['APEX', -23.0058, -67.7592, 5100]]
scopelen = len(scope_data)
'''

sgra_ra = 17.0 + ((45.0 + (40.0409/60.0))/60.0)
m87_ra = 12.0 + ((30.0 + (49.42338/60.0))/60.0)
sgra_dec = -(29.0 + ((0.0 + (28.118/60.0))/60.0)) #deg
m87_dec = 12.0 + ((23.0 + (28.0439/60.0))/60.0) #deg

s_RA_deg = sgra_ra *15
m_RA_deg = m87_ra  *15
'''
def j2000(year, month, day, hour, minute):
    days = []
    if year < 2000:
        a = date(year, month, day)
        b = date(2000, 1, 1)
    if year >= 2000:
        b = date(year, month, day)
        a = date(2000, 1, 1)
    

    for dt in rrule(DAILY, dtstart=a, until=b):
        days.append(int(dt.strftime("%Y%m%d")))
    if year < 2000:
        j2000_days = -(len(days)+1.5)+2 #there is still some bug here from the date indexing for negative dates from j2000
    if year >= 2000:
        j2000_days = len(days)-1.5
    j2000_minute = minute/60.0
    j2000_hour = (float(hour) + j2000_minute)/24.0
    d = j2000_hour + j2000_days
    return d

def LST(j2000_d , longi, UT):
    lst = 100.46 + (0.985647* j2000_d) + longi + (15.0*UT)
    while lst < 0.0:
        lst = lst + 360.0
    while lst >= 360.0:
        lst = lst - 360.0
    return lst

def elevation(y, m, d, h, minu, longi, lat, ra_deg, decli):
    d = j2000(y, m, d, h, minu)
    lst = LST(d, longi, float(h) + (minu/60.0))
    HA = lst - ra_deg
    longi = longi * np.pi/180.0
    lat = lat * np.pi/180.0
    decli = decli * np.pi/180.0
    while HA < 0.0:
        HA = HA + 360.0
    while HA > 360.0:
        HA = HA - 360.0
    HA = HA * np.pi/180.0
    alt = np.arcsin((np.sin(decli)*np.sin(lat))+ (np.cos(decli)*np.cos(lat)*np.cos(HA)))
    alt = alt * 180.0/np.pi
    return alt


def sgra_elevation(y,m,d,h,minu,scope_info):
    sgra_ra = 17.0 + ((45.0 + (40.0409/60.0))/60.0)
    sgra_dec = -(29.0 + ((0.0 + (28.118/60.0))/60.0)) #deg
    s_RA_deg = sgra_ra *15
    ele = elevation(y, m, d, h, minu, scope_info[2], scope_info[1], s_RA_deg, sgra_dec)
    return ele


def m87_elevation(y,m,d,h,minu,scope_info):
    m87_ra = 12.0 + ((30.0 + (49.42338/60.0))/60.0)
    m87_dec = 12.0 + ((23.0 + (28.0439/60.0))/60.0) #deg
    m_RA_deg = m87_ra  *15
    ele = elevation(y, m, d, h, minu, scope_info[2], scope_info[1], m_RA_deg, m87_dec)
    return ele


def day(y1, m1, d1, y2, m2, d2, yf):
    days = []
    i = 0
    while i <= (yf -y1):
        a = date(y1+i, m1, d1)
        b = date(y2+i, m2, d2)

        for dt in rrule(DAILY, dtstart=a, until=b):
            days.append(int(dt.strftime("%Y%m%d")))
        i = i+1
    return np.array(days)


#days = day(args.year, args.mon, args.day, args.year2, args.mon2, args.day2, args.fyear)


'''
SMT10 = Table.read('/Users/rodrigo/Documents/EHT/SMT10YR_8-8-17.dat',format = 'ascii')
LMT10 = Table.read('/Users/rodrigo/Documents/EHT/LMT10YR_8-8-17.dat',format = 'ascii')
JCMT10 = Table.read('/Users/rodrigo/Documents/EHT/JCMT10YR_8-8-17.dat',format = 'ascii')
SMA10 = Table.read('/Users/rodrigo/Documents/EHT/SMA10YR_8-8-17.dat',format = 'ascii')
Iram3010 = Table.read('/Users/rodrigo/Documents/EHT/Iram3010YR_8-8-17.dat',format = 'ascii')
SPT10 = Table.read('/Users/rodrigo/Documents/EHT/SPT10YR_8-8-17.dat',format = 'ascii')
ALMA10 = Table.read('/Users/rodrigo/Documents/EHT/ALMA10YR_8-8-17.dat',format = 'ascii')
APEX10 = Table.read('/Users/rodrigo/Documents/EHT/APEX10YR_8-8-17.dat',format = 'ascii')
'''

def sgr_el_1day(yr, mnth, dy, site):
    minutes = np.linspace(0,60, 11)[:-1]
    hours = np.arange(0,24)
    a = len(minutes)*len(hours)
    hs = np.zeros(a)
    ms = np.zeros(a)
    elevation = np.zeros(a)
    Date = np.zeros(a)
    b = 0
    for j in np.arange(len(hours)):
        for i in np.arange(len(minutes)):
             elevation[b]= sgra_elevation(yr, mnth, dy, hours[j], minutes[i], scope_data[site])
             hs[b] = hours[j]
             ms[b] = minutes[i]
             if (mnth >= 10) & (dy >= 10):
                Date[b] = int(str(yr)+str(mnth)+str(dy))
             if (mnth < 10) & (dy >= 10):
                Date[b] = int(str(yr)+str('0' + str(int(mnth)))+str(int(dy)))
             if (mnth >= 10) & (dy < 10):
                Date[b] = int(str(yr)+str(int(mnth))+str('0' + str(int(dy))))
             if (mnth < 10) & (dy < 10):
                Date[b] = int(str(yr)+str('0' + str(int(mnth)))+str('0' + str(int(dy))))
             Date[b] = int(Date[b])
             b = b +1
    #els = Table([Date, hs + (ms/60.0), elevation], names = ('Date', 'Hour', 'SgrA* Elevation'))
    #els = els[np.where(els['SgrA* Elevation'] > 15.0)]

    return elevation

def sgr_ele_year(year, site):
    days = day(year,1,1,year,12,31,year)
    year_ele = []
    for i in np.arange(len(days)):
        print days[i]
        a = []
        a.append(jd_calc(year, int(str(days[i])[4:6]), int(str(days[i])[6:8])))
        ele = sgr_el_1day(year, int(str(days[i])[4:6]), int(str(days[i])[6:8]), site)
        for i in np.arange(len(ele)):
            a.append(ele[i])
        year_ele.append(a)
    return year_ele

def m87_el_1day(yr, mnth, dy, site_no):
    minutes = np.linspace(0,60, 11)[:-1]
    hours = np.arange(0,24)
    a = len(minutes)*len(hours)
    hs = np.zeros(a)
    ms = np.zeros(a)
    elevation = np.zeros(a)
    Date = np.zeros(a)
    b = 0
    for j in np.arange(len(hours)):
        for i in np.arange(len(minutes)):
             elevation[b]= m87_elevation(yr, mnth, dy, hours[j], minutes[i], scope_data[site_no])
             b = b +1
    return elevation

def m87_ele_year(year, site):
    days = day(year,1,1,year,12,31,year)
    year_ele = []
    for i in np.arange(len(days)):
        print days[i]
        a = []
        a.append(jd_calc(year, int(str(days[i])[4:6]), int(str(days[i])[6:8])))
        ele = m87_el_1day(year, int(str(days[i])[4:6]), int(str(days[i])[6:8]), site)
        for i in np.arange(len(ele)):
            a.append(ele[i])
        year_ele.append(a)
    return year_ele

#M87 year elevations
'''
SMT_year_87 = m87_ele_year(2018,0)
np.save('SMT_2018_m87_ele.npy', SMT_year_87)

SMT_year_sgra19 = sgr_ele_year(2019,0)
np.save('SMT_2019_sgra_ele.npy', SMT_year_sgra19)
#########################################################
ALMA_year_87 = m87_ele_year(2018,6)
np.save('ALMA_2018_m87_ele.npy', ALMA_year_87)

ALMA_year_sgra = sgr_ele_year(2018,6)
np.save('ALMA_2018_sgra_ele.npy', ALMA_year_sgra)
ALMA_year_sgra19 = sgr_ele_year(2019,6)
np.save('ALMA_2019_sgra_ele.npy', ALMA_year_sgra19)
#########################################################
LMT_year_87 = m87_ele_year(2018,1)
np.save('LMT_2018_m87_ele.npy', LMT_year_87)

LMT_year_sgra = sgr_ele_year(2018,1)
np.save('LMT_2018_sgra_ele.npy', LMT_year_sgra)
LMT_year_sgra19 = sgr_ele_year(2019,1)
np.save('LMT_2019_sgra_ele.npy', LMT_year_sgra19)
#########################################################
Iram_year_87 = m87_ele_year(2018,4)
np.save('Iram_2018_m87_ele.npy', Iram_year_87)

Iram_year_sgra = sgr_ele_year(2018,4)
np.save('Iram_2018_sgra_ele.npy', Iram_year_sgra)
Iram_year_sgra19 = sgr_ele_year(2019,4)
np.save('Iram_2019_sgra_ele.npy', Iram_year_sgra19)
#########################################################
SPT_year_87 = m87_ele_year(2018,5)
np.save('SPT_2018_m87_ele.npy', SPT_year_87)

SPT_year_sgra = sgr_ele_year(2018,5)
np.save('SPT_2018_sgra_ele.npy', SPT_year_sgra)
SPT_year_sgra19 = sgr_ele_year(2019,5)
np.save('SPT_2019_sgra_ele.npy', SPT_year_sgra19)

#########################################################
SMA_year_87 = m87_ele_year(2018,3)
np.save('SMA_2018_m87_ele.npy', SMA_year_87)

SMA_year_sgra = sgr_ele_year(2018,3)
np.save('SMA_2018_sgra_ele.npy', SMA_year_sgra)
SMA_year_sgra19 = sgr_ele_year(2019,3)
np.save('SMA_2019_sgra_ele.npy', SMA_year_sgra19)


'''
'''

def ut_offset():
    lons = np.zeros(scopelen)
    diff = np.zeros(scopelen)
    for i in np.arange(scopelen):
        lons[i] = scope_data[i][2]
        diff[i] = lons[i] * (4.0/60.0)
    for i in np.arange(scopelen):
        lons[i]
sgralmt = Table.read('elevationSgr.dat', format = 'ascii')
m87lmt = Table.read('elevationM87.dat', format = 'ascii')
'''
