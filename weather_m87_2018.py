# This script reads-in historical data files of opacity (tau) at
# zenith at various observational sites and plots statistical either
# percentile of tau or probability of tau<given_value as a function
# of time of a day for a given date range.
# Historical data files can be either from GFS or real EHT sites. 
# Hotaka Shiokawa, 12/29/2017

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import percentileofscore
import os
import math
from jd_calc import jd_calc
plt.rc('font',family='serif')
plt.rcParams['interactive'] = True
mpl.rcParams['figure.figsize'] = [10., 10.]


elev_cutoff = 15. # deg

#----- Function to read-in GFS data -----#
# Data files and data processing scripts are provided by Rodrigo Cordova.
def get_data_gfs( fname, site ):
    global time, tau, elev, elev_t, site_g

    # Use Pandas to read in csv data. Fast.
    d = pd.read_csv( fname, sep='\s+', header=None, names=['date','hr','site','freq','tau','tra'], skiprows=[0])
    d['yr'] = d['date'].astype(str).str[:4].astype(int)
    d['mon'] = d['date'].astype(str).str[4:6].astype(int)
    d['day'] = d['date'].astype(str).str[6:8].astype(int)
    d['minu'] = d['date']*0.
    d['sec'] = d['date']*0.

    # Construct "time" dataframe, time=['yr','mon','day','utc','mjd'] 
    mjd = jd_calc( d['yr'], d['mon'], d['day'] )
    utc = d['hr'] + d['minu']/60. + d['sec']/60./60.
    time = d.loc[:,'yr':'day']
    time['utc'] = utc
    time['mjd'] = mjd+utc/24.
    time['tau'] = d['tau']
    # Copy tau values to 'tau' for convenience
    tau = d['tau']

    # Memorize the site you just read-in
    site_g = site
    # Read-in m87* elevation time series at the site
    elev = np.load('elev_RC/%s_2018_m87_ele.npy'%(site)) # elevation
    elev_t = np.loadtxt('elev/elev_time.txt') # corresponding time in UTC



#----- Function to read-in LMT real data -----#
# CAUTION: The time zone in the read-in file is likely to be UTC, but could be local.
# The file is already processed by Dimitrios, but he doesn't quite remember the time zone either.
# We assume it is UTC for now, but need to be confirmed.
def get_data_LMT():
    global time, tau, elev, elev_t, site_g

    print('Reading in LMT data...')
    fname = '../real/LMT/lmtdata.dat'
    # Use Pandas to read in csv data. Fast.
    d = pd.read_csv( fname, sep='\s+', header=None,\
        names=['yr','mon','day','hr','minu','sec','i1','i2','i3','i4','i5','i6','i7','tau'])

    # Construct "time" dataframe, time=['yr','mon','day','utc','mjd'] 
    mjd = jd_calc( d['yr'], d['mon'], d['day'] )
    utc = d['hr'] + d['minu']/60. + d['sec']/60./60.
    time = d.loc[:,'yr':'day']
    time['utc'] = utc
    time['mjd'] = mjd+utc/24.

    # Add 'tau' column to the dataframe to filter out rows based on tau values
    time['tau'] = d['tau']
    time = time.loc[(time['tau']>0.) & (time['tau']<100.)] # filter out rows with strange tau values
    # Separately store 'tau' for convenience
    tau = np.array(time['tau'])

    # Memorize the site you just read-in
    site_g = 'LMT'
    # Read-in SgrA* elevation time series at the site
    elev = np.loadtxt('elev/elev_LMT_2018.txt') # elevation
    elev_t = np.loadtxt('elev/elev_time.txt') # corresponding time in UTC



#----- Function to read-in MaunaKea real data -----#
# CAUTION: MaunaKea's final result obtained by this script doesn't look correct;
# tau doesn't vary much over a day while it's supposed to.
# There could be a bug in this function or misunderstanding of the read-in data file format/unit..
# This definitely need to be checked, probably starting from plotting some example days' tau vs. utc.
# If the problem exists in "get_binned_data()", that's fatal for many other functions.
# Data is provided by Simon Radford. (The 'mk' directory in the link below).
# https://drive.google.com/drive/folders/1W0ugANyElV7gskl6jf1a4_sanTg_JFzX
def get_data_Maunakea():
    global time, tau, elev, elev_t, site_g

    print('Reading in Mauna Kea data...')
    cnt = 0
    for y in range(1950,2050+1,1): # year by year
        for m in range(1,12+1,1): # month by month
            fname = '../real/Maunakea/Opacity.%s%02d.CSO'%(str(y)[2:],m)
            if os.path.isfile(fname):
                # Use Pandas to read in csv data. Fast.
                d = pd.read_csv( fname, sep='\s+', header=None, skiprows=5, names=['yr','mon','day','ut','i1','i2','i3','tau','unc'])
                if cnt==0:
                    yr = d['yr']; mon = d['mon']; day = d['day']; utc = d['ut']; tau = d['tau']
                    mjd = jd_calc( d['yr'], d['mon'], d['day'] ) + d['ut']/24.
                else:
                    # Append each file's data to grand arrays
                    yr  = np.append( yr,d[ 'yr'])
                    mon = np.append(mon,d['mon'])
                    day = np.append(day,d['day'])
                    utc = np.append(utc,d[ 'ut'])
                    tau = np.append(tau,d['tau'])
                    mjd = np.append(mjd,jd_calc(d['yr'],d['mon'],d['day']) + d['ut']/24.)
                cnt += 1
            else: continue

    # Construct "time" dataframe.
    d = {'yr':yr, 'mon':mon, 'day':day, 'utc':utc, 'mjd':mjd, 'tau':tau}
    time = pd.DataFrame(data=d)
    # Use 'tau' column to filter out rows with strange tau values
    time = time.loc[(time['tau']>0.) & (time['tau']<100.)]
    # Separately store tau for convenience
    tau = np.array(time['tau'])

    # Memorize the site you just read-in
    site_g = 'MaunaKea'
    # Read-in SgrA* elevation time series at the site
    elev = np.loadtxt('elev/elev_SMA_2018.txt') # elevation
    elev_t = np.loadtxt('elev/elev_time.txt') # corresponding time in UTC



#----- Function to read-in ALMA 1995-2004 real data -----#
# Downloaded from 225GHz Transparency data column (text) at
# http://legacy.nrao.edu/alma/site/Chajnantor/data.c.html
def get_data_ALMA_1():
    global time, tau, elev, elev_t, site_g
    cnt = 0
    print('Reading in ALMA 1995-2004 data...')
    for y in range(1950,2050+1,1): # year by year
        for m in range(1,12+1,1): # month by month
            fname = '../real/ALMA/Opacity.%s%02d.C'%(str(y)[2:],m)
            if os.path.isfile(fname):
                # Use Pandas to read-in files. Fast
                d = pd.read_csv( fname, sep='\s+', header=None, skiprows=5, names=['yr','mon','day','ut','i1','i2','i3','tau'])
                if y<2000: yr_tmp = d['yr']+1900 # for the annoying 2-digits year....
                else: yr_tmp = d['yr']+2000
                if cnt==0:
                    yr = yr_tmp; mon = d['mon']; day = d['day']; utc = d['ut']; tau = d['tau']
                    mjd = jd_calc( d['yr'], d['mon'], d['day'] ) + d['ut']/24.
                else:
                    # Append each file's data to grand arrays
                    yr  = np.append(yr,yr_tmp)
                    mon = np.append(mon,d['mon'])
                    day = np.append(day,d['day'])
                    utc = np.append(utc,d[ 'ut'])
                    tau = np.append(tau,d['tau'])
                    mjd = np.append(mjd,jd_calc(d['yr'],d['mon'],d['day']) + d['ut']/24.)
                cnt += 1
            else: continue

    # Construct "time" dataframe
    d = {'yr':yr, 'mon':mon, 'day':day, 'utc':utc, 'mjd':mjd, 'tau':tau}
    time = pd.DataFrame(data=d)
    # Use 'tau' column to filter out rows with strange tau values
    time = time.loc[(time['tau']>0.) & (time['tau']<100.)] # filter out rows with strange tau values
    # Separately store tau for convenience
    tau = np.array(time['tau'])

    # Memorize the site you just read-in
    site_g = 'ALMA_1'
    # Read-in SgrA* elevation time series at the site
    elev = np.loadtxt('elev/elev_ALMA_2018.txt') # elevation
    elev_t = np.loadtxt('elev/elev_time.txt') # corresponding time in UTC

    # Return "time" dataframe so that it can be combined with ALMA_2 if necessary
    return time



#----- Function to read-in ALMA 2007-2017 real data -----#
# Downloaded from
# http://archive.eso.org/wdb/wdb/eso/meteo_apex/form
# PWV - > Tau relation at the site is calculated by Scott Paine:
# tau = 0.0375*d['pwv'] + 0.0115
def get_data_ALMA_2():
    global time, tau, elev, elev_t, site_g
    fname = '../real/ALMA/2007-2017.csv'
    print('Reading in ALMA 2007-2017 data...')
    # Use Pandas. Fast
    d = pd.read_csv( fname, header=None, skiprows=2, names=['time','pwv'], comment='#')

    d['yr'] = d['time'].str[:4].astype(int)
    d['mon'] = d['time'].str[5:7].astype(int)
    d['day'] = d['time'].str[8:10].astype(int)
    d['hr'] = d['time'].str[11:13].astype(int)
    d['minu'] = d['time'].str[14:16].astype(int)
    d['sec'] = d['time'].str[17:19].astype(int)
    # Convert PWV to Tau. The relation is calculated by Scott Paine
    tau = 0.0375*d['pwv'] + 0.0115

    # Construct "time" dataframe, time=['yr','mon','day','utc','mjd','tau'] 
    utc = d['hr'] + d['minu']/60. + d['sec']/60./60.
    mjd = jd_calc( d['yr'], d['mon'], d['day'] ) + utc/24.
    time_tmp = {'yr':d['yr'], 'mon':d['mon'], 'day':d['day'], 'utc':utc, 'mjd':mjd, 'tau':tau}
    time = pd.DataFrame(data=time_tmp)
    # Separately store tau for convenience
    tau = np.array(time['tau'])

    # Memorize the site you just read-in
    site_g = 'ALMA_2'
    # Read-in SgrA* elevation time series at the site
    elev = np.loadtxt('elev/elev_ALMA_2018.txt') # elevation
    elev_t = np.loadtxt('elev/elev_time.txt') # corresponding time in UTC

    # Return "time" dataframe so that it can be combined with ALMA_1 if necessary
    return time



#----- Function to read-in ALMA 1995-2017 real data -----#
# Run get_data_ALMA_1() and get_data_ALMA_2(), and combine.
def get_data_ALMA():
    global time, tau, elev, elev_t, site_g

    # Construct "time" dataframe, time=['yr','mon','day','utc','mjd','tau'] 
    time = get_data_ALMA_1() # 1995-2004
    time = time.append(get_data_ALMA_2(),ignore_index=True) # 2007-2017

    # Separately store tau for convenience
    tau = time['tau']
    # Memorize the site you just read-in
    site_g = 'ALMA'
    # Read-in SgrA* elevation time series at the site
    elev = np.loadtxt('elev/elev_ALMA_2018.txt') # elevation
    elev_t = np.loadtxt('elev/elev_time.txt') # corresponding time in UTC



def get_data_general(site,GFS=False):
    if GFS:
        if site=='LMT': get_data_gfs('../GFS/10yr_data/LMT10YR_8-8-17.dat','LMT')
        elif site=='ALMA': get_data_gfs('../GFS/10yr_data/ALMA10YR_8-8-17.dat','ALMA')
        elif site=='SMA': get_data_gfs('../GFS/10yr_data/SMA10YR_8-8-17.dat','SMA')
        elif site=='Maunakea': get_data_gfs('../GFS/10yr_data/SMA10YR_8-8-17.dat','SMA')
        elif site=='SMT': get_data_gfs('../GFS/10yr_data/SMT10YR_8-8-17.dat','SMT')
        elif site=='SPT': get_data_gfs('../GFS/10yr_data/SPT10YR_8-8-17.dat','SPT')
        elif site=='PV': get_data_gfs('../GFS/10yr_data/Iram3010YR_8-8-17.dat','PV')
        #elif site=='PDB': get_data_gfs('../GFS/data_april/PBD_hourly_20070401-20170431.dat','PDB')
        else:
            print('Site data not available: %s'%site)
            return False
        return True
    else:
        if site=='LMT': get_data_LMT()
        elif site=='ALMA': get_data_ALMA()
        elif site=='SMA': get_data_Maunakea()
        elif site=='Maunakea': get_data_Maunakea()
        else:
            print('Site data not available: %s'%site)
            return False
        return True



#----- Function to return time binned tau values for a specified date range -----#
# Input parameters:
# m1,d1 = month and date of starting date.
# m2,d2 = month and date of ending date.
# dt = time bin size in hr
# get_data_*() must be called in advance, i.e.. and 'time', 'tau', 'elev', 'elev_t' need to be all ready. 
def get_binned_data(m1,d1,m2,d2,dt):

    # Get a list of row index in 'time' for the specified dates
    if m1==m2:
        idx = time.query('mon==%d and day>=%d and day<=%d'%(m1,d1,d2)).index
    elif m2==(m1+1):
        idx = time.query('(mon==%d and day>=%d) or (mon==%d and day<=%d)'%(m1,d1,m2,d2)).index
    elif m2>(m1+1):
        print('Time range is set too wide')
        return
    elif m2<m1:
        print('m2 > m1 required')
        return
    if len(idx)==0:
        print('No matching dates')
        return

    # For convenience (cumbersome to write 'time['*']' everytime)
    yr  = time['yr' ]
    mon = time['mon']
    day = time['day']
    utc = time['utc']


    # The routine below calculates the following arrays.
    # t=[[list of t of day 1],[list of t of day 2],...]
    # ta=[[list of tau of day 1],[list of tau of day 2],...]
    # moda=[[yr of day1, mon of day1, day1],[yr of day2, mon of day2, day2],....]

    t = []
    ta = []
    i = idx[0]
    t_tmp = [utc[i]]
    ta_tmp = [tau[i]]
    d = day[i]
    moda = [[yr[i],mon[i],d]]
    for i in idx[1:]:
        if d!=day[i]: # if date differs from previous step...
            t.append(t_tmp)
            ta.append(ta_tmp)
            moda.append([yr[i],mon[i],day[i]])
            t_tmp = []
            ta_tmp = []
            d = day[i]
        t_tmp.append(utc[i])
        ta_tmp.append(tau[i])
    # Don't forget to append the final *_tmp lists
    t.append(t_tmp)
    ta.append(ta_tmp)


    # The routine below calculates the following arrays.
    # t2=[[list of t of day 1 every dt],[list of t of day 2 every dt],...] = [np.arange(0,24,dt),np.arange(0,24,dt),...]
    # ta2_orig=[[time-binned list of tau of day 1],[time-binned list of tau of day 2],...]
    # ta2=[[time-binned list of tau/sin(elev) of day 1],[time-binned list of tau/sin(elev) of day 2],...]
    # Here, time-binning is done by averaging tau values that fall in each time bin.
    # For example, say
    # t  = [[0,2,6,14,15],[10,11,15,20],...]
    # ta = [[3,4,5, 4, 3],[ 5, 6, 7, 6],...] (note these numbers are too high for realistic tau)
    # then, for dt=8
    # t2  = [[0,   8,   16],[   0, 8, 16],...]
    # ta2 = [[4, 3.5, None],[None, 6,  6],...]

    t2 = [] # 2d time array ('total # of days' x '# of time bins in a day')
    ta2 = [] # 2d array of tau/sin(elev) corresponding to t2
    ta2_orig = [] # 2d array of tau corresponding to t2 (orig=original, meaning just the original tau value)
    elev_list = [] # 2d array of elevation of SgrA* corresponding to t2
    nbin = int(np.round(24./dt)) # number of binns per day

    for i in range(len(t)): # day by day loop
        td = np.array(t[i]) # list of time of day-i
        taud = np.array(ta[i]) # list of tau of day-i corresponding to td
        tau_day = [None]*nbin # prepare empty list for time-binned tau/sin(elev) for day-i
        tau_day_orig = [None]*nbin # prepare empty list for time-binned tau for day-i
        # Calculate MJD of day-i to find index to refer to in the 'elev' array
        mjd_2018 = jd_calc( 2018, moda[i][1], moda[i][2] )
        i_elev = int(mjd_2018 - elev[0][0]) # index of elev that corresponds to mjd_2018
        elev_day = [] # for storing interpolated SgrA* elevation to time-bin's center for day-i
        for u_s in np.arange(0.,24.,dt): # for each time bin
            ibin = int(u_s/dt) # bin index
            u_e = u_s + dt # end point of this bin
            # Interpolate elevation value to this bin's center and append to elev_day
            th = np.interp(u_s+dt/2.,elev_t,elev[i_elev][1:])
            elev_day.append(th)
            li = (td>u_s)*(td<=u_e) # index list of data points that falls in this bin
            if u_s<=1.e-10: li = (td>=u_s)*(td<=u_e)
            if li.sum()>0: # 'if data points exist in this time bin'. Otherwise left as 'None'
                # Average all tau/sin(elev) values in this bin
                tau_day[ibin] = np.average(taud[li]/np.sin(th*np.pi/180.))
                # Average all tau values in this bin
                tau_day_orig[ibin] = np.average(taud[li])
        t2.append(np.arange(0.,24.,dt))
        # Append the binned data of day-i to the grand arrays
        ta2.append(tau_day)
        ta2_orig.append(tau_day_orig)
        # Append the elevation of SgrA* interpolated to bin's center of day-i to the grand elev array
        elev_list.append(elev_day)


    # The routine below calculates the following arrays.
    # t3 = np.arange(0,24,dt)
    # ta3 = np.transpose(ta2) but None removed
    #  = [[list of tau of all dates in 1st time bin],[... for 2nd time bin],[],[]...[... for 24/dt time bin]]

    # For tau/sin(elev)...
    ta2 = np.transpose(ta2)
    ta3 = [[ ta2[i][j] for j in range(len(ta2[i])) if ta2[i][j] is not None] for i in range(len(ta2))] # remove None
    # For tau
    ta2_orig = np.transpose(ta2_orig)
    ta3_orig = [[ ta2_orig[i][j] for j in range(len(ta2_orig[i])) if ta2_orig[i][j] is not None] for i in range(len(ta2_orig))] # remove None

    t3 = t2[0] # = np.arange(0,24,dt)

    # For each time bin, find the minimum elevation of SgaA* among all the dates in consideration
    elev_list = np.transpose(elev_list)
    elev_min = [min(el) for el in elev_list] # list of minimum possible SgrA* elevation for every time bin

    return t3,ta3,ta3_orig,moda,elev_min



#----- Function to customize plot tick parameters and frame border width-----#
def tick_customize(f, ma_size, ma_width, ma_labelsize, mi_size, mi_width, mi_labelsize, border_width):

    f.tick_params(axis='both', size=ma_size, width=ma_width, which='major', labelsize=ma_labelsize)
    f.tick_params(axis='both', size=mi_size, width=mi_width, which='minor', labelsize=mi_labelsize)

    plt.gca().spines['bottom'].set_linewidth( border_width )
    plt.gca().spines['top'   ].set_linewidth( border_width )
    plt.gca().spines['left'  ].set_linewidth( border_width )
    plt.gca().spines['right' ].set_linewidth( border_width )



#----- Function to plot t vs. tau for given percentile values and a date range -----#
# Input parameters:
# m1,d1 = month and date of starting date.
# m2,d2 = month and date of ending date.
# dt = time bin size in hr
# site = site of your choice, 'ALMA','LMT','Maunakea'(or 'SMA') for GFS=False,
#  additional choices of 'SMT','SPT','PV' for GFS=True
# Choose list of percentile values to be plotted in per_list below.
# Choose elevation cutoff below: tau/sin(elev) with elev<elev_cutoff won't be plotted
per_list = np.arange(25,100,25) # percent

def plot_percentile(m1,d1,m2,d2,dt,site,GFS=False):

    get_data_general(site,GFS)

    # Get time binned data tau (ta3_orig) and tau/sin(elev) (ta3), where time bins are np.arange(0,24,dt)
    # See the comments in get_binned_data() for details of the variables.
    t3,ta3,ta3_orig,moda,elev_min = get_binned_data(m1,d1,m2,d2,dt)


    # The routine below calculates the following arrays.
    # ta4=[[time series of tau/sin(elev) for percentile value 1 in per_list], [...for percentile value 2...],...,[]]
    # t4=[[list of t corresponding to contents of ta4[0]],[...of ta4[1]],...,[]]
    # ta4_orig=[[time series of tau for percentile value 1], [...for percentile value 2...],...,[]]
    # t4_orig=[[list of t corresponding to contents of ta4_orig[0]],[...of ta4_orig[1]],...,[]]
    # Note, any time-bin in ta3 that contains even a single tau/sin(elev) with elev<elev_cutoff are excluded to be added to ta4

    t4_tmp_tmp = t3+dt/2.
    # Prepare the empty grand arrays... their len() will be equal to len(per_list)
    t4 = []
    t4_orig = []
    ta4 = []
    ta4_orig = []
    # Go over each percentile value in per_list
    for per in per_list:
        ta4_tmp = []
        t4_tmp = []
        for i in range(len(ta3)):
            if len(ta3[i])>0 and elev_min[i]>elev_cutoff:
                ta4_tmp.append( np.percentile(ta3[i],per) )
                t4_tmp.append(t4_tmp_tmp[i])
        ta4.append(np.array(ta4_tmp))
        t4.append(np.array(t4_tmp))

        ta4_orig_tmp = []
        t4_orig_tmp = []
        for i in range(len(ta3_orig)):
            if len(ta3_orig[i])>0:
                ta4_orig_tmp.append( np.percentile(ta3_orig[i],per) )
                t4_orig_tmp.append(t4_tmp_tmp[i])
        ta4_orig.append(np.array(ta4_orig_tmp))
        t4_orig.append(np.array(t4_orig_tmp))

    # find time range
    modat = np.transpose(np.array(moda))
    modat_reduce = modat[1]*31+modat[2]
    date_min = moda[np.argmin(modat_reduce)] # actual earliest date (among any year)
    date_max = moda[np.argmax(modat_reduce)] # actual latest date (among any year)
    yr_min = min(modat[0]) # earliest year among all the dates
    yr_max = max(modat[0]) # latest year among all the dates

    # plotting
    xmin = 0
    xmax = 24
    ymax = 1.

    fig = plt.figure()
    f1 = fig.add_subplot(211)
    for i, per in enumerate(per_list):
        f1.plot(t4[i],ta4[i],'o-',lw=2,label='%d%%'%per)
    f1.grid()
    f1.set_xlim(xmin,xmax)
    f1.set_ylim(0,ymax)
    f1.set_ylabel(r'$\tau$ / sin( elevation )',fontsize=15)
    f1.set_xlabel('UTC',fontsize=15)
    #f1.set_xlabel('Time [hr]',fontsize=15)
    f1.set_title(r'%d/%d - %d/%d (%d-%d), %s'\
        %(date_min[1],date_min[2],date_max[1],date_max[2],yr_min,yr_max,site_g),fontsize=20)
    plt.legend(loc='upper left',frameon=False)
    tick_customize(f1, 8.,2.,15.,5.,2.,15.,2.)

    f2 = fig.add_subplot(212)
    for i, per in enumerate(per_list):
        f2.plot(t4_orig[i],ta4_orig[i],'o-',lw=2,label='%d%%'%per)
    f2.grid()
    f2.set_xlim(xmin,xmax)
    f2.set_ylim(0,ymax)
    f2.set_ylabel(r'$\tau$',fontsize=30)
    f2.set_xlabel('UTC',fontsize=15)
    #f2.set_xlabel('Time [hr]',fontsize=15)
    plt.legend(loc='upper left',frameon=False)
    tick_customize(f2, 8.,2.,15.,5.,2.,15.,2.)

    f12 = f1.twinx()
    mjd_2018_s = int(jd_calc( 2018, date_min[1], date_min[2] ))
    mjd_2018_e = int(jd_calc( 2018, date_max[1], date_max[2] ))
    for mjd_2018 in range(mjd_2018_s,mjd_2018_e+1):
        i_elev = int(mjd_2018 - elev[0][0])
        f12.plot(elev_t,elev[i_elev][1:],'m-',alpha=0.2)
    f12.set_xlim(xmin,xmax)
    f12.set_ylim(0,90)
    f12.set_ylabel('Elevation [deg]',fontsize=15)
    f12.plot([30,30],[0,1],'m-',label='M87* Elevation') # just for legend
    plt.legend(loc='upper right',frameon=False)
    tick_customize(f12, 8.,2.,15.,5.,2.,15.,2.)
    plt.savefig(r'M87*_Plot_Percentile_for_'+str(site)+'_during_'+str(m1)+'-'+str(d1)+'_to_'+str(m2)+'-'+str(d2))


def plot_percentile_all(m1,d1,m2,d2,dt,GFS=False):

    plot_percentile(m1,d1,m2,d2,dt,'LMT',GFS)
    plot_percentile(m1,d1,m2,d2,dt,'ALMA',GFS)
    plot_percentile(m1,d1,m2,d2,dt,'Maunakea',GFS)

    if GFS==True:
        plot_percentile(m1,d1,m2,d2,dt,'SMT',GFS)
        plot_percentile(m1,d1,m2,d2,dt,'SPT',GFS)
        plot_percentile(m1,d1,m2,d2,dt,'PV',GFS)



def plot_probability(m1,d1,m2,d2,dt,val,GFS=False):

    if GFS:
        sites = ['LMT','ALMA','SMA','SMT','PV','SPT']
    else:
        sites = ['LMT','ALMA','SMA']

    da = []; da_orig = []; moda = []; elev = []; elev_min = []
    for site in sites:
        avail = get_data_general(site,GFS=GFS)
        if avail:
            tbins,da_temp,da_orig_temp,moda_temp,elev_min_temp = get_binned_data(m1,d1,m2,d2,dt)
            elev_temp = np.load('elev_RC/%s_2018_m87_ele.npy'%site)
            da.append(da_temp)
            da_orig.append(da_orig_temp)
            moda.append(moda_temp)
            elev.append(elev_temp)
            elev_min.append(elev_min_temp)
        else:
            return

    # Plots probability that elevation taken into account
    prob = [[] for i in range(len(sites))]
    for i_t in range(len(tbins)):
        for i_d in range(len(da)):
            prob_temp = -1.
            if len(da[i_d][i_t])>0 and elev_min[i_d][i_t]>elev_cutoff:
                prob_temp = percentileofscore(da[i_d][i_t],val)/100.
            prob[i_d].append(prob_temp)

    c = ['r','b','g','m','y','c','k']
    fig = plt.figure()
    f1 = fig.add_subplot(211)
    tbins = np.array(tbins)
    for i, site in enumerate(sites):
        li = np.array(prob[i])>0.
        f1.plot(tbins[li],np.array(prob[i])[li],'%so-'%c[i],lw=2.,label=site)

    prob = np.transpose(prob)
    prob_tot = []
    t_tot = []
    for i, p in enumerate(prob):
        if any(np.array(p)>0.):
            prob_tot.append( np.prod(np.array(p)[np.array(p)>0.]) )
            t_tot.append(tbins[i])
    f1.plot(t_tot,prob_tot,'k--',lw=1.,label='Net Prob.')

    xmax = 27.
    ymax = 1.05
    f1.set_xlim(0.,xmax)
    f1.set_ylim(0.,ymax)
    f1.set_xlabel('UTC [hr]',fontsize=17.)
    f1.set_ylabel(r'Probability: $\tau$ / sin( elevation ) < %1.1f'%val,fontsize=15)
    tick_customize(f1, 8.,2.,15.,5.,2.,15.,2.)
    #f1.plot([0,24],[-1000,-1000],'k--',lw=1.,label='Net Prob.')
    f1.plot([0,xmax],np.array([elev_cutoff,elev_cutoff])/90.*ymax,'k-.',lw=1.,label='Elev. cutoff')
    plt.legend(loc='upper right',frameon=False)

    # Plots pure probability
    prob_orig = [[] for i in range(len(sites))]
    for i_t in range(len(tbins)):
        for i_d in range(len(da_orig)):
            prob_temp = -1.
            if len(da_orig[i_d][i_t])>0:
                prob_temp = percentileofscore(da_orig[i_d][i_t],val)/100.
            prob_orig[i_d].append(prob_temp)

    f2 = fig.add_subplot(212)
    for i, site in enumerate(sites):
        li = np.array(prob_orig[i])>=0.
        f2.plot(tbins[li],np.array(prob_orig[i])[li],'%so-'%c[i],lw=2.,label=site)

    prob_orig = np.transpose(prob_orig)
    prob_orig_tot = []
    t_tot = []
    for i, p in enumerate(prob_orig):
        if any(np.array(p)>0.):
            prob_orig_tot.append( np.prod(np.array(p)[np.array(p)>0.]) )
            t_tot.append(tbins[i])
    f2.plot(t_tot,prob_orig_tot,'k--',lw=1.,label='Net Prob.')

    f2.set_xlim(0.,xmax)
    f2.set_ylim(0.,ymax)
    f2.set_xlabel('UTC [hr]',fontsize=17.)
    f2.set_ylabel(r'Probability: $\tau$ < %1.1f'%val,fontsize=17)
    tick_customize(f2, 8.,2.,15.,5.,2.,15.,2.)
    plt.legend(loc='lower right',frameon=False)


    # find time range
    elev_t = np.loadtxt('elev/elev_time.txt')
    f12 = f1.twinx()
    for i, m in enumerate(moda):
        modat = np.transpose(np.array(m))
        modat_reduce = modat[1]*31+modat[2]
        date_min = m[np.argmin(modat_reduce)] # actual earliest date (among any year)
        date_max = m[np.argmax(modat_reduce)] # actual latest date (among any year)

        mjd_2018_s = int(jd_calc( 2018, date_min[1], date_min[2] ))
        mjd_2018_e = int(jd_calc( 2018, date_max[1], date_max[2] ))
        for mjd_2018 in range(mjd_2018_s,mjd_2018_e+1):
            i_elev = int(mjd_2018 - elev[i][0][0])
            if mjd_2018==mjd_2018_s:
                li = elev[i][i_elev][1:]>elev_cutoff
            else:
                li *= elev[i][i_elev][1:]>elev_cutoff
        for mjd_2018 in range(mjd_2018_s,mjd_2018_e+1):
            i_elev = int(mjd_2018 - elev[i][0][0])
            f12.plot(elev_t[li],elev[i][i_elev][1:][li],color=c[i],alpha=0.2)
        #f12.plot([0,24],[-1000,-1000],color=c[i],label=sites[i])

    #f12.plot([0,24],[-1000,-1000],'k--',lw=1.,label='Net Prob.')
    #f12.plot([0,24],[elev_cutoff,elev_cutoff],'k-.',lw=1.,label='Elev. cutoff')
    f12.set_xlim(0.,xmax)
    f12.set_ylim(0,90)
    f12.set_ylabel('Elevation [deg]',fontsize=15)
    #plt.legend(loc='upper right',frameon=False)
    tick_customize(f12, 8.,2.,15.,5.,2.,15.,2.)

    f1.set_title(r'%d/%d - %d/%d'%(m1,d1,m2,d2),fontsize=20)
    plt.savefig(r'M87_Plot_Probability_for_%d-%d_to_%d-%d'%(m1,d1,m2,d2))

"""
def plot_elev(m1,d1,m2,d2,site):

    elev = np.loadtxt('elev/elev_%s_2018.txt'%site)
    elev_t = np.loadtxt('elev/elev_time.txt')

    fig = plt.figure()
    f = fig.add_subplot(111)
    mjd_2018_s = int(jd_calc( 2018, m1, d1 ))
    mjd_2018_e = int(jd_calc( 2018, m2, d2 ))
    for mjd_2018 in range(mjd_2018_s,mjd_2018_e+1):
        i_elev = int(mjd_2018 - elev[0][3])
        f.plot(elev_t,elev[i_elev][4:],'m-',alpha=0.2)
    f.set_xlim(0,24)
    f.set_ylim(0,90)
    f.set_ylabel('Elevation [deg]',fontsize=15)
    #plt.legend(loc='upper right',frameon=False)
    tick_customize(f, 8.,2.,15.,5.,2.,15.,2.)


# Inspired by answer 2 at:
# https://stackoverflow.com/questions/10997577/python-timezone-conversion
# List of time zones can be found by pytz.all_timezones.
# The problem of this function is it's not vectorized and takes infinite
# amount of time to process real data...
from datetime import datetime
import pytz
def local_to_utc(tz, dt_in):
    tz1 = pytz.timezone(tz)
    tz2 = pytz.timezone('UTC')

    dt_in = np.transpose(dt_in)
    for i, dt in enumerate(dt_in):
        t = tz1.localize(datetime(dt[0],dt[1],dt[2],dt[3],dt[4],dt[5]))
        t = t.astimezone(tz2)
        dt_tmp = np.array([int(t.strftime('%Y')),int(t.strftime('%m')),int(t.strftime('%d')),\
            int(t.strftime('%H')),int(t.strftime('%M')),int(t.strftime('%S'))])
        if i==0: dt_out=dt_tmp
        else: dt_out = np.vstack((dt_out,dt_tmp))
    dt_out = np.transpose(dt_out)
    return dt_out
"""

