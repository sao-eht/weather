#!/usr/bin/env python
#
# gfs2am.py - For a given site latitude, longitude, and altitude, this
# script will download an appropriately-subsetted GFS forecast file in
# grib2 format from the NOAA Operational Model Archive Distribution
# System (NOMADS), then generate a corresponding set of am layers
# interpolated to the site position.

import argparse
import datetime
import dateutil.parser as dparser
import math
import pygrib
import requests
from pydap.client import open_url
import numpy as np
import os


PATH = '/Volumes/My_Book_Pro/EHT/GRIB_Files/'


# Numerical and physical constants
BADVAL = -99999. # placeholder data that is missing or not defined on level
BADVAL_TEST = -99998.
G_STD = 9.80665  # standard gravitational acceleration in m / s^2
M_AIR = 28.964   # average dry air mass [g / mole]
M_O3  = 47.997   # O3 mass [g / mole]
H2O_SUPERCOOL_LIMIT = 238. # ice assumed below this temperature [K]
PASCAL_ON_MBAR = 100. # conversion from mbar (hPa) to Pa
R_Earth = 6371* (10**3) #meters

# Ignore H2O above the pressure level defined here.  This is
# needed because the vertical grid used for RH skips the 20 mbar
# level.  In any case, if the small amount of water vapor in
# the stratosphere (about 2 precipitable microns) needs to be
# modeled, the best way is with a climatological profile instead
# of GFS.
H2O_TOP_PLEVEL = 29. # Ignore H2O above this pressure level.

#
# Below are constants and format strings for constructing the data
# request URL.  These include the base URL for the NOMADS CGI interface,
# and various strings for formatting the arguments given to it.  Note
# that some information (e.g. grid, forecast production cycle) gets
# used more than once to construct the CGI request.
#


# GFS variables to be requested, and the format string for adding them
# to the CGI request URL.  The variables are
#   CLWMR - Cloud liquid water mass mixing ratio [kg liquid / kg air]
#   HGT   - Geopotential height [m]
#   O3MR  - Ozone mass mixing ratio [kg O3 / kg air]
#   RH    - Relative Humidity [%]
#   TMP   - Temperature [K]



# Comment string to be printed above model layers:
LAYER_HEADER = """
    #
    # Layer data below were derived from NCEP GFS model data obtained
    # from the NOAA Operational Model Archive Distribution System
    # (NOMADS).  See http://nomads.ncep.noaa.gov for more information.
    #
    #         Production date: {0}
    #                   Cycle: {1:02d} UT
    #                 Product: {2}
    #
    # Interpolated to
    #
    #                latitude: {3} deg. N
    #               longitude: {4} deg. E
    #   Geopotential altitude: {5} m
    #
    """

#
# Function for bilinear grid interpolation.  a[i_lat][i_lon] is a
# 2x2 array of adjacent points on the lon,lat grid.  u,v are the
# fractional distances in grid spacing in grid spacing units from
# the "bottom left" grid point to the interpolation point.
#
def grid_interp(a, u, v):
    return ( a[0][0] * (1.0 - u) * (1.0 - v) + a[1][0] * u * (1.0 - v)
            + a[0][1] * (1.0 - u) * v         + a[1][1] * u * v        )


#
# Parse the command line and validate arguments.
#
parser = argparse.ArgumentParser()
parser.add_argument("lat",      help="site latitude [deg], (-90 to 90)",
                    type=float)
parser.add_argument("lon",      help="site longitude [deg], (-180 to 180)",
                    type=float)
parser.add_argument("altitude", help="site altitude [m]",
                    type=float)
parser.add_argument("gfsdate",  help="GFS production date (YYYYMMDD)",
                    type=str)
parser.add_argument("gfscycle", help="GFS production cycle (0, 6, 12, 18)",
                    type=int)
parser.add_argument("gfsprod",  help="GFS product: anl or f000 - f384 in multiples of 3",
                    type=str)
args = parser.parse_args()

if (args.lat < -90. or args.lat > 90.):
    parser.error("invalid latitude")
if (args.lon < -360. or args.lon > 360.):
    parser.error("invalid longitude")
if (args.lon < 0.):
    args.lon = args.lon + 360.0

if (args.altitude < -500.):
    parser.error("invalid altitude")
try:
    gfsdatetime = dparser.parse(args.gfsdate)
except:
    parser.error("bad GFS production date")
if (gfsdatetime < datetime.datetime(2004, 3, 2)):
    parser.error("GFS production date too early")
if (args.gfscycle not in (0, 6, 12, 18)):
    parser.error("invalid GFS production cycle")
if (args.gfsprod != "anl"):
    if (args.gfsprod[0:1] == "f"):
        forecast_hour = int(args.gfsprod[1:])
        if (forecast_hour < 0 or forecast_hour > 384):
            parser.error("invalid GFS product (hour out of range)")
    else:
        parser.error("invalid GFS product name")

if (args.gfsprod == "anl"):
    
    URL = 'https://nomads.ncdc.noaa.gov/thredds/fileServer/gfs-004-anl'
    
    URL += '/' + str(args.gfsdate[0:-2]) +'/' + str(args.gfsdate) + '/'  + 'gfsanl_4_' + str(args.gfsdate) +'_'
    
    if (args.gfscycle in (0,6)):
        URL += '0' + str(args.gfscycle)+ '00_000.grb2'
    
    if (args.gfscycle in (12, 18)):
        URL +=  str(args.gfscycle)+ '00_000.grb2'
    fname = URL[-31:]

if (args.gfsprod != "anl"):
    if (gfsdatetime < datetime.datetime(2015, 8, 18)):
        parser.error("GFS production date too early")

    URL = 'https://nomads.ncdc.noaa.gov/thredds/fileServer/gfs-004'

    URL += '/' + str(args.gfsdate[0:-2]) +'/' + str(args.gfsdate) + '/'  + 'gfs_4_' + str(args.gfsdate) +'_'
    
    if (args.gfscycle in (0,6)):
        URL += '0' + str(args.gfscycle)+ '00_' + str(args.gfsprod[1:]) +'.grb2'

    if (args.gfscycle in (12, 18)):
        URL +=  str(args.gfscycle)+ '00_' + str(args.gfsprod[1:]) +'.grb2'
    fname = URL[-28:]



if (os.path.isdir(str(PATH)+ '/' + str(args.gfsdate[0:-2])) == False):
    os.system('mkdir ' + str(PATH) +'/'+ str(args.gfsdate[0:-2]))

if (os.path.isdir(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate)) == False):
    os.system('mkdir ' + str(PATH) +'/' + str(args.gfsdate[0:-2]) + '/' + str(args.gfsdate))

if (os.path.isdir(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)) == False):
    os.system('mkdir ' + str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle))

if (os.path.isdir(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle) + '/' + str(args.lat) + 'N' + str(args.lon) + 'E') == False):
    os.system('mkdir ' + str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle) + '/' + str(args.lat) + '_N_' + str(args.lon) + '_E')


'''
if os.path.isfile(str(PATH) + '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsprod) + '/'+ str(fname)) == False:
    
    os.system('sudo wget ' + str(URL) + ' -P ' + str(PATH) + '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsprod) )
'''
################################################################################################################################################################

if os.path.isfile(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_levels'+'.npy') == False:

    grbs = pygrib.open(str(PATH) + '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsprod) + '/'+ str(fname))
    
    LEVELS = []
    geoh_data = []
    tmp_data  = []
    o3m_data = []
    relh_data = []
    clm_data = []
    for grb in grbs:
        if grb.parameterName == 'Temperature' and (grb.level == 10 or grb.level ==  20 or grb.level == 30 or grb.level ==  50 or grb.level ==   70 or grb.level == 100 or grb.level ==  150 or grb.level == 200 or grb.level == 250 or grb.level ==  300 or grb.level == 350 or grb.level == 400 or grb.level == 450 or grb.level ==  500 or grb.level == 550 or grb.level == 600 or grb.level ==  650 or grb.level == 700 or grb.level ==  750 or grb.level == 800 or grb.level ==  850 or grb.level == 900 or grb.level == 925 or grb.level ==  950 or grb.level == 975 or grb.level == 1000):
            if grb.level > 1000:
                continue
            if grb.level < 10:
                continue
            d = (grb.data(lat1 = (args.lat-0.5), lat2 = (args.lat+0.5), lon1 = (args.lon-0.5), lon2 = (args.lon+0.5)))
            tmp_data.append(d[0])
            LEVELS.append(grb.level)
            if LEVELS[0] != 1000:
                if grb.level == 1000:
                    break
        elif grb.parameterName == 'Relative humidity'and (grb.level == 10 or grb.level ==  20 or grb.level == 30 or grb.level ==  50 or grb.level ==   70 or grb.level == 100 or grb.level ==  150 or grb.level == 200 or grb.level == 250 or grb.level ==  300 or grb.level == 350 or grb.level == 400 or grb.level == 450 or grb.level ==  500 or grb.level == 550 or grb.level == 600 or grb.level ==  650 or grb.level == 700 or grb.level ==  750 or grb.level == 800 or grb.level ==  850 or grb.level == 900 or grb.level == 925 or grb.level ==  950 or grb.level == 975 or grb.level == 1000):
            if grb.level > 1000:
                continue
            if grb.level < 10:
                continue
            d = (grb.data(lat1 = (args.lat-0.5), lat2 = (args.lat+0.5), lon1 = (args.lon-0.5), lon2 = (args.lon+0.5)))
            relh_data.append(d[0])
        elif grb.parameterName == 'Geopotential height'and (grb.level == 10 or grb.level ==  20 or grb.level == 30 or grb.level ==  50 or grb.level ==   70 or grb.level == 100 or grb.level ==  150 or grb.level == 200 or grb.level == 250 or grb.level ==  300 or grb.level == 350 or grb.level == 400 or grb.level == 450 or grb.level ==  500 or grb.level == 550 or grb.level == 600 or grb.level ==  650 or grb.level == 700 or grb.level ==  750 or grb.level == 800 or grb.level ==  850 or grb.level == 900 or grb.level == 925 or grb.level ==  950 or grb.level == 975 or grb.level == 1000):
            if grb.level > 1000:
                continue
            if grb.level < 10:
                continue
            d = (grb.data(lat1 = (args.lat-0.5), lat2 = (args.lat+0.5), lon1 = (args.lon-0.5), lon2 = (args.lon+0.5)))
            geoh_data.append(d[0])
        elif grb.name == "Ozone mixing ratio"and (grb.level == 10 or grb.level ==  20 or grb.level == 30 or grb.level ==  50 or grb.level ==   70 or grb.level == 100 or grb.level ==  150 or grb.level == 200 or grb.level == 250 or grb.level ==  300 or grb.level == 350 or grb.level == 400 or grb.level == 450 or grb.level ==  500 or grb.level == 550 or grb.level == 600 or grb.level ==  650 or grb.level == 700 or grb.level ==  750 or grb.level == 800 or grb.level ==  850 or grb.level == 900 or grb.level == 925 or grb.level ==  950 or grb.level == 975 or grb.level == 1000):
            if grb.level > 1000:
                continue
            if grb.level < 10:
                continue
            d = (grb.data(lat1 = (args.lat-0.5), lat2 = (args.lat+0.5), lon1 = (args.lon-0.5), lon2 = (args.lon+0.5)))
            o3m_data.append(d[0])
        elif grb.parameterName == 'Cloud mixing ratio'and (grb.level == 10 or grb.level ==  20 or grb.level == 30 or grb.level ==  50 or grb.level ==   70 or grb.level == 100 or grb.level ==  150 or grb.level == 200 or grb.level == 250 or grb.level ==  300 or grb.level == 350 or grb.level == 400 or grb.level == 450 or grb.level ==  500 or grb.level == 550 or grb.level == 600 or grb.level ==  650 or grb.level == 700 or grb.level ==  750 or grb.level == 800 or grb.level ==  850 or grb.level == 900 or grb.level == 925 or grb.level ==  950 or grb.level == 975 or grb.level == 1000):
            if grb.level > 1000:
                continue
            if grb.level < 10:
                continue
            d = (grb.data(lat1 = (args.lat-0.5), lat2 = (args.lat+0.5), lon1 = (args.lon-0.5), lon2 = (args.lon+0.5)))
            clm_data.append(d[0])

    tmp_data = np.array(tmp_data)
    LEVELS = np.array(LEVELS)
    geoh_data = np.array(geoh_data)
    o3m_data = np.array(o3m_data)
    relh_data = np.array(relh_data)
    clm_data = np.array(clm_data)

    np.save(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_' + str(args.gfsprod)+ '_T'+'.npy',tmp_data)
    np.save(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_RH'+'.npy', relh_data)
    np.save(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_z'+'.npy', geoh_data)
    np.save(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_o3_vmr'+'.npy', o3m_data)
    np.save(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_cloud_mr'+'.npy', clm_data)
    np.save(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_levels'+'.npy', LEVELS)


tmp_data  = np.load(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod) + '_T'+'.npy')
relh_data = np.load(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_RH'+'.npy')
geoh_data = np.load(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_z'+'.npy')
o3m_data  = np.load(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_o3_vmr'+'.npy')
clm_data  = np.load(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_cloud_mr'+'.npy')
LEVELS  = np.load(str(PATH)+ '/' + str(args.gfsdate[0:-2])+ '/' + str(args.gfsdate) + '/' + str(args.gfsdate)+str(args.gfscycle)+ '/' + str(args.lat) + '_N_' + str(args.lon) + '_E' + '/' + str(args.gfsdate)+str(args.gfscycle)+ '_' + str(args.lat) + '_N_' + str(args.lon) + '_E_' + str(args.altitude)+ 'm_'+ str(args.gfsprod)+'_levels'+'.npy')


Pbase = np.zeros(len(LEVELS))
z = np.zeros(len(LEVELS))
T = np.zeros(len(LEVELS))
o3_vmr = np.zeros(len(LEVELS))
RH = np.zeros(len(LEVELS))
cloud_mr = np.zeros(len(LEVELS))


if LEVELS[0] == 1000:
    geoh_data = geoh_data[::-1]
    tmp_data = tmp_data[::-1]
    o3m_data = o3m_data[::-1]
    relh_data = relh_data[::-1]
    clm_data = clm_data[::-1]
    LEVELS = LEVELS[::-1]

for i in np.arange(len(LEVELS)):
    Pbase[i] = LEVELS[i]
    T[i] = np.mean(tmp_data[i])

for i in np.arange(len(geoh_data)):
    z[i] = np.mean(geoh_data[i])

for i in np.arange(len(o3m_data)):
    o3_vmr[i] = np.mean(o3m_data[i])

for i in np.arange(len(relh_data)):
    RH[i] = np.mean(relh_data[i])

for i in np.arange(len(clm_data)):
    cloud_mr[i] = np.mean(clm_data[i])





#######################################################################################################################



# Print out the layer descriptions.  On a layer, mixing ratios and
# RH are set to their averages over the two levels bounding the layer.
if (args.gfsprod == "anl"):
    product_str = "analysis"
else:
    product_str = args.gfsprod[1:] + " hour forecast"
print LAYER_HEADER.format(
                          args.gfsdate,
                          args.gfscycle,
                          product_str,
                          args.lat,
                          args.lon,
                          args.altitude)


for i in np.arange(len(LEVELS)):
    if (z[i] > args.altitude):
        print "layer"
        print "Pbase {0:.1f} mbar  # {1:.1f} m".format(Pbase[i], z[i])
        print "Tbase {0:.1f} K".format(T[i])
        print "column dry_air vmr"
        if (i > 0):
            o3_vmr_mid   = 0.5 * (  o3_vmr[i-1] +   o3_vmr[i])
            RH_mid       = 0.5 * (      RH[i-1] +       RH[i])
            cloud_mr_mid = 0.5 * (cloud_mr[i-1] + cloud_mr[i])
            T_mid        = 0.5 * (       T[i-1] +        T[i])
        else:
            o3_vmr_mid   = o3_vmr[i]
            RH_mid       = RH[i]
            cloud_mr_mid = cloud_mr[i]
            T_mid = T[i]
        if (o3_vmr_mid > 0.0):
            print "column o3 vmr {0:.3e}".format(o3_vmr_mid)
        if (RH_mid > 0.0):
            if (T_mid > H2O_SUPERCOOL_LIMIT):
                print "column h2o RH {0:.2f}%".format(RH_mid)
            else:
                print "column h2o RHi {0:.2f}%".format(RH_mid)
        if (cloud_mr_mid > 0.0):
            # Convert cloud mixing ratio [kg / kg] to cloud total water
            # across the layer [kg / m^2].
            dP = PASCAL_ON_MBAR * (Pbase[i] - Pbase[i-1])
            m = dP / G_STD
            ctw = m * cloud_mr_mid
            if ctw <= 0:
                print ""
                continue
            if (T_mid > H2O_SUPERCOOL_LIMIT):
                print "column lwp_abs_Rayleigh {0:.3e} kg*m^-2".format(ctw)
            
            else:
                print "column iwp_abs_Rayleigh {0:.3e} kg*m^-2".format(ctw)
        print ""



# The base layer and base level of the model are special cases.  First,
# we find the pressure and temperature of the base level by linearly
# interpolating (or extrapolating) log P and T in z.
if (i == 0):
    print "User-specified altitude exceeds top GFS level"
    exit()
if (z[i] == args.altitude): # exact coincidence with model level
    exit()
u = (args.altitude - z[i-1]) / (z[i] - z[i-1])
logP_s = u * math.log(Pbase[i]) + (1.0 - u) * math.log(Pbase[i-1])
P_s = math.exp(logP_s)
T_s = u * T[i] + (1.0 - u) * T[i-1]

if P_s <= Pbase[i]:
    exit()

# Other variables are interpolated or extrapolated linearly in P
# to the base level and clamped at zero.
u = (P_s - Pbase[i-1]) / (Pbase[i] - Pbase[i-1])
o3_vmr_s   =   u * o3_vmr[i] + (1.0 - u) *   o3_vmr[i-1]
RH_s       =       u * RH[i] + (1.0 - u) *       RH[i-1]
cloud_mr_s = u * cloud_mr[i] + (1.0 - u) * cloud_mr[i-1]
if (o3_vmr_s < 0.0):
    o3_vmr_s = 0.0
if (RH_s < 0.0):
    RH_s = 0.0
if (cloud_mr_s < 0.0):
    cloud_mr_s = 0.0
o3_vmr_mid   = 0.5 * (  o3_vmr[i-1] +   o3_vmr_s)
RH_mid       = 0.5 * (      RH[i-1] +       RH_s)
cloud_mr_mid = 0.5 * (cloud_mr[i-1] + cloud_mr_s)
print "layer"
print "Pbase {0:.1f} mbar  # {1:.1f} m".format(P_s, args.altitude)
print "Tbase {0:.1f} K".format(T_s)
print "column dry_air vmr"
if (o3_vmr_mid > 0.0):
    print "column o3 vmr {0:.3e}".format(o3_vmr_mid)
if (RH_mid > 0.0):
    if (T_mid > H2O_SUPERCOOL_LIMIT):
        print "column h2o RH {0:.2f}%".format(RH_mid)
    else:
        print "column h2o RHi {0:.2f}%".format(RH_mid)
if (cloud_mr_mid > 0.0):
    dP = PASCAL_ON_MBAR * (Pbase[i] - Pbase[i-1])
    m = dP / G_STD
    ctw = m * cloud_mr_mid
    if (T_mid > H2O_SUPERCOOL_LIMIT):
        print "column lwp_abs_Rayleigh {0:.3e} kg*m^-2".format(ctw)
    else:
        print "column iwp_abs_Rayleigh {0:.3e} kg*m^-2".format(ctw)
exit()

