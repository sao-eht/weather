from __future__ import division
from __future__ import print_function

import numpy as np
import ehtim as eh
from   ehtim.calibrating import self_cal as sc
from jd_calc import jd_calc


def get_elev(site,dt,mjd):

    # Load the image and the array
    im = eh.image.load_txt('avery_sgra_eofn.txt')
    im.mjd = mjd
    eht = eh.array.load_txt('SITES.txt')

    tint_sec = dt*3600
    tadv_sec = dt*3600
    tstart_hr = 0.
    tstop_hr = 24.
    bw_hz = 4e9
    obs = im.observe(eht, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz,
                     sgrscat=False, ampcal=True, phasecal=False)

    #sites = ['PDB','PV','SMT','SMA','LMT','ALMA','SPT']
    #for site in sites:
    #    plotdata = obs.unpack_bl(site, site, 'el1', ang_unit='deg', debias=True)
    #    plt.plot(plotdata['time'][:,0], plotdata['el1'][:,0],label=site)
    #plt.legend()

    plotdata = obs.unpack_bl(site, site, 'el1', ang_unit='deg', debias=True)
    elev = [plotdata['time'][:,0], plotdata['el1'][:,0]]
    elev[0]=[elev[0][i] for i in range(0,len(elev[0]),2)]
    elev[1]=[elev[1][i] for i in range(0,len(elev[1]),2)]

    return elev


def record_elev_apr_2018(site):
    f1 = open('elev_%s_apr_2018.txt'%(site),'w')
    f2 = open('elev_apr_2018_time.txt','w')
    #for m in range(1,13):
    for m in range(4,5):
        #if m==1: me=31
        #if m==2: me=28
        #if m==3: me=31
        if m==4: me=30
        #if m==5: me=31
        #if m==6: me=30
        #if m==7: me=31
        #if m==8: me=31
        #if m==9: me=30
        #if m==10: me=31
        #if m==11: me=30
        #if m==12: me=31
        for d in range(1,me+1):
            print('%s: %d/%d'%(site,m,d))
            mjd = jd_calc( 2018, m, d )
            elev = get_elev(site,0.1,mjd)
            #if m*d==1:
            if m*d==4:
                for t in elev[0]: f2.write('%f '%t)
                f2.close()
            f1.write('2018 %d %d %d '%(m,d,mjd))
            for el in elev[1]: f1.write('%f '%el)
            f1.write('\n')
    f1.close()


def record_elev_2018(site):
    f1 = open('elev_%s_2018.txt'%(site),'w')
    f2 = open('elev_2018_time.txt','w')
    for m in range(1,13):
        if m==1: me=31
        if m==2: me=28
        if m==3: me=31
        if m==4: me=30
        if m==5: me=31
        if m==6: me=30
        if m==7: me=31
        if m==8: me=31
        if m==9: me=30
        if m==10: me=31
        if m==11: me=30
        if m==12: me=31
        for d in range(1,me+1):
            print('%s: %d/%d'%(site,m,d))
            mjd = jd_calc( 2018, m, d )
            elev = get_elev(site,0.1,mjd)
            if m*d==1:
                for t in elev[0]: f2.write('%f '%t)
                f2.close()
            f1.write('3018 %d %d %d '%(m,d,mjd))
            for el in elev[1]: f1.write('%f '%el)
            f1.write('\n')
    f1.close()

#record_elev_apr_2018('LMT')
#record_elev_apr_2018('SMT')
#record_elev_apr_2018('SMA')
#record_elev_apr_2018('ALMA')
#record_elev_apr_2018('SPT')
#record_elev_apr_2018('PV')
#record_elev_apr_2018('PDB')

#record_elev_2018('LMT')
record_elev_2018('SMT')
#record_elev_2018('SMA')
#record_elev_2018('ALMA')
#record_elev_2018('SPT')
#record_elev_2018('PV')
#record_elev_2018('PDB')
