#!/usr/bin/env python

from numpy import *
import climt

#--- initialise T,q
# Surface temperature
Ts = 273.15 + 30.                         
# Strospheric temp
Tst = 273.15 - 80.                         
# Surface pressure
ps = 1000.
# Equispaced pressure levels
nlev = climt.get_nlev()
p = ( arange(nlev)+ 0.5 )/nlev * ps
# Return moist adiabat with 70% rel hum
T = nlev * [Ts]
q = nlev * [0.]
o3 = nlev * [0.]

cldf = q
# cldf[len(cldf)/3] = 0. # 0.5

ciwp = q
# ciwp[len(cldf)/3] = 10.

r_liq = q

for scheme in ['rrtm']:
    r = climt.radiation(scheme=scheme, o3=o3, T=T, p=p, ps=ps, q=q, Ts=Ts, co2=400)
    # r(p=p, ps=ps, T=T, Ts=Ts, q=q, cldf=cldf, ciwp=ciwp, r_liq=r_liq, o3 = o3)
    print '\n%s' % scheme
    print r['SwToa'],r['LwToa'] #,r['SwSrf'],r['LwSrf']
    # print r['SwToaCf'],r['LwToaCf'],(r['solin']-r['SwToa'])/r['solin'],r['asdir']

    # print r['ciwp']
