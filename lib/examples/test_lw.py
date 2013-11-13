import climt
import pylab as pl
import numpy as np
from numpy import arange

pl.ioff()

co2 = 0.

#--- instantiate radiation module
r = climt.radiation(scheme='cam3')

#--- initialise T,q
# Surface temperature
Ts = 273.15 + 30.                         
# Strospheric temp
Tst = 273.15 - 80.                         
# Surface pressure
ps = 1000.
# Equispaced pressure levels
p = ( arange(r.nlev)+ 0.5 )/r.nlev * ps
# Return moist adiabat with 70% rel hum
(T,q) = climt.thermodyn.moistadiabat(p, Ts, Tst, 1.)


cam3 = climt.radiation(scheme='cam3')
cam3(co2=co2, p=p, ps=ps, T=T, Ts=Ts, q=q)
# pl.plot(cam3['lwdflx'], cam3['p'])

chou = climt.radiation(scheme='chou')
chou(co2=co2, p=p, ps=ps, T=T, Ts=Ts, q=q)
# pl.plot(chou['lwflx'] - chou['LwToa'], chou['p'])

rrtm = climt.radiation(scheme='rrtm')
rrtm(co2=co2, p=p, ps=ps, T=T, Ts=Ts, q=q * 0.)
f = rrtm['lwdflx']
f = (f[1:]+f[:-1])/2.
# pl.plot(f,rrtm['p'])
# 
# pl.legend(['cam3','chou','rrtm'])
# 
# pl.ylim([1000,0])
# pl.title('no CO2, H2O at 1%')
# pl.show()

pl.plot([(b[1]+b[0])/2. for b in rrtm['lwbands']], rrtm['lwuflxband'][0], 'b-')
pl.plot([(b[1]+b[0])/2. for b in rrtm['lwbands']], rrtm['lwuflxband'][-1], 'r-')
pl.title('OLR')
pl.show()