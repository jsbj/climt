import climt
import pylab as pl
import numpy as np

pl.ioff()

r=climt.radiation(scheme='cam3')
r(co2=300., o3=r['p']*0.+1.e-9, q=r['p']*0.+1.e-9)
print r['LwToa']+r['Ts']**4*climt.Parameters()['stebol']
pl.plot(r['lwdflx'],r['p'])

r=climt.radiation(scheme='chou')
r(co2=300., o3=r['p']*0.+1.e-9, q=r['p']*0.+1.e-9)
print r['LwToa']+r['Ts']**4*climt.Parameters()['stebol']
pl.plot(r['lwflx']-r['LwToa'],r['p'])

r=climt.radiation(scheme='rrtm')
r(co2=300., p=r['p'], lev=r['lev'], ps=r['ps'], o3=r['p']*0.+1.e-9, q=r['p']*0.+1.e-9)
print r['LwToa']+r['Ts']**4*climt.Parameters()['stebol']
f = r['lwdflx'][::-1]
f = (f[1:]+f[:-1])/2.
pl.plot(f,r['p'])

pl.legend(['cam3','chou','rrtm'])

pl.ylim([1000,0])
pl.show()
