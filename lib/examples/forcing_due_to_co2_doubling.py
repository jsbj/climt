#!/usr/bin/env python

#
# Set up realistic tropical temperature and moisture profiles
# and compute radiative fluxes
#

from pylab import *
import climt
# import PyTran
ioff()
r = climt.radiation()


globalTs = 299.7
globalT = [
293.7,
287.7,
283.7,
277.0,
270.3,
263.6,
257.0,
250.3,
243.6,
237.0,
230.1,
223.6,
217.0,
210.3,
203.7,
197.0,
194.8,
198.8,
202.7,
206.7,
210.7,
214.6,
217.0,
219.2,
221.4,
232.3,
243.1,
254.0,
264.8,
270.2,
218.9,
190.7]



globalT.reverse()






































def getModtranOlrVsCo2():
  co2 = array([280.*2.**i for i in arange(-1,12)])
  olr = array([361.7, 357.6, 353.3, 348.9, 344.1, 339.4, 333.8, 327.2, 319.7, 310.7, 300.5, 289.0, 276.6]) # tropical, No O3, no H2O
  #olr = array([355.448, 351.366, 347.284, 342.888, 338.492, 333.782, 328.444, 322.478, 315.256]) # tropical, with O3, no H2O
  return co2,olr
    
def getModtranProfile(isothermal_stratosphere=True):
  if isothermal_stratosphere:
    jim = """\
    1013.000  299.7
     904.000  293.7
     805.000  287.7
     715.000  283.7
     633.000  277.0
     559.000  270.3
     492.000  263.6
     432.000  257.0
     378.000  250.3
     329.000  243.6
     286.000  237.0
     247.000  230.1
     213.000  223.6
     182.000  217.0
     156.000  210.3
     132.000  203.7
     111.000  197.0
      93.700  197.0
      78.900  197.0
      66.600  197.0
      56.500  197.0
      48.000  197.0
      40.900  197.0
      35.000  197.0
      30.000  197.0
      25.700  197.0
      12.200  197.0
       6.000  197.0
       3.050  197.0
       1.590  197.0
       0.854  197.0
       0.058  197.0
       0.000  197.0"""

  else:
    jim = """\
    1013.000  299.7
     904.000  293.7
     805.000  287.7
     715.000  283.7
     633.000  277.0
     559.000  270.3
     492.000  263.6
     432.000  257.0
     378.000  250.3
     329.000  243.6
     286.000  237.0
     247.000  230.1
     213.000  223.6
     182.000  217.0
     156.000  210.3
     132.000  203.7
     111.000  197.0
      93.700  194.8
      78.900  198.8
      66.600  202.7
      56.500  206.7
      48.000  210.7
      40.900  214.6
      35.000  217.0
      30.000  219.2
      25.700  221.4
      12.200  232.3
       6.000  243.1
       3.050  254.0
       1.590  264.8
       0.854  270.2
       0.058  218.9
       0.000  190.7"""
  values = [float(i) for i in jim.split()]
  p,T = array(values[::2]), array(values[1::2])
  p = p[::-1]
  T = T[::-1]
  ps = p[-1]
  Ts = T[-1]
  dp = p[1:]-p[:-1]
  p = (p[1:]+p[:-1])/2.
  T = (T[1:]+T[:-1])/2.
  return ps,p,dp,T,Ts

def getProfiles(Ts,rh=0.4,pcld=800.):
  #--- initialise T,q
  # Surface temperature
  #Ts = 273.15 + 30.                         
  # Strospheric temp
  Tst = 273.15 - 80.                         
  # Surface pressure
  ps = 1000.
  # Equispaced pressure levels
  p = ( arange(r.nlev)+ 0.5 )/r.nlev * ps
  # Return moist adiabat with given rel hum
  if rh == 0.: rh=1.e-21
  (T,q) = climt.thermodyn.moistadiabat(p, Ts, Tst, rh)
  cldf = q*0. 
  clwp = q*0.
  kcld = argmin(abs(p-pcld))
  #cldf[kcld] = 1.
  #clwp[kcld] = 100.
  z = -climt.Parameters()['Rd']*T/9.8*log(p/ps)
  dp = p*0. + p[1]-p[0]
  return T,q,cldf,clwp,p,dp,ps,z

def getOlrVsCo2(scheme='cam3',profile='modtran_tropical', isothermal_stratosphere=False):
  if profile == 'modtran_tropical':
    ps,p,dp,T,Ts = getModtranProfile(isothermal_stratosphere)
    cldf = p*0.
    clwp = p*0.
    q = p*0.+1.e-19
  else:
    Ts = 273. + 10 
    T,q,cldf,clwp,p,dp,ps,z = getProfiles(Ts,rh=0.)
  co2 = 280.*2.**arange(-1,12)
  olr=[]
  if scheme == 'pytran':
    r = PyTran.RadiativeTransfer(dWave=0.02)
    for c in co2:
      q = c * 1.e-6 + p*0.
      olr.append( r.olr(p,dp,T,Ts,q,r.wn).sum()*r.dWave )
      print olr[-1]
  else:
    r = climt.radiation(scheme=scheme,do_sw=0)
    for c in co2:
      r(p=p, ps=ps, T=T, Ts=Ts, q=q, cldf=cldf, clwp=clwp, co2=c, o3=(p*0+1.e-19))
      olr.append(-r['LwToa'])      
  olr=array(olr)
  return co2,olr

def getIPCC(co2,form=1):
  if form==1:
    alpha=5.35
    DeltaF = alpha*log(co2[1:]/co2[0:-1]) 
  if form==2:
    alpha = 4.841
    beta = 0.0906
    DeltaF = alpha*log(co2[1:]/co2[0:-1]) + beta*(co2[1:]**.5 - co2[0:-1]**.5)
  if form==3:
    def g(c): return log(1. + 1.2*c + 0.005*c**2 + 1.4e-6*c**3)
    alpha = 3.35
    DeltaF = alpha*(g(co2[1:]) - g(co2[0:-1]))
  return DeltaF

if __name__ == "__main__":
  import matplotlib.pyplot as pyplot
  profile = 'modtran_tropical'
  fig = pyplot.figure()
  fig_left = fig.add_subplot(1,2,1)
  fig_right = fig.add_subplot(1,2,2)

  # co2,olr = getOlrVsCo2(scheme='rrtm',profile=profile,isothermal_stratosphere=True)
  # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'g-s')
  # fig_right.semilogx(co2,olr,'g-s')
  
  co2,olr = getOlrVsCo2(scheme='rrtm',profile=profile)
  # import pdb; pdb.set_trace()
  # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'g-o')
  # fig_right.semilogx(co2,olr,'g-o')
  r = climt.radiation(scheme = 'rrtm')
  olr = []
  for c in co2:
      r(Ts = 400, T = 400 * (r['p']/1000.) ** (2/7.), co2 = c)
      olr.append(r['LwToa'])
  olr = array(olr)
  fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'g-s')
  fig_right.semilogx(co2,abs(olr),'g-s')
  z = climt.radiation(scheme = 'cam3')
  olr = []
  for c in co2:
      c = array(c)
      print(c)
      z(Ts = 400, T = 400 * (z['p']/1000.) ** (2/7.), co2 = c)
      olr.append(z['LwToa'])
  olr = array(olr)
  fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'r-s')
  fig_right.semilogx(co2,abs(olr),'r-s')
  print(str(globalT))
  print(str(globalTs))

  # 
  # co2,olr = getOlrVsCo2(scheme='cam3',profile=profile)
  # # import pdb; pdb.set_trace()
  # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'b-o')
  # fig_right.semilogx(co2,olr,'b-o')
  # # co2,olr = getOlrVsCo2(scheme='cam3',profile=profile,isothermal_stratosphere=True)
  # # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'b-s')
  # # fig_right.semilogx(co2,olr,'b-s')
  # 
  # co2,olr = getOlrVsCo2(scheme='chou',profile=profile)
  # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'r-o')
  # fig_right.semilogx(co2,olr,'r-o')
  # 
  # # co2,olr = getOlrVsCo2(scheme='chou',profile=profile,isothermal_stratosphere=True)
  # # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'r-s')
  # # fig_right.semilogx(co2,olr,'r-s')
  # 
  # 
  # #semilogx(co2[:-1],getIPCC(co2,form=1),'y-o')
  # # semilogx(co2[:-1],getIPCC(co2,form=2),'y-s')
  # #semilogx(co2[:-1],getIPCC(co2,form=3),'y-^')
  # 
  # #co2,olr = getOlrVsCo2(scheme='pytran',profile=profile)
  # #semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'g-o')
  # 
  # co2,olr = getModtranOlrVsCo2()
  # fig_left.semilogx(co2[:-1],abs(olr[1:]-olr[:-1]),'k-o')
  # fig_right.semilogx(co2,olr,'k-o')
  
  
  
  fig.show()
