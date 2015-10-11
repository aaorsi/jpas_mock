import numpy as np
import pylab as pl

from readmock import *
from load_laedata import *

GalArr,Lya0,NGals = readmock()

h = 0.73
facl = -2*np.log10(h)
facn = +3*np.log10(h)

hayes = laedata('hayes_z2')
blanc = laedata('blanc_z19_38')
ciardullo = laedata('ciardullo')

gronwall = laedata('gronwall_z3_cum')
ouchi	 = laedata('ouchi_z3_cum')

zb = 0.02
zrange2 = [2.2-zb/2.,2.2+zb/2.]
zrange3 = [3.09,3.13]

izr2 = [(GalArr['redshift'] < zrange2[1]) and \
	  (GalArr['redshift'] > zrange2[0])]

d1 = np.max(GalArr['distance'][izr2])
d0 = np.min(GalArr['distance'][izr2])

Volume_z2 = 4./3*np.pi*(d1**3 - d0**3)/8.0
Lyaz2 = Lya0[izr2]


#facl = facn = 0

print 'computing mock LF in redshiftrange ',zrange2
#lf0_z2 = lae.get_lf(np.log10(Lyaz2)+40.0,Volume=Volume_z2,binsize=0.1,minsample=42.0)


lf0_z2 = np.loadtxt('LF_mock0',dtype=lftype)

izr3 = np.where(GalArr['redshift'] < zrange3[1]) and \
	  np.where(GalArr['redshift'] > zrange3[0])

d1 = np.max(GalArr['distance'][izr3])
d0 = np.min(GalArr['distance'][izr3])

Volume_z3 = 4./3*np.pi*(d1**3 - d0**3)/8.0
Lyaz3 = Lya0[izr3]


#facl = facn = 0

print 'computing mock LF in redshift range ',zrange3
lf0_z3 = get_lf(np.log10(Lyaz3)+40.0,Volume=Volume_z3,binsize=0.1,minsample=42.0)

rev_cum = np.cumsum(10**lf0_z3['dn_dL'][::-1])*0.1
cum = np.log10(rev_cum[::-1])/0.1

nbb = len(lf0_z3['Lum'])
cum = np.zeros(nbb)
diff = 10**lf0_z3['dn_dL']
for i in range(nbb):
	cum[i] = np.log10(np.sum(diff[i:]))


print 'differential arr:',np.log10(diff)
print 'cummulative arr:',cum


xh = hayes['Lum'] 
yh = hayes['dn_dL'] 
eh = hayes['err'] 

xc = ciardullo['Lum']
yc = ciardullo['dn_dL']
ec = ciardullo['err'] 

xg = gronwall['Lum']
yg = gronwall['dn_dL']
eg = gronwall['err']

xo = ouchi['Lum']
yo = ouchi['dn_dL']
eo = ouchi['err']

#pl.plot(xh+facl,yh+facn,'gD')
#pl.errorbar(xh+facl,yh+facn,yerr = (yh - (np.log10(10**yh-eh)),(np.log10(eh + 10**yh)-yh)))
#pl.plot(blanc['Lum'] + facl,blanc['dn_dL'] +facn,'bx')

pl.figure(1)

pl.subplot(121)

pl.plot(xc+facl,yc+facn,'o',markersize=10,label='Ciardullo et al. (2011)',color='0.75')
pdata=pl.errorbar(xc+facl,yc+facn,yerr = (yc - (np.log10(10**yc-ec)),(np.log10(ec + 10**yc)-yc)),
ecolor='0.75',linewidth=2,fmt=None)

pl.plot(lf0_z2['Lum'],lf0_z2['dn_dL'],label=r'$f_{\rm esc}(Ly\alpha) = 1$',color='Green',
linewidth=2)

pl.xlabel(r'$\log(L_{Ly\alpha}[{\rm erg \ s^{-1} h^{-2}}])$',fontsize=15)
pl.ylabel(r'$\log({\rm d}n/{\rm d}L_{Ly\alpha} (\Delta \log(L_{Ly\alpha}))^{-1} [{\rm (Mpc/h)^{-3}}])$',fontsize=15)

fesc = 0.08
lflya_z2 = lf0_z2['Lum'] + np.log10(fesc)

pmock=pl.plot(lflya_z2,lf0_z2['dn_dL'],'-b',label=r'$f_{esc}(Ly\alpha)=0.08$',linewidth=2)
pl.ylim([-10,-2])
pl.xlim([41.5,46])
#pl.axes().set_aspect('equal','datalim')

pl.legend(loc=3,fontsize=10)
pl.text(41.7,-8,r'$z=2.2 \pm 0.01$',fontsize=15)

pl.subplot(122)

pl.plot(xg+facl,yg+facn,'o',markersize=10,label='Gronwall et al. (2007)',color='0.75')
pdata = pl.errorbar(xg+facl,yg+facn,yerr = (yg - (np.log10(10**yg-eg)),(np.log10(eg + 10**yg)-yg)),
ecolor='0.75',linewidth=2,fmt=None)

pl.plot(xo+facl,yo+facn,'*',markersize=10,label='Ouchi et al. (2008)',color='1')
pdata = pl.errorbar(xo+facl,yo+facn,yerr = (yo - (np.log10(10**yo-eo)),(np.log10(eo + 10**yo)-yo)),
ecolor='1',linewidth=2,fmt=None)

pl.plot(lf0_z3['Lum'],cum,label=r'$f_{\rm esc}(Ly\alpha) = 1$',color='Green',
linewidth=2)

pl.xlabel(r'$\log(L_{Ly\alpha}[{\rm erg \ s^{-1} h^{-2}}])$',fontsize=15)
pl.ylabel(r'$\log(n>L_{Ly\alpha}[{\rm (Mpc/h)^{-3}}])$',fontsize=15)

fesc = 0.5
lflya_z3 = lf0_z3['Lum'] + np.log10(fesc)

pmock=pl.plot(lflya_z3,cum,'-b',label=r'$f_{esc}(Ly\alpha)=0.5$',linewidth=2)
pl.ylim([-10,-2])
pl.xlim([41.5,46])
#pl.axes().set_aspect('equal','datalim')

pl.legend(loc=3,fontsize=10)
pl.text(41.7,-8,r'$z=3.06 \pm 0.07$',fontsize=15)


