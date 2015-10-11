import numpy as np
import pylab as pl

from readmock import *
from load_oiidata import *
from read_photoion import *
from get_lf import *
from astropy.cosmology import FlatLambdaCDM 
import multiprocessing as mp

def compute_oiilines(lineinfo,linesarr,qpar,zcold,llyc,lname):
	 l1 = integ_line(lineinfo,linesarr,qpar,zcold,llyc,'OII_3727')
	 l2 =  integ_line(lineinfo,linesarr,qpar,zcold,llyc,'OII_3729')
	 loii = np.log10(10**l1 + 10**l2)
	 return loii


cosmo = FlatLambdaCDM(H0=73, Om0=0.25)

zrArr = comparat_zranges()

zminlist = map(float,zrArr['zmin'])
zmaxlist = map(float,zrArr['zmax'])
minz = np.min(zminlist)
maxz = np.max(zmaxlist)
#maxz = np.max(zrArr['zmax'])

print minz,maxz
minsfr = 10.0

print 'loading lightcone data'
GalArr,lLyc,ngg,DistNz,zspace = readmock(zmin=minz,zmax=maxz, sfrmin=minsfr,zspace=1)

print 'now loading photo-ionisation grid'
lineinfo,linesarr = read_photoion()

qpar = qZrelation(GalArr['Zcold'])  # assuming default pars.


#loii_1 = np.zeros(ngg,dtype=np.float32)
#loii_2 = np.zeros(ngg,dtype=np.float32)

print 'Computing emission line(s)'

def get_all_oii(i):
#for i in range(ngg):
	loii_1 = integ_line(lineinfo,linesarr,qpar[i],GalArr['Zcold'][i],lLyc[i],'OII_3727') 
	loii_2 = integ_line(lineinfo,linesarr,qpar[i],GalArr['Zcold'][i],lLyc[i],'OII_3729') 
	if i == ngg/10:
		print '.10%'
	if i == ngg/5:
		print '..20%'
	if i == ngg/2:
		print '...50%'
	if i == ngg*0.9:
		print '....90%'
	return np.log10(10**loii_1 + 10**loii_2)

pool = mp.Pool(processes = 8)
loii = [pool.map(get_all_oii, range(ngg))]

#	if i == ngg/10:
#		print '.10%'
#	if i == ngg/5:
#		print '..20%'
#	if i == ngg/2:
#		print '...50%'

#loii = np.log10(10**loii_1 + 10**loii_2)
print 'emission line calculated for all galaxies'

# simple attenuation to match observed data.

att 	= [0.1,1.0]
zratt = [0.2,1.6]
att_gal = np.interp(GalArr['redshift'],zratt,att)
loii = loii + np.log10(att_gal)


h = 0.73
facl = -2*np.log10(h)
facn = +3*np.log10(h)

nLfs = len(zrArr['zmean'])


# loop over each Obs. plot

font = {'family' : 'TeX Gyre Heros',
	        'weight' : 'normal',
	        'size'   : 6}

pl.rc('font', **font)

pl.figure(1)
for i in range(nLfs):

	data = oiidata('comparat',zrArr['zmean'][i])
	nbb = len(data['Lum'])

	d1 = cosmo.comoving_distance(zmaxlist[i])
	d0 = cosmo.comoving_distance(zminlist[i])
	Vol_z = 4./3*np.pi*(d1.value**3 - d0.value**3)/8.0

	izr = np.where((GalArr['redshift'] < zmaxlist[i]) &
	(GalArr['redshift'] >= zminlist[i]))

	loii_i = loii[0][izr]

	print 'computing mock LF for redshift range ',zrArr['zmean'][i]
	lf0 = get_lf(loii_i+40.0,Vol=Vol_z,binsize=0.1,
	minsample=42.0)


	pl.subplot(240 + (i+1),xlim=[40,45.9],ylim=[-10,-1])

	pl.plot(data['Lum']+facl,np.log10(data['dn_dL'])+facn,'-o',markersize=3,
	label='Comparat+14, z='+np.str(zrArr['zmean'][i]),color='0.75')

	errhi = np.log10(data['err'] + data['dn_dL']) - np.log10(data['dn_dL'])
	errlo = - np.log10( - data['err'] + data['dn_dL']) + np.log10(data['dn_dL'])

	pl.errorbar(data['Lum']+facl,np.log10(data['dn_dL'])+facn,[errlo,errhi],
	color='0.75',ecolor='0.75')

#pdata=pl.errorbar(xc+facl,yc+facn,yerr = (yc - (np.log10(10**yc-ec)),(np.log10(ec + 10**yc)-yc)),
#ecolor='0.75',linewidth=2,fmt=None)

	pl.plot(lf0['Lum'],lf0['dn_dL'],color='Green',linewidth=2)

	pl.xlabel(r'$\log(L_{[O II]\lambda 3727}[{\rm erg \ s^{-1} h^{-2}}])$',fontsize=8)
	if i == 0 or i == 4:
		pl.ylabel(r'$\log({\rm d}n/{\rm d}L(\Delta \log(L))^{-1} [{\rm (Mpc/h)^{-3}}])$',fontsize=8)

	pl.legend(loc=3,fontsize=6)



pl.savefig('oii_lfs.png',bbox_inches='tight',dpi=350)
pl.close("all")
