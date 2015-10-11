"Script to calculate magnitudes for LAE sample, given their redshifts and Ly-alpha luminosity"

import numpy as np
import pylab as pl
import weighted as wp
import scipy.integrate as integral
import read_filters as rf
from get_lya import *
import get_jpaslims as jp
from astropy.io import ascii

def find_nearest_vector(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return str(f)[:slen]

Make_Figure = False	#Creates a plot with the LFs at different bands
Make_LF		= False
Make_Output = True	#Outputs galaxy catalogues that can be used for clustering analysis.
Make_Random	= True	
Make_FigCompareMags = False
CheckMock = False

zfesc = [2.2,3.1]
fescfix = [0.08,0.3]

zb = 0.02
zrange = [2.2-zb/2.,2.2+zb/2.]

zrange = [1.5,6.0]

nu_Lya0 = 2.47e15  # Hertz
L_Lya0  = 1215.67  # Angstroms

Lyaz,izr,GalArr,Vol = get_lya(zrange)

ngals = len(GalArr)
fesc_gal = np.interp(GalArr['redshift'],zfesc,fescfix)
for i in np.arange(ngals):
	if fesc_gal[i] > 1.0 :
		fesc_gal[i] = 1.0
	if fesc_gal[i] < 0.0:
		fesc_gal[i] = 0.0
	if GalArr['redshift'][i] < zrange[0]:
		fesc_gal[i] = fescfix[0]


FiltArr,wrange = rf.read_filters()
jpaslims = jp.get_jpaslims()

app_type = 1	# 0 means per 3 arcsec diameter, 1 means per square arcsec


nFil = len(FiltArr['w'])
CentralW = np.zeros(nFil)
int_ab = np.zeros(nFil)
p10p90W = np.zeros([nFil,2])

for i in range(nFil):
	CentralW[i] = wp.quantile(FiltArr['w'][i],FiltArr['t'][i],0.5)
	nz = np.where(FiltArr['w'][i] != 0)
	int_ab[i] = (integral.simps(FiltArr['t'][i][nz]/FiltArr['w'][i][nz],FiltArr['w'][i][nz]))
	p10p90W[i,0] = wp.quantile(FiltArr['w'][i][nz],FiltArr['t'][i][nz],0.1)
	p10p90W[i,1] = wp.quantile(FiltArr['w'][i][nz],FiltArr['t'][i][nz],0.9)
	print 'int_ab[i]',int_ab[i]

mpc2cm_24 = 3.0856778 # mpc-to-cm /1e24
distance = np.sqrt(GalArr['pos'][:,0]**2 + GalArr['pos'][:,1]**2+GalArr['pos'][:,2]**2)

log_lya0 = np.log10(Lyaz) + 40.0 + np.log10(fesc_gal)

for i in np.arange(ngals):
	if np.isinf(log_lya0[i]):
		print 'gal ',i,' lya luminosity',Lyaz[i]

log_dist = np.log10(distance) + np.log10(mpc2cm_24)+24.0 + np.log10(1+GalArr['redshift'])

log_Flya0 = log_lya0 - (np.log10(4*np.pi)+2*(log_dist))
FLya0 = 10**log_Flya0
fnu = FLya0/(nu_Lya0/(1+GalArr['redshift']))
wgal = L_Lya0*(1 + GalArr['redshift'])
flam = FLya0/wgal

mag_app = np.zeros(ngals)
id_filter = np.zeros(ngals)
mag2_app = np.zeros(ngals)

for i in np.arange(ngals):
	idFilter,Val = find_nearest_vector(CentralW,wgal[i])
	id_filter[i] = idFilter

print 'id filters:',id_filter[0:100]

minF = np.min(id_filter)
maxF = np.max(id_filter)

lfbin = 0.25
lfmax = 26.0
lfmin = 15.0
nlbins = (lfmax - lfmin)/(lfbin) + 1.0

allmag_lf = np.zeros(nFil,dtype=lftype*int(nlbins))

#for i in range(ngals):
print 'computing narrow-band magnitudes for each galaxy'
print 'tlam0 wgal int_ab log_FLya0'
for i in np.arange(ngals):
	idFilter = id_filter[i]
	nz = np.where(FiltArr['w'][idFilter] != 0)
	lamb = FiltArr['w'][idFilter][nz]
	tlamb = FiltArr['t'][idFilter][nz]
	tlam0 = np.interp(wgal[i],lamb,tlamb)
	mag_app[i] = -2.5*(log_Flya0[i]+ np.log10((wgal[i]*tlam0)/int_ab[idFilter]) - 18-np.log10(3.0)) - 48.60
		

if Make_LF:
	nfilters = 0
	for j in np.arange(nFil):
		id_filter_i = np.where((id_filter == j) & 
					  (~np.isinf(mag_app)))
		gals_filt = len(id_filter_i[0])
		if gals_filt == 0:
			continue
		
		allmag_lf[j] = get_lf(mag_app[id_filter_i],
		Volume = Vol, binsize=lfbin,minsample = lfmin,maxsample=lfmax)
		print 'LF computed for filter ',j
		nfilters += 1


if Make_Figure:
	pl.figure()
	pl.subplots_adjust(hspace=0,wspace=0)
	xtloc = [21,23,25,27]
	for ip in np.arange(nfilters):
		ax = pl.subplot(3,np.ceil(nfilters/3.),ip+1,xlim=[20,24],ylim=[-12,-4],aspect=1.75)
		if (ip ==0) or (ip == np.ceil(nfilters/3)+1) or (ip == 2*(np.ceil(nfilters/3)+1)):
			pl.setp(ax.get_yticklabels(), visible=True)	
			pl.ylabel(r'$\log({\rm d}n/{\rm d}m \ (\Delta m)^{-1} [{\rm (Mpc/h)^{-3}}])$',fontsize=5)
		else:
			pl.setp(ax.get_yticklabels(), visible=False) 

		pl.xticks(xtloc)
		pl.plot(allmag_lf['Lum'][ip],allmag_lf['dn_dL'][ip],'-b')
		pl.plot([jpaslims[ip][app_type],jpaslims[ip][app_type]],
		[-12,-4],'--',linewidth=0.75,color='black')
		pl.fill_between,([19,jpaslims[ip][app_type]],[-12,-12],[-4,-4],'r')
		flabel = str(wrange[ip])
		pl.xlabel(r'$m_{'+flabel+'}$',fontsize=5)
		pl.tick_params(axis='both',labelsize=6)	
		zmed = trunc(CentralW[ip]/L_Lya0 - 1,2)
		zminus = trunc(CentralW[ip]/L_Lya0 -1 - (p10p90W[ip,0]/L_Lya0 - 1),2)
		zplus  = trunc(p10p90W[ip,1]/L_Lya0 - 1 - (CentralW[ip]/L_Lya0 - 1) ,2)
		pl.text(20.,-5,r'$z='+zmed+'^{+'+zplus+'}_{-'+zminus+'}$',fontsize=3,
		bbox=dict(facecolor='0.85', alpha=0.5)) 
	
	pl.savefig('../plots/test.png',bbox_inches='tight',dpi=750)	


galslim = np.zeros(ngals) 
if Make_Output:
	RandomFrac = 10.0	# N_random = RandomFrac * Ngal
	OutputDir = '/home/aaorsi/work/jpas_mock/out/'
	FileAll   = OutputDir + 'jpas_laes_all'	# A file with ALL LAEs in jpas mocks
	RandomFileAll = FileAll +'.random'
#	f = open(OutputDir + FileAll,'w')

	# Catalogue containing all Ly-alpha emitters:
	

	for iM in np.arange(1):
		k = 0
		done = 0
		for ip in np.arange(nFil):
			gals_ip = np.where(( id_filter == ip) &
					(mag_app <= (iM * 3) +jpaslims[ip][app_type]) & (mag_app > 0))
			ng = len(gals_ip[0])
			
			if ((ng > 0) & (done == 0)):
				minfilter = ip
				done = 1		
	
			if ((ng == 0) & (done == 1)):
				maxfilter = ip
				done = 2
	
			galslim[k:k+ng] = gals_ip[0]
			k += ng
			print 'number of galaxies in filter '+np.str(ip)+':'+str(ng)

		print 'sel '+str(iM)


		nng = list(galslim[:k])
#	ascii.write([GalArr['pos'][nng], GalArr['redshift'][nng]],f,formats={'%12.1f','%12.1f','%12.1f','%12.1f'})	
		X = np.zeros((k,4))
		X[:,0:3] = GalArr['pos'][nng]
		X[:,3] = GalArr['redshift'][nng]
		ascii.write(X,FileAll)
	
#	Here the selection function is computed

		print 'minfilter:',minfilter,' maxfilter:',maxfilter
		
		z_min = p10p90W[minfilter,0]/L_Lya0 - 1
		z_max = p10p90W[maxfilter,1]/L_Lya0 - 1

		medfilter = (minfilter+maxfilter)/2.
		zbin = p10p90W[medfilter,1]/L_Lya0 - p10p90W[medfilter,0]/L_Lya0

		print 'sorting redshifts and distances'

		zsort = np.sort(GalArr['redshift'])
		dsort = np.sort(distance)
		
		print 'interpolating...'

		rmin = np.interp(z_min,zsort,dsort)
		rmax = np.interp(z_max,zsort,dsort)

		rdata = distance[nng]
		zdata = np.interp(rdata,dsort,zsort)

#		zbin = 0.02

		NDistanceBins = (z_max - z_min)/zbin + 1.0

		FileSnapshots = '../data/Millennium_outputs.txt'
		datasnap = np.loadtxt(FileSnapshots,dtype=np.dtype('f4','f4'),usecols=(0,1))
		nsnap = len(datasnap[:,0])

		print 'building selection function...'
		nghist,zarr = np.histogram(zdata,bins=np.sort(datasnap[:,1]))
		lenr = len(zarr)

		zarr2 = zarr[:lenr-1]
		
		TotalVol = (1./6)*np.pi*(rmax**3 - rmin**3)	
		AreaRad = 4*np.pi/8.

		selec_function_all = np.zeros(len(zarr2))
		deg_factor = (180.0**2)/(2*np.pi)
		for isf in np.arange(len(zarr2)):
			selec_function_all[isf] = nghist[isf]/deg_factor/(zarr[isf+1] - zarr[isf])

#		selec_function_all = nghist/((180.0**2)/(2*np.pi))/zbin

		pl.plot(zarr2,np.log10(selec_function_all),'-',linewidth=1,color=pl.cm.rainbow(iM*30))


#		if iM == 0:
#			nghist_all,zarr_aux = np.histogram(GalArr['redshift'],bins=NDistanceBins,range = [z_min-zbin/2.,z_max+zbin/2.])
#			pl.plot(zarr2,np.log10(nghist_all/TotalVol/zbin),'-',linewidth=2,color='green',label='All galaxies')


		pl.ylabel(r'$\log ({\rm d}N/{\rm d}z \  [\rm deg^{-2}])$',fontsize=20)
		pl.xlabel(r'$z$',fontsize=20)
		pl.xlim([1.5,5.5])
		pl.ylim([-2,2.5])
		pl.legend(loc=1,fontsize=10)

		maxsf = np.max(selec_function_all)
		maxft = np.max(FiltArr['t'])
		
		col1 = 'darkgray'
		j = 0
		for i in np.arange(minfilter,maxfilter):
			col2 = pl.cm.rainbow(i*10)
			if j == 0:
				col = col1
				j = 1
			else:
				col = col2
				j = 0

			zw = FiltArr['w'][i]/L_Lya0 - 1
			pl.plot(zw,3*((FiltArr['t'][i]/maxft)) - 2 ,color=col,linewidth=1)


#		for iz in np.arange(nsnap):
#			pl.plot([data[iz,1],data[iz,1]],[-9,-3],'--',color='black')
		
		pl.savefig('../plots/selection_function_laes_dz.png',bbox_inches='tight',dpi=400)	
		print 'done.'
		

# Checking that the sel function is well computed

	ngals_sf = integral.simps(selec_function_all,zarr2) * 5000.
	print 'integral of selec function vs. total n. of galaxies',ngals_sf,len(nng)	

	
	

	FileFilters = OutputDir + 'jpas_laes_'
	k = 0
	for ip in np.arange(nFil):
		gals_ip = np.where((id_filter == ip) & 
				  (mag_app <= jpaslims[ip][app_type]) &
				  (mag_app > 0))
		ng = len(gals_ip[0])
		X = np.zeros((ng,4))
		X[:,0:3] = GalArr['pos'][gals_ip]
		X[:,3] = GalArr['redshift'][gals_ip]
		FileOut = FileFilters + str(ip) 
		ascii.write(X,FileOut)


#	np.savetxt(FileAll,X,fmt= 'Iteration %10.5f')
#	f.close()	


if CheckMock:

	bright = np.where(mag_app < 23 - 5*np.log10(0.72))

	nx = 50
	ny = 50
	nz = 10.0

	mean_sfr = np.zeros([nx,ny])
	mean_z   = np.zeros([nx,ny])
	mean_zcold   = np.zeros([nx,ny])

	rad2deg = (180.0/np.pi)
	phi = np.arctan2(GalArr['pos'][bright,1],GalArr['pos'][bright,0]) * rad2deg
	theta = np.arccos(GalArr['pos'][bright,2]/distance[bright]) * rad2deg

	pmin = np.min(phi)
	pmax = np.max(phi)
	tmin = np.min(theta)
	tmax = np.max(theta)

	pbin = (pmax - pmin)/(nx-1)
	tbin = (tmax - tmin)/(ny-1)

#	xmin = 0
#	xmax = np.max(GalArr['pos'][izr[0][bright],0])
	
#	ymin = 0
#	ymax = np.max(GalArr['pos'][izr[0][bright],1])

	zmin = 3.0
	zmax = 4.5

#	xbin = (xmax - xmin)/(nx-1)
#	ybin = (ymax - ymin)/(ny-1)
	zbin = (zmax - zmin)/(nz-1)


	xall     = GalArr['pos'][bright,0]
	yall     = GalArr['pos'][bright,1]
	sfrall   = GalArr['sfr'][bright]
	zall     = GalArr['redshift'][bright]
	zcoldall = GalArr['Zcold'][bright]

	for iz in np.arange(nz):
		for ix in np.arange(nx):
			for iy in np.arange(ny):
				_x = ix*pbin + pmin
				_y = iy*tbin + tmin
				_z = iz*zbin + zmin
				sel = np.where( (phi > _x-pbin/2.) &
							(phi < _x+pbin/2.) &
							(theta > _y-tbin/2.) &
							(theta < _y+tbin/2.) & 
							(zall > _z-zbin/2.) &
							(zall < _z+zbin/2.))
	
				mean_sfr[ix][iy] = np.mean(np.log10(sfrall[sel]))
				mean_z[ix][iy] = np.mean(zall[sel])
				mean_zcold[ix][iy] = np.mean(np.log10(zcoldall[sel]))
	
#		pl.matshow(mean_sfr)
		ssel = np.where((zall > _z-zbin/2.) &
						(zall < _z+zbin/2.))

		pl.close()
		a1 =  pl.figure(1)
		pl.imshow(mean_sfr,extent=[pmin,pmax,tmin,tmax],interpolation='nearest')
		pl.xticks(np.arange(pmin,pmax,10))
		pl.yticks(np.arange(tmin,tmax,10))
		pl.xlabel('deg')
		pl.ylabel('deg')		
		pl.autoscale(False)
		pl.plot(theta[ssel],phi[ssel],'.',color='gray',alpha=0.5)
		cbar = pl.colorbar()
		cbar.set_label(r'$\log({\rm SFR} [M_{\odot}/yr)]$',rotation=270)
		pl.text(10,10,r'$'+trunc(_z-zbin/2.,1)+'<z<'+trunc(_z+zbin/2.,1)+'$',bbox=dict(facecolor='0.85', alpha=0.85),fontsize=18)
		pl.savefig('sfr_iz'+trunc(iz,0)+'.png',bbox_inches='tight',dpi=350)
		pl.close(a1)
#		pl.matshow(mean_zcold)

		a2 = pl.figure(2)
		pl.imshow(mean_zcold,extent=[pmin,pmax,tmin,tmax],interpolation='nearest')
		pl.xticks(np.arange(pmin,pmax,10))
		pl.yticks(np.arange(tmin,tmax,10))
		pl.xlabel('deg')
		pl.ylabel('deg')		
		pl.autoscale(False)
		pl.plot(theta[ssel],phi[ssel],'.',color='gray',alpha=0.5)
		cbar = pl.colorbar()
		cbar.set_label(r'$\log(Z_{\rm cold})$',rotation=270)
		pl.text(10,10,r'$'+trunc(_z-zbin/2.,1)+'<z<'+trunc(_z+zbin/2.,1)+'$',bbox=dict(facecolor='0.85', alpha=0.85),fontsize=18)
		pl.savefig('zcold_iz'+trunc(iz,0)+'.png',bbox_inches='tight',dpi=350)
		pl.close(a2)
			

if Make_FigCompareMags:
	fig = pl.figure()
	pl.subplot(1,1,1,xlim=[20,35],ylim=[20,35],aspect=1.0)
	pl.plot(mag2_app,mag_app,',')
	pl.plot([20,35],[20,35],'--',color='gray')
	pl.xlabel(r'$m_{AB}(f_{\lambda})$',fontsize=15)
	pl.ylabel(r'$m_{AB}(\langle f_{\lambda} \rangle)$',fontsize=15)	
#	pl.axis([20,35,20,35]0i)
#	pl.axis('equal')
#	pl.setp(fig,aspect='equal')
	pl.savefig('mag_comp.png',dpi=400)

