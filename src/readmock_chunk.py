import numpy as np
import pylab as pl

def readmock_chunk(ip, nproc, props_array = ['All'], zmin = 0, zmax = 1000, sfrmin = 0.0, 
mstellarmin = 0.0, zspace = None):
	"""
	This function reads the MXXL lightcone of J-PAS and returns a catalogue of galaxies
	with limits given by the keywords. If no arguments are defined, then
	the catalogue will return ALL galaxies and ALL properties

	Parallelisation with multiprocessing should be easy (?), but it's not implemented yet.

	props_array = list of strings containing the name of the properties in the output. If not specified, returns all properties.
	zmin, zmax: redshift range to extract the data from the lightcone
	sfrmin: minimum SFR
	mstellarmin: minimum stellar mass
	zspace: compute redshift space. This includes s coordinate, i.e. distance in redshift space and z_s, which includes the redshift distortion.

	"""

	datadir = '/home/CEFCA/aaorsi/work/jpas_mock/data/'
	filename = 'gals_mxxl_v0_5.sav'
	datafile = datadir + filename

	fd = open(datafile,'rb')

	ngals = np.fromfile(fd,dtype='i4',count=1)
	print 'number of galaxies in catalogue:', ngals

	ngals_chunk = (ngals / nproc)
	print 'number of galaxies in chunk '+str(ip)+'='+str(ngals_chunk)
	
#	dtype used to read the data:
	gal = np.dtype([('redshift',np.float32),('stellarmass',np.float32),('pos',np.float32,3),('vel',np.float32,3),\
	('ObsMag',np.float32,5),('Mvir',np.float32),('Type',np.int32),('sfr',np.float32),('Mz',np.float32),('Mcold',np.float32)])

# dtype used to store the data:
	len_proparr = len(props_array)
	if len_proparr == 1:
		print 'All data stored'
		galstore = gal
	else:
		for i in range(len_proparr):
			iprop = props_array[i]
			varlen = 1
			if iprop == 'pos' or iprop == 'vel':
				varlen = 3
			if iprop == 'ObsMag':
				varlen = 5

			if (i == 0):
				galstore = [(iprop,np.float32,varlen)]
			else:
				galstore.append((iprop,np.float32,varlen))
		
	GalArr = np.empty(ngals_chunk,dtype=galstore)

	i = j = k = nchunks = 0
	while i < ip:
		data = np.fromfile(fd,dtype=gal,count=ngals_chunk)
		i += 1

	data = np.fromfile(fd,dtype=gal,count=ngals_chunk)
	# All conditions are applied below:
	karr = np.where((data['redshift'] >= zmin) &
	(data['redshift'] < zmax) & 
	(data['sfr'] >= sfrmin) &
	(data['stellarmass'] >= mstellarmin))
	nk = len(karr[0])
	if len_proparr == 1:
		GalArr[k:k+nk] = data[karr]
	else:
		for ip in range(len_proparr):
			iprop = props_array[ip]
			GalArr[iprop][0:nk] = data[iprop][karr]

	del data
	del karr
	
	fd.close
	gsel = np.arange(nk)

	distance = np.sqrt(GalArr['pos'][gsel,0]**2 + GalArr['pos'][gsel,1]**2+GalArr['pos'][gsel,2]**2)
	Nlyc_ = np.log10(1.35) + np.log10(GalArr['sfr'][gsel]) + 53.0

#	nozero = np.where(GalArr['sfr'][gsel] > 0)
	ngg = nk
	print 'after selection, ngals = ',ngg
	if zspace:
		vr = ((GalArr['pos'][gsel,0]*GalArr['vel'][gsel,0])/distance +
		(GalArr['pos'][gsel,1]*GalArr['vel'][gsel,1])/distance +
		(GalArr['pos'][gsel,2]*GalArr['vel'][gsel,2])/distance)
		a_exp = 1 / (1 + GalArr['redshift'][gsel])
		Om = 0.25
		Ol = 0.75
		E_z = np.sqrt(Om * (1 + GalArr['redshift'][gsel])**3 + Ol)

		s = distance + vr / (100*a_exp*E_z) 


	if zspace:
#	return GalArr[nozero], Lya0[nozero], ngg, distance[nozero]
		return GalArr[0:ngg], Nlyc_, ngg, distance,s
	else:
		return GalArr[0:ngg], Nlyc_, ngg, distance
		
if __name__ == "__readmock_chunk__":
	readmock()

