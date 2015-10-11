import numpy as np
import pylab as pl
import sys

def readmock(props_array = ['All'], zmin = 0, zmax = 1000, sfrmin = 1e-10, mstellarmin = 1e-10, zspace = None):
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
	filename = 'gals_mxxl_v0_4.sav'
	datafile = datadir + filename

	fd = open(datafile,'rb')
	SizeChunk = 1.0e10

	ngals = np.fromfile(fd,dtype='i4',count=1)
	print 'number of galaxies in catalogue:', ngals
	sys.stdout.flush()

#	gal = np.dtype([('redshift',np.float32),('stellarmass',np.float32),('pos',np.float32,3),('vel',np.float32,3),('imag',np.float32),\
#	('ObsMag',np.float32,5),('distance',np.float32),('sfr',np.float32)])

#	dtype used to read the data:
	gal = np.dtype([('redshift',np.float32),('stellarmass',np.float32),('pos',np.float32,3),('vel',np.float32,3),\
	('ObsMag',np.float32,5),('Mvir',np.float32),('Type',np.int32),('sfr',np.float32),('Zcold',np.float32)])

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
		
	
	size_gal = gal.itemsize
	nGals_step = int(SizeChunk / size_gal)
	print 'number of galaxies to be kept in memory:',nGals_step
	GalArr = np.empty(2e8,dtype=galstore)

	i = j = k = nchunks = 0
	while i <= ngals:
		data = np.fromfile(fd,dtype=gal,count=nGals_step)
		nchunks += 1
		print 'istep ',nchunks
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
				GalArr[iprop][k:k+nk] = data[iprop][karr]

		k += nk
		print 'scanned...ok'
		i += nGals_step
		print k
		del data
		del karr
	
	fd.close

	gsel = np.arange(k)

	distance = np.sqrt(GalArr['pos'][gsel,0]**2 + GalArr['pos'][gsel,1]**2+GalArr['pos'][gsel,2]**2)
	Nlyc_ = np.log10(1.35) + np.log10(GalArr['sfr'][gsel]) + 53.0

	nozero = np.where(GalArr['sfr'][gsel] > 0)
	ngg = len(nozero[0])
	if zspace:
		vr = ((GalArr['pos'][gsel,0]*GalArr['vel'][gsel,0])/distance +
		(GalArr['pos'][gsel,1]*GalArr['vel'][gsel,1])/distance +
		(GalArr['pos'][gsel,2]*GalArr['vel'][gsel,2])/distance)
		a_exp = 1 / (1 + GalArr['redshift'][gsel])

		Om = 0.25
		Ol = 0.75
		E_z = np.sqrt(Om * (1 + GalArr['redshift'][gsel])**3 + Ol)

		s = distance + vr / (100*a_exp*E_z) 
		ss = s[nozero]

	GalArr = GalArr[nozero]
	distance = distance[nozero]
	Nlyc_ = Nlyc_[nozero]

	if zspace:
		s_nz	 = s[nozero]
#	return GalArr[nozero], Lya0[nozero], ngg, distance[nozero]
		return GalArr, Nlyc_, ngg, distance,ss
	else:
		return GalArr, Nlyc_, ngg, distance
		
if __name__ == "__readmock__":
	readmock()

