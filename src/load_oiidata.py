"obs. data functions"
import numpy as np

lftype = np.dtype([('Lum',np.float32),('dn_dL',np.float32
),('err',np.float32)])

def qZrelation(Zgas, q0 = 2.8e7, g0 = -1.3):
	return q0 * (Zgas/0.012)**g0


def comparat_zranges():
	"""
	This function defines the redshift ranges
	of Comparat+14 LFs of [OII] 3727+3729 and
	a mean value to point at them
	"""
	zmin = ['1.340', '1.090', '0.880', '0.695', '0.100','0.240','0.500' ] 
	zmax = ['1.65', '1.34', '1.09', '0.88', '0.24','0.4','0.695']
	zmean= [1.4 , 1.2 , 1.0 , 0.7  , 0.2  ,0.3, 0.6 ]

	zrange_comparat = {'zmin':zmin,'zmax':zmax,'zmean':zmean}

	return zrange_comparat


def oiidata( name, z=None, quiet=False ):
	data ='False'
	if name == 'comparat':
		if (z == None):
			print 'comparat data set needs z as an argument. Exiting'
			exit()
		if (not quiet):
			print 'loading comparat LF for z=',z
		
		zcom = comparat_zranges()
		_iz = np.where(np.array(zcom['zmean']) == np.array(z))
		print zcom['zmean'],z,_iz
		iiz = int(_iz[0])
#		_iz = _iz[0]
	
		print _iz
		DataDir = '/home/CEFCA/aaorsi/work/jpas_mock/data/comparat/'
#		zmin = '%.'+str(zcom['mind'][iiz]-2)+'f' % zcom['zmin'][iiz]
#		zmin = '%.'+str(zcom['maxd'][iiz]-2)+'f' % zcom['zmax'][iiz]
		zmin = zcom['zmin'][iiz]
		zmax = zcom['zmax'][iiz]
		File = DataDir+'comparat2014-LF-'+zmin+'-z-'+zmax+'.dat'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,2]
		data['dn_dL'] = dread[:,3]
		data['err'] = dread[:,4]

	return data
	

