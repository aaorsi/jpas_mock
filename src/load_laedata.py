"obs. data functions"
import numpy as np

lftype = np.dtype([('Lum',np.float32),('dn_dL',np.float32),('err',np.float32)])

def laedata( name ):
	DataDir = '/home/aaorsi/work/laes_obsdata/LF/'
	data ='False'
	if name == 'hayes_z2':
		File = DataDir+'hayes_z2_data.txt'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
		ErrFile = DataDir + 'hayes_z2_hierr.txt'
		derr = np.genfromtxt(ErrFile,comments='#')
		data['err'] = 10**derr[:,1] - 10**data['dn_dL']
	
	if name == 'blanc_z19_38':
		File = DataDir+'blanc_z19_38.txt'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
		ErrFile = DataDir + 'blanc_z19_38_hierr.txt'
		derr = np.genfromtxt(ErrFile,comments='#')
		data['err'] = 10**derr[:,1] - 10**data['dn_dL']

	if name == 'blanc_z19_28':
		File = DataDir+'blanc_z19_28.txt'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
 		data['err'] = zeros(nbins)
		
	if name == 'ciardullo':
		File = DataDir+'ciardullo_z2.txt'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
 		data['err'] = 10**(dread[:,1]+dread[:,3])-10**(dread[:,1])

	if name == 'gronwall_z3':
		File = DataDir + 'lf_musyc_3.1'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
 		data['err'] = dread[:,2]

	if name == 'ouchi_z3':
		File = DataDir + 'lf_ouchi_3.1'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
 		data['err'] = dread[:,2]

	if name =='gronwall_z3_cum':	# Cumulative LF
		File = DataDir + 'cum_gronwall07_z3p1'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
 		data['err'] = dread[:,2]

	if name =='ouchi_z3_cum':	# Cumulative LF
		File = DataDir + 'clf_ouchi_3.1'
		dread = np.genfromtxt(File,comments='#')
		nbins = len(dread)
		data = np.empty(nbins,dtype=lftype)
		data['Lum'] = dread[:,0]
		data['dn_dL'] = dread[:,1]
 		data['err'] = dread[:,2]


	return data

	
	

			
			
				
	



	






