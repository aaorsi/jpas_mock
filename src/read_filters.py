"Script to read JPAS filters"

import numpy as np

def read_filters():

	NBinMax = 200
	fctype = np.dtype([('w',np.float32),('t',np.float32)])

	DataDir = '/home/CEFCA/aaorsi/work/jpas_mock/data/'
	wrange_ = np.arange(3900,9100,100)
	nn = len(wrange_)
	wrange = np.append([3518,3785],wrange_)
	nbfile = np.append([535,165],np.zeros(nn)+145)

	FilterData = np.zeros([nn+2,NBinMax],dtype=fctype)

	for i in np.arange(nn+2):
		File = DataDir + 'JPAS'+ str(wrange[i]) + '_' + \
		str(int(nbfile[i]))+'.res'
		print File
		data = np.loadtxt(File,dtype=fctype)
		ndata = len(data['w'])
		FilterData['w'][i][0:ndata] = data['w']
		FilterData['t'][i][0:ndata] = data['t']

	return FilterData, wrange
#	Data read. todo: compute filter magnitude for each LAE from jpas mock catalogue.



