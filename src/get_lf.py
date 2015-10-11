import numpy as np

lftype = np.dtype([('Lum',np.float32),('dn_dL',np.float32
),('err',np.float32)])

def get_lf( data, Vol=1.0,binsize=1.0,minsample=10.0,maxsample=100):
	
#	maxsample = np.max(data)
	nbins = np.ceil((maxsample-minsample)/binsize) +1.0
	larray = np.linspace(minsample,maxsample,num=nbins)
	lf = np.zeros(nbins)	

	j = 0
	ndata = len(data)
#	print 'nbins = ',nbins
#	print 'ndata = ',ndata 

#	data = np.sort(data)
#	data = list(data).sort)
	data.sort()


	while data[j] < larray[0] :
		j += 1

	esc = 0
	for i in np.arange(nbins):
		while (data[j] < larray[i]+binsize/2. and data[j] >=
		larray[i]-binsize/2. and j < ndata-1):
			lf[i] += 1
			j +=1
			if j == ndata-1:
				esc = 1
				break
#			print j,i,lf[i],data[j],larray[i]

	
	lf = np.log10(lf / Vol / binsize)

	lfresult = np.zeros(nbins,dtype=lftype)
	lfresult['Lum'] = larray
	lfresult['dn_dL'] = lf
	
	return lfresult

if __name__ == "__get_lf__":
	get_lf()	


