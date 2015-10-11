"Script to return Ly-alpha luminosities with scaled f_esc"

import numpy as np

from readmock import *
from load_laedata import *

def get_lya(zrange):
	GalArr,Lya0,NGals,distance = readmock()
	h = 0.73
	facl = -2*np.log10(h)
	facn = +3*np.log10(h)
	hayes = laedata('hayes_z2')
	blanc = laedata('blanc_z19_38')
	ciardullo = laedata('ciardullo')
	gronwall = laedata('gronwall_z3_cum')
	ouchi	 = laedata('ouchi_z3_cum')
	izr = np.where((GalArr['redshift'] < zrange[1]) &
	  	  (GalArr['redshift'] > zrange[0]))
#	d1 = np.max(GalArr['distance'][izr[0]])
#	d0 = np.min(GalArr['distance'][izr[0]])
	d1 = np.max(distance)
	d0 = np.min(distance)
	print 'distances between ',d0,d1
	Volume_z = 4./3*np.pi*(d1**3 - d0**3)/8.0
	Lyaz = Lya0[izr]

	GalArrzr = GalArr[izr]

#	return Lyaz,izr,GalArr[izr],Volume_z
	return Lyaz,izr,GalArrzr,Volume_z


