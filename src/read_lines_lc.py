
#!/usr/bin/python

import pdb
import numpy as np
import multiprocessing as mp
import pylab as pl
from get_lf import *
import read_filters as rf
import get_jpaslims as jp
import weighted as wp
import scipy.integrate as integral
import sys
import struct
from read_photoion import *

FiltArr,wrange = rf.read_filters()
print 'filters read'
jpaslims = jp.get_jpaslims()
print 'jpas lims set'

datafile = '/home/CEFCA/aaorsi/work/jpas_mock/out/'
#name = 'emlines_MXXL'
#name = 'test'
name = 'lines_mags_sfr_mstellar_z'


nFil = len(FiltArr['w'])
CentralW = np.zeros(nFil)
int_ab = np.zeros(nFil)
p10p90W = np.zeros([nFil,2])

  
lineinfo,linesarr = read_photoion()
nline = lineinfo['nlines']

for i in range(nFil):
  CentralW[i] = wp.quantile(FiltArr['w'][i],FiltArr['t'][i],0.5)
  nz = np.where(FiltArr['w'][i] != 0)
  int_ab[i] = (integral.simps(FiltArr['t'][i][nz]/FiltArr['w'][i][nz],FiltArr['w'][i][nz]))
  p10p90W[i,0] = wp.quantile(FiltArr['w'][i][nz],FiltArr['t'][i][nz],0.1)
  p10p90W[i,1] = wp.quantile(FiltArr['w'][i][nz],FiltArr['t'][i][nz],0.9)
  print 'int_ab[i]',int_ab[i]


def read_chunks(ip,Mags=True):

  filename = datafile + name +'.' + str(ip)
#  print 'reading '+filename
  nf    = open(filename,"rb")
  ngg   = np.fromfile(nf,dtype=np.int64,count=1)
  nline = np.fromfile(nf,dtype=np.int32,count=1)
  nline = nline[0]
  ngg = ngg[0]
  magtype = np.dtype('f4, i4')
#  print ngg, nline
  LinesLumArr = np.zeros([ngg,nline])
#   LinesLumArr = nf.read(LinesLumArr.nbytes)
  LinesLumArr = np.fromfile(nf,dtype='('+str(ngg)+','+str(nline)+')f8',count=1)
#  print 'lines data read.'

  MagsLumArr = np.zeros((ngg,nline),dtype=magtype)
  if Mags is True:
    for _i in range(ngg):
      MagsLumArr[_i,:] = np.fromfile(nf,dtype=magtype,count=nline)
#    nf.close()
#    return LinesLumArr, MagsLumArr_
#  else:
  zArr = np.fromfile(nf,dtype=np.double,count=ngg)
  SfrArr = np.fromfile(nf,dtype=np.double,count=ngg)
  MStellarArr = np.fromfile(nf,dtype=np.double,count=ngg)



#    print 'zdata read.'
  nf.close()
  return {'lines':LinesLumArr, 'z':zArr,'Mstellar':MStellarArr,'SFR':SfrArr,'Mag':MagsLumArr}
#   MagsLumArr = np.memmap(nf,dtype=magtype, shape=(ngg,nline))
  

def get_lfmags(LinesLumArr,MagsLumArr_):

  allmag_lf = np.zeros(nFil,dtype=lftype*int(nlbins))
  MagsLumArr = np.reshape(MagsLumArr_,(ngg,nline))

#   nf.read(struct.unpack('i',nline))
#   LinesLumArr = np.zeros((ngg,nline))
#   MagsLumArr = np.zeros((ngg,nline),dtype=magtype)
#   LinesLumArr.fromfile(nf)
#   MagsLumArr.fromfile(nf)

  print 'file '+filename+' read succesfully.'
  for j in range(25,nFil):
    print 'computing the LF for filter '+np.str(j)
    mag_sel = np.where(MagsLumArr[:,iline]['f1'] == j)
    if len(mag_sel[0]) == 0:
      continue

    mags = MagsLumArr[mag_sel,iline]['f0']
    allmag_lf[j] = get_lf(mags[0],
    Vol = 1.0, binsize=lfbin,minsample = lfmin,maxsample=lfmax)
    print 'LF computed for filter ',j

  return allmag_lf  



def main(iline):
  """
  Read light-cone emission line files. Used to compute LFs
  """

  lfbin = 0.25
  lfmax = 26.0
  lfmin = 15.0
  nlbins = (lfmax - lfmin)/(lfbin) + 1.0


# output = mp.Queue()

#     output.put(allmag_lf)
    
  nproc = 1
  pool = mp.Pool(processes=nproc)
  lfarray  = pool.map(read_chunks, range(nproc))

# mp.Process(target=read_chunks, args=(ip,)) for ip in range(nproc)]
# mags = read_chunks(0)
# pdb.set_trace()
# for p in processes:
#   p.start()
# for p in processes:
#   p.join()
# lfarray = [output.get() for p in processes]

  pdb.set_trace()

if __name__ == "__main__":
  main(sys.argv[1])


