# Compute the LF of a given line at a given redshift in the lightcone. It requires reading a lot of data  

import read_lines_lc as read
import multiprocessing as mp
import numpy as np 

filename = '/home/CEFCA/aaorsi/work/jpas_mock/out/lines_z.'
nfiles = 20
from functools import partial

# Read all files:
def read_lines(zrange,ip):  
  fname = filename + str(ip)
  data = read.read_chunks(ip)
  print 'file '+fname+' read.'
  sel = np.where((data['z'] >= zrange[0]) & (data['z'] <= zrange[1]))
  print len(sel[0])
  print data['lines'][sel,:]
  import pdb ; pdb.set_trace()
  return {'lines':data['lines'][sel,:],'z':data['z'][sel]}

zrange = [1.2,2.7]
func = partial(read_lines, zrange)
nproc = 1
pool = mp.Pool(processes=nproc)
#data = pool.map(func, range(nfiles))
read_lines(zrange,0)



