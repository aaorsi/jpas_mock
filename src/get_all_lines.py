#!/usr/bin/python
"Script to calculate magnitudes and properties for all emission lines in J-PAS lightcones."
" Based in get_oii.py, which is based in get_mags.py"

import pdb
import numpy as np
import pylab as pl
import weighted as wp
import scipy.integrate as integral
import read_filters as rf
from readmock_chunk import *
import get_jpaslims as jp
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM 
import multiprocessing as mp
from load_oiidata import *
from read_photoion import *
from get_lf import *
import struct
import sys

def find_nearest_vector(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return str(f)[:slen]

mpc2cm_24 = 3.0856778 # mpc-to-cm /1e24

def main(nproc,outfile):
  """
  Creates files containing all emission line luminosities and
  magnitudes for lightcone galaxies.
  Calculation is split in a given number of processes, where
  each computes the emission lines of a chunk of data independently
  and outputs the results to a file.

  Arguments:

  nproc = number of processes
  outfile = Output file name. (final name will be followed by a number)

  """
# nproc = argv[0]
# outfile = argv[1]

# if (len(argv) < 2):
#   print 'you need to provide 2 arguments: nproc outfile.'
#   print 'Exiting..'
#   return 0

  nproc = int(nproc)

  IncludeMags = True

  print 'nproc:',nproc,' outfile:',outfile
  sys.stdout.flush()
  #L_line  = 1215.67  # Angstroms
  L_line  = 3727.0  # Angstroms
  #Lyaz,izr,GalArr,Vol = get_lya(zrange)

  FiltArr,wrange = rf.read_filters()
  jpaslims = jp.get_jpaslims()

  app_type = 1  # 0 means per 3 arcsec diameter, 1 means per square arcsec
  cosmo = FlatLambdaCDM(H0=73, Om0=0.25)
  zrArr = comparat_zranges()

  zminlist = map(float,zrArr['zmin'])
  zmaxlist = map(float,zrArr['zmax'])
  minz = np.min(zminlist)
  maxz = np.max(zmaxlist)
  #maxz = np.max(zrArr['zmax'])

  print minz,maxz
  minsfr = 0.1

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

  print 'loading lightcone data'
  props = ['redshift','pos','Zcold','sfr','vel','stellarmass']
  print 'now loading photo-ionisation grid'
  lineinfo,linesarr = read_photoion()
  nline = lineinfo['nlines']
  usezspace = 1
  
  # magtype stores both the value of the magnitude (float) and the
  # filter id (int) where that magnitude was computed.

  magtype = np.dtype('f4, i4')
  def read_lightcone_chunk(ip, nproc):
    GalArr,lLyc,ngg,DistNz,zdist = readmock_chunk(ip,nproc,props_array = props,
                                   zspace=usezspace)
    qpar = qZrelation(GalArr['Zcold'])  # assuming default pars.

    # This stores all line luminosities for all galaxies
    # LinesLumArr[i,j] = Luminosity of line j for galaxy i.
    LinesLumArr = np.zeros((ngg,nline))
    redshiftArr = np.zeros(ngg)
    SfrArr = np.zeros(ngg)
    MStellarArr = np.zeros(ngg)
    print 'ip '+str(ip)+', computing lines for  ngals=',ngg,' galaxies...'
    sys.stdout.flush()
    for ig in range(ngg):
      LinesLumArr[ig,:] = integ_line(lineinfo,linesarr,qpar[ig],GalArr['Zcold'][ig],lLyc[ig],all_lines=True) 
      redshiftArr[ig] = GalArr['redshift'][ig]
      SfrArr[ig] = GalArr['sfr'][ig]
      MStellarArr[ig] = GalArr['stellarmass'][ig]
#     LinesLumArr[ig,il] = iline
      
    print 'redshiftArr', redshiftArr[0:10]
    print 'Lines[0,:]' , LinesLumArr[0,:]
    print 'Lines[10,:]' , LinesLumArr[10,:]
    print 'Lines[100,:]' , LinesLumArr[100,:]
    print 'Proc ip:'+str(ip)+' done with that.'

    # This stores the AB observed-frame apparent magnitude for all luminosities of all galaxies, indicating also to which filter
    # the stored magnitudes corresponds to.
    # MagsLumArr[i,j] = (AB magnitude, filter id) of line j of galaxy i.

    MagsLumArr = np.zeros((ngg,nline),dtype=magtype)
    realdist = np.sqrt(GalArr['pos'][:,0]**2 + GalArr['pos'][:,1]**2+GalArr['pos'][:,2]**2)
    log_dist = np.log10(realdist) + np.log10(mpc2cm_24)+24.0 + np.log10(1+GalArr['redshift'])

    if IncludeMags is True:
      for il in range(nline):
        L_line = lineinfo['lambda0'][il]
        wgal = L_line*(1 + GalArr['redshift'])
        log_Fline = LinesLumArr[:,il] - (np.log10(4*np.pi)+2*(log_dist))
        FLine = 10**log_Fline

        for i in range(ngg):
          idFilter,Val = find_nearest_vector(CentralW,wgal[i])
          nz = np.where(FiltArr['w'][idFilter] != 0)
          lamb = FiltArr['w'][idFilter][nz]
          tlamb = FiltArr['t'][idFilter][nz]
          tlam0 = np.interp(wgal[i],lamb,tlamb)
    # mag_app[i] = -2.5*(log_Fline[i]+ np.log10((wgal[i]*tlam0)/int_ab[idFilter]) - 18-np.log10(3.0)) - 48.60
          mag_app = -2.5*(log_Fline[i]+ np.log10((wgal[i]*tlam0)/int_ab[idFilter]) - 18-np.log10(3.0)) - 48.60

          MagsLumArr[i,il] = (mag_app, idFilter)
  
    filename = outfile + '.' + str(ip)
    nf = open(filename,"wb")
    print 'writing file '+filename
    nf.write(struct.pack('l',ngg))
    nf.write(struct.pack('i',nline))
    LinesLumArr.tofile(nf)
    if IncludeMags is True:
      MagsLumArr.tofile(nf)
    redshiftArr.tofile(nf)
    SfrArr.tofile(nf)
    MStellarArr.tofile(nf)
    nf.close()
    print 'file '+filename+' written succesfully.'
    
    return 1

  if nproc == 1:
    processes = read_lightcone_chunk(0,1)
  else:  
    print 'start processes'
    processes = [mp.Process(target=read_lightcone_chunk, args=(ip,nproc)) for ip in range(nproc)]
    for p in processes:
      p.start()

    for p in processes:
      p.join()


  return 1

if __name__ == "__main__":
  main(sys.argv[1],sys.argv[2])
