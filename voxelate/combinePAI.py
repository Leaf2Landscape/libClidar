import numpy as np
import argparse
from sys import exit
import csv
import h5py



###############################
# Combines output of voxelTLS #
# from multiple scans         #
#                             #
# svenhancock@gmail.com       #
# November 2018               #
###############################


##############################################
# read the command line
##############################################

def getCmdArgs():
    '''
    Get commdnaline arguments
    '''
    p = argparse.ArgumentParser(description=("Combines voxel PAI estimates from multiple scans"))
    p.add_argument("--input", dest ="inNamen", default=" ",type=str, help=("Input of combined scans\nNo default"))
    p.add_argument("--inList", dest ="listNamen",default=" ",type=str, help=("Input list of single scans\nNo default"))
    p.add_argument("--output", dest ="outNamen",default="meanVox.vox", type=str,help=("Outpt filename\nDefault = meanVox.vox"))
    p.add_argument("--minGap", dest ="minGap", type=float, default=0.0, help=("Minimum trustable gap fraction to a voxel\nDefault = 0"))
    p.add_argument("--hdf", dest ="writeHDF", action='store_true',default=False,help=("Write output in HDF5\nDefault is ASCII"))
    cmdargs = p.parse_args()
    return cmdargs


##############################################

class voxelData(object):
  '''
  A class to contain voxel
  gap fractions and LAIs
  '''

  ##############################################

  def __init__(self,inNamen,listNamen):
    '''Read the data'''
    # do we have a single file
    if(inNamen!=" "):
      self.readAllScans(inNamen)
    elif(listNamen!=" "):   # multiple files
      self.readIndScans(listNamen)
    else:   # no input. Fail program
      print("No input specified")
      exit(1)
    # determine resolutions and numbers
    if(self.hdf==False):
      self.findDims()


  ##############################################

  def findDims(self):
    '''Find dimensions of voxel file'''
    # scene bounds
    self.minX=np.min(self.x)
    self.maxX=np.max(self.x)
    self.minY=np.min(self.y)
    self.maxY=np.max(self.y)
    self.minZ=np.min(self.z)
    self.maxZ=np.max(self.z)
    # number per
    self.nX=len(set(self.x))
    self.nY=len(set(self.y))
    self.nZ=len(set(self.z))
    # resolutions
    self.xRes=(self.maxX-self.minX)/self.nX
    self.yRes=(self.maxY-self.minY)/self.nY
    self.zRes=(self.maxZ-self.minZ)/self.nZ


  ##############################################

  def readIndScans(self,listNamen):
    '''read each scan from separate files'''
    # read list and count scans
    files=np.genfromtxt(listNamen,unpack=True,usecols=(0,),dtype=str,comments='#',delimiter=" ")
    nScans=len(files)
    i=0
    # loop over files
    for inNamen in files:
      # is this an ASCII or a HDF5 file?
      ending=inNamen.strip().split('.')[-1]
      if((ending=="h5")|(ending=="hdf")):  # HDF5
        self.hdf=True
        f=h5py.File(inNamen.strip(),'r')
        self.nX=list(f['NX'])[0]
        self.nY=list(f['NY'])[0]
        self.nZ=list(f['NZ'])[0]
        self.nVox=self.nX*self.nY*self.nZ
        # if first, allocate space
        if(i==0):
          self.x=np.full(self.nVox,-1,dtype=float)
          self.y=np.full(self.nVox,-1,dtype=float)
          self.z=np.full(self.nVox,-1,dtype=float)
          self.bounds=np.array(f['BOUNDS'])
          self.minX=self.bounds[0]
          self.minY=self.bounds[1]
          self.minZ=self.bounds[2]
          self.maxX=self.bounds[3]
          self.maxY=self.bounds[4]
          self.maxZ=self.bounds[5]
          self.res=np.array(f['RES'])
          self.xRes=self.res[0]
          self.yRes=self.res[1]
          self.zRes=self.res[2]
          place=0
          for x in np.arange(self.minX,self.maxX,self.xRes):
            for y in np.arange(self.minY,self.maxY,self.yRes):
              for z in np.arange(self.minZ,self.maxZ,self.zRes):
                self.x[place]=x
                self.y[place]=y
                self.z[place]=z
                place=place+1
          self.gapIn=np.full((self.nVox,nScans),-1,dtype=float)
          self.gapTo=np.full((self.nVox,nScans),-1,dtype=float)
          self.nBeam=np.full((self.nVox,nScans),-1,dtype=float)
          self.PAIb=np.full((self.nVox,nScans),-1,dtype=float)
          self.PAIrad=np.full((self.nVox,nScans),-1,dtype=float)
          self.meanRefl=np.full((self.nVox,nScans),-1,dtype=float)
          self.meanZen=np.full((self.nVox,nScans),-1,dtype=float)
          self.totVol=np.full((self.nVox,nScans),-1,dtype=float)
          self.sampVol=np.full((self.nVox,nScans),-1,dtype=float)
          self.meanRange=np.full((self.nVox,nScans),-1,dtype=float)
        self.gapIn[:,i]=np.array(f['GAP'])
        self.gapTo[:,i]=np.array(f['GAPTO'])
        self.nBeam[:,i]=np.array(f['HITS'])
        self.PAIb[:,i]=np.array(f['PAIB'])
        self.PAIrad[:,i]=np.array(f['PAIRAD'])
        self.meanRefl[:,i]=np.array(f['MEANREFL'])
        self.meanZen[:,i]=np.array(f['MEANZEN'])
        self.totVol[:,i]=np.array(f['TOTVOL'])
        self.sampVol[:,i]=np.array(f['SAMPVOL'])
        self.meanRange[:,i]=np.array(f['MEANRANGE'])
        f.close()
      else:   # ASCII file
        self.hdf=False
        temp=np.loadtxt(inNamen,unpack=True, dtype=float,comments='#',delimiter=" ")
        # if first, allocate space
        if(i==0):
          self.nVox=len(temp[0])
          self.x=temp[0]
          self.y=temp[1]
          self.z=temp[2]
          self.gapIn=np.full((self.nVox,nScans),-1,dtype=float)
          self.gapTo=np.full((self.nVox,nScans),-1,dtype=float)
          self.nBeam=np.full((self.nVox,nScans),-1,dtype=float)
          self.PAIb=np.full((self.nVox,nScans),-1,dtype=float)
          self.PAIrad=np.full((self.nVox,nScans),-1,dtype=float)
          self.meanRefl=np.full((self.nVox,nScans),-1,dtype=float)
          self.meanZen=np.full((self.nVox,nScans),-1,dtype=float)
          self.totVol=np.full((self.nVox,nScans),-1,dtype=float)
          self.sampVol=np.full((self.nVox,nScans),-1,dtype=float)
          self.meanRange=np.full((self.nVox,nScans),-1,dtype=float)
        # copy data
        baseInd=3
        self.gapIn[:,i]=temp[baseInd,:]
        self.gapTo[:,i]=temp[baseInd+3,:]
        self.nBeam[:,i]=temp[baseInd+4,:]
        self.PAIb[:,i]=temp[baseInd+6,:]
        self.PAIrad[:,i]=temp[baseInd+9,:]
        self.meanRefl[:,i]=temp[baseInd+10,:]
        self.meanZen[:,i]=temp[baseInd+11,:]
        self.totVol[:,i]=temp[baseInd+8,:]
        self.sampVol[:,i]=temp[baseInd+7,:]
      # keep count
      i=i+1


  ##############################################

  def readAllScans(self,inNamen):
    '''read all scans from a single file'''
    # is this an ASCII or a HDF5 file?
    ending=inNamen.strip().split('.')[-1].lower()
    if((ending=="h5")|(ending=="hdf")):  # HDF5
      self.hdf=True
      f=h5py.File(inNamen.strip(),'r')
      nScans=list(f['NSCANS'])[0]
      self.nX=list(f['NX'])[0]
      self.nY=list(f['NY'])[0]
      self.nZ=list(f['NZ'])[0]
      self.nVox=self.nX*self.nY*self.nZ
      self.gapIn=np.reshape(np.array(f['GAP']),(self.nVox,nScans))
      self.gapTo=np.reshape(np.array(f['GAPTO']),(self.nVox,nScans))
      self.nBeam=np.reshape(np.array(f['HITS']),(self.nVox,nScans))
      self.PAIb=np.reshape(np.array(f['PAIB']),(self.nVox,nScans))
      self.PAIrad=np.reshape(np.array(f['PAIRAD']),(self.nVox,nScans))
      self.meanRefl=np.reshape(np.array(f['MEANREFL']),(self.nVox,nScans))
      self.meanZen=np.reshape(np.array(f['MEANZEN']),(self.nVox,nScans))
      self.totVol=np.reshape(np.array(f['TOTVOL']),(self.nVox,nScans))
      self.sampVol=np.reshape(np.array(f['SAMPVOL']),(self.nVox,nScans))
      self.meanRange=np.reshape(np.array(f['MEANRANGE']),(self.nVox,nScans))
      self.x=np.full((self.nVox,nScans),-1,dtype=float)
      self.y=np.full((self.nVox,nScans),-1,dtype=float)
      self.z=np.full((self.nVox,nScans),-1,dtype=float)
      self.bounds=np.array(f['BOUNDS'])
      self.minX=self.bounds[0]
      self.minY=self.bounds[1]
      self.minZ=self.bounds[2]
      self.maxX=self.bounds[3]
      self.maxY=self.bounds[4]
      self.maxZ=self.bounds[5]
      self.res=np.array(f['RES'])
      self.xRes=self.res[0]
      self.yRes=self.res[1]
      self.zRes=self.res[2]
      place=0
      for i in range(0,self.nX):
        for j in range(0,self.nY):
          for k in range(0,self.nZ):
            place=i+j*self.nX+k*self.nX*self.nY
            self.x[place]=self.minX+i*self.xRes
            self.y[place]=self.minY+j*self.yRes
            self.z[place]=self.minZ+k*self.zRes
      f.close()
    else:    # ASCII
      self.hdf=False
      # read in to array and get 1D arrays
      temp=np.loadtxt(inNamen,unpack=True, dtype=float,comments='#',delimiter=" ")
      self.x=temp[0]
      self.y=temp[1]
      self.z=temp[2]
      nScans=int((temp.shape[0]-3)/12)
      self.nVox=len(self.x)
      # create blank arrays
      self.gapIn=np.full((self.nVox,nScans),-1,dtype=float)
      self.gapTo=np.full((self.nVox,nScans),-1,dtype=float)
      self.nBeam=np.full((self.nVox,nScans),-1,dtype=float)
      self.PAIb=np.full((self.nVox,nScans),-1,dtype=float)
      self.PAIrad=np.full((self.nVox,nScans),-1,dtype=float)
      self.meanRefl=np.full((self.nVox,nScans),-1,dtype=float)
      self.meanZen=np.full((self.nVox,nScans),-1,dtype=float)
      self.totVol=np.full((self.nVox,nScans),-1,dtype=float)
      self.sampVol=np.full((self.nVox,nScans),-1,dtype=float)
      self.meanRange=np.full((self.nVox,nScans),-1,dtype=float)
      # loop over scans
      for i in range(0,nScans):
        baseInd=i*12+3
        self.gapIn[:,i]=temp[baseInd,:]
        self.gapTo[:,i]=temp[baseInd+3,:]
        self.nBeam[:,i]=temp[baseInd+4,:]
        self.PAIb[:,i]=temp[baseInd+6,:]
        self.PAIrad[:,i]=temp[baseInd+9,:]
        self.meanRefl[:,i]=temp[baseInd+10,:]
        self.meanZen[:,i]=temp[baseInd+11,:]
        self.totVol[:,i]=temp[baseInd+8,:]
        self.sampVol[:,i]=temp[baseInd+7,:]


  ##############################################

  def combineScans(self,minGap,outNamen,writeHDF):
    '''combine scans, making decisons on the way'''
    # open output and write header
    if(writeHDF==False):   # ASCII header
      f=open(outNamen,'w')
      line="# 1 x, 2 y, 3 z, 4 mGapTo, 5 maxGapTo, 6 mGapIn, 7 mPAIb, 8 mPAIrad, 9 mPAIw, 10 totBeam, 11 bPAIb\n"
      f.write(line)
    else:     # HDF5 header
      f=h5py.File(outNamen,'w')
      f.create_dataset('NSCANS',(1,),data=1,dtype='i4')
      f.create_dataset('NX',(1,),data=self.nX,dtype='i4')
      f.create_dataset('NY',(1,),data=self.nY,dtype='i4')
      f.create_dataset('NZ',(1,),data=self.nZ,dtype='i4')
      f.create_dataset('RES',data=self.res,dtype='f4')
      f.create_dataset('BOUNDS',data=self.bounds,dtype="f8")
      # space to store data
      arrmGapTo=np.full((1,self.nVox),-1,dtype='f4')
      arrmaxGapTo=np.full((1,self.nVox),-1,dtype='f4')
      arrmGapIn=np.full((1,self.nVox),-1,dtype='f4')
      arrmPAIb=np.full((1,self.nVox),-1,dtype='f4')
      arrtotBeam=np.full((1,self.nVox),-1,dtype='f4')
      arrmPAIrad=np.full((1,self.nVox),-1,dtype='f4')
      arrtotVol=np.full((1,self.nVox),-1,dtype='f4')
      arrsampVol=np.full((1,self.nVox),-1,dtype='f4')
      arrmeanRange=np.full((1,self.nVox),-1,dtype='f4')
    # loop over voxels
    for i in range(0,len(self.x)):
      if((i%int(len(self.x)/100))==0):
        print("Voxel",int(i/int(len(self.x)/100)),"%")
      useInd=np.where(self.gapTo[i,:]>=minGap)[0]
      if(len(useInd)>0):  # usable scans?
        if(np.max(self.gapTo[i,useInd])>0.0):
          # mean across all usable scan locations, weighted by visibility
          mGapTo=np.average(self.gapTo[i,useInd],weights=self.gapTo[i,useInd])
          maxGapTo=np.max(self.gapTo[i,useInd])
          mGapIn=np.average(self.gapIn[i,useInd],weights=self.gapTo[i,useInd])
          if(np.sum(self.sampVol[i,useInd])>0.0):
            mPAIb=np.average(self.PAIb[i,useInd],weights=self.sampVol[i,useInd])
          else:
            mPAIb=-1.0
          mPAIrad=np.average(self.PAIrad[i,useInd],weights=self.gapTo[i,useInd])
          totBeam=int(np.sum(self.nBeam[i,useInd]))
          # fit a function for waveform based PAI
          mPAIw=self.fitLAD(i,useInd) #gapTo[i,useInd],meanRefl[i,useInd],meanZen[i,useInd])
          totVol=np.sum(self.totVol[i,useInd])
          sampVol=np.sum(self.sampVol[i,useInd])
          meanRange=np.min(self.meanRange[i,useInd])
        elif(np.max(self.gapTo[i,useInd])<=-2.0):
          mGapIn=mGapTo=maxGapTo=mPAIb=mPAIrad=mPAIw=-2
          totBeam=totVol=sampVol=meanRange=0
        else:
          mGapIn=mGapTo=maxGapTo=mPAIb=mPAIrad=mPAIw=-1
          totBeam=totVol=sampVol=meanRange=0
      elif(np.max(self.gapTo[i,:])<=-2.0):
        mGapIn=mGapTo=maxGapTo=mPAIb=mPAIrad=mPAIw=-2
        totBeam=totVol=sampVol=meanRange=0
      else:
        mGapIn=mGapTo=maxGapTo=mPAIb=mPAIrad=mPAIw=-1
        totBeam=totVol=sampVol=meanRange=0
      # Beland's maximum visibility PAI
      bPAIb=self.PAIb[i,np.argmax(self.gapTo[i,:])]
      # write to file
      if(writeHDF==False):   # ASCII data
        line=str(self.x[i])+" "+str(self.y[i])+" "+str(self.z[i])+" "+str(mGapTo)+" "+str(maxGapTo)+" "+str(mGapIn)+" "+str(mPAIb)+" "+str(mPAIrad)+" "+str(mPAIw)+" "+str(totBeam)+" "+str(bPAIb)+" "+str(meanRange)+"\n"
        f.write(line)
      else:
        arrmGapTo[0,i]=mGapTo
        arrmaxGapTo[0,i]=maxGapTo
        arrmGapIn[0,i]=mGapIn
        arrmPAIb[0,i]=mPAIb
        arrtotBeam[0,i]=totBeam
        arrmPAIrad[0,i]=mPAIrad
        arrtotVol[0,i]=totVol
        arrsampVol[0,i]=sampVol
        arrmeanRange[0,i]=meanRange
    if(writeHDF):  # write the HDF data
      f.create_dataset('GAP',data=arrmGapIn,dtype='f4',compression="gzip")
      f.create_dataset('GAPTO',data=arrmGapTo,dtype='f4',compression="gzip")
      f.create_dataset('MAXGAP',data=arrmaxGapTo,dtype='f4',compression="gzip")
      f.create_dataset('HITS',data=arrtotBeam,dtype='f4',compression="gzip")
      f.create_dataset('PAIB',data=arrmPAIb,dtype='f4',compression="gzip")
      f.create_dataset('PAIRAD',data=arrmPAIrad,dtype='f4',compression="gzip")
      f.create_dataset('TOTVOL',data=arrtotVol,dtype='f4',compression="gzip")
      f.create_dataset('SAMPVOL',data=arrsampVol,dtype='f4',compression="gzip")
      f.create_dataset('MEANRANGE',data=arrmeanRange,dtype='f4',compression="gzip")
    f.close()
    print("Written to",outNamen)
    return


  ##############################################

  def fitLAD(self,i,useInd):# gapTo,meanRefl,meanZen):
    '''fit leaf angle distribution'''
    mPAIw=-1  # blank for now
    return(mPAIw)


##############################################
# main
##############################################

if __name__ == '__main__':
  # read command line
  cmdargs = getCmdArgs()
  # transfer to local variables, to ease debugging
  inNamen=cmdargs.inNamen
  listNamen=cmdargs.listNamen
  minGap=cmdargs.minGap
  outNamen=cmdargs.outNamen

  # read data
  voxels=voxelData(inNamen,listNamen)

  # combine data
  voxels.combineScans(minGap,outNamen,cmdargs.writeHDF)


##############################################
# end
##############################################

