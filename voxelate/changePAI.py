import numpy as np
import argparse
from sys import exit
import csv



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
    p.add_argument("--changeFile", dest ="changeNamen",default=None,type=str, help=("Detect change against this file"))
    p.add_argument("--output", dest ="outNamen",default="changeVox.vox", type=str,help=("Outpt filename\nDefault = changeVox.vox"))
    p.add_argument("--minGap", dest ="minGap", type=float, default=0.0, help=("Minimum trustable gap fraction to a voxel\nDefault = 0"))
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
    self.readVoxel(inNamen)
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

  def readVoxel(self,inNamen):
    '''read all scans from a single file'''
    # read in to array and get 1D arrays
    temp=np.loadtxt(inNamen,unpack=True, dtype=float,comments='#',delimiter=" ")
    self.x=temp[0]
    self.y=temp[1]
    self.z=temp[2]
    self.gapTo=temp[3]
    self.gapIn=temp[5]
    self.maxGapTo=temp[4]
    self.nBeam=temp[9]


  ##############################################

  def findChange(self,cVox,outNamen):
    '''Detect change between two voxel files'''
    # open output
    f=open(outNamen,'w')
    f.write("# 1 x, 2 y, 3 z, 4 diff, 5 gapTo1, 6 gapTo2, 7 gapIn1, 8 gapIn2\n")
    # loop over voxels
    for i in range(0,self.nX):
      for j in range(0,self.nY):
        for k in range(0,self.nZ):
          p=k*self.nY*self.nX+j*self.nX+i
          # find nearest other voxel
          xBin=int((self.x[p]-cVox.minX)/cVox.xRes)
          yBin=int((self.y[p]-cVox.minY)/cVox.yRes)
          zBin=int((self.z[p]-cVox.minZ)/cVox.zRes)
          # does this voxel exist?
          if((xBin>=0)&(xBin<cVox.nX)&(yBin>=0)&(yBin<cVox.nY)&(zBin>=0)&(zBin<cVox.nZ)):
            cp=zBin*cVox.nY*cVox.nX+yBin*cVox.nX+xBin
            # find and write difference
            diff=cVox.gapIn[cp]-self.gapIn[p]
            line=str(self.x[p])+" "+str(self.y[p])+" "+str(self.z[p])+" "+str(diff)+" "+str(self.maxGapTo[p])\
                               +" "+str(cVox.maxGapTo[cp])+" "+str(self.gapIn[p])+" "+str(cVox.gapIn[cp])+"\n"
            f.write(line)
    print("Written to ",outNamen)


##############################################
# main
##############################################

if __name__ == '__main__':
  # read command line
  cmdargs = getCmdArgs()
  # transfer to local variables, to ease debugging
  inNamen=cmdargs.inNamen
  minGap=cmdargs.minGap
  outNamen=cmdargs.outNamen
  changeNamen=cmdargs.changeNamen

  # read data
  voxels=voxelData(inNamen,None)
  changeVox=voxelData(changeNamen,None)

  # test change between two files
  voxels.findChange(changeVox,outNamen)


##############################################
# end
##############################################

