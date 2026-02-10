#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libDEMhandle.h"
#include "libLasProcess.h"
#include "libLidVoxel.h"
#include "tiffRead.h"
#include "libLidarHDF.h"




/*#############################*/
/*# Voxelises TLS data       #*/
/*# Reads readRXP.cpp output #*/
/*# S Hancock, 2016          #*/
/*############################*/


/*#######################################*/
/*# Copyright 2014-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of voxelate, a library for voxelising lidar.       #*/
/*#                                                                      #*/
/*# voxelate is free software: you can redistribute it and/or modify     #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# voxelate is distributed in the hope that it will be useful,          #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with voxelate.  If not, see <http://www.gnu.org/licenses/>. #*/
/*########################################################################*/


/*##################################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  int nScans;        /*number of scans*/
  char **inList;     /*input filename list*/
  char outNamen[200];/*output filename*/
  char writeHDF;     /*write output as a HDF5 file*/
  /*beam description*/
  double origin[3];  /*origin of beam*/
  float grad[3];     /*beam vector*/
  char doVoxCov;     /*cover for all voxels switch*/
  char gaussFoot;    /*Gaussian footprint switch*/
  float fSigma;      /*footprint sigma if Gaussian, width of not*/
  float pSigma;      /*pulse width*/
  float rRes;        /*range resolution*/
  int nBins;         /*number of bins*/
  int nWaves;        /*number of waves*/
  float buffLen;     /*buffer length*/
  /*lidar parameters for voxelisation*/
  lidVoxPar lidPar;
  /*switches*/
  char writeWave;    /*write waveform switch*/
  char readBounds;   /*read bounds from file rather that prescribing*/
  char silouhetteVox;/*silouhette voxels or gap voxels*/
  char onlyPrintFilled; /*only print out filled voxel coords*/
  /*options*/
  double bounds[6];  /*area of interest, minX minY minZ, maxX, maxY maxZ*/
  float vRes[3];     /*voxel resolution in each axis*/
  double maxZen;     /*maximum absolute zenith to allow. To filter tilt mount*/
  /*dem*/
  demStruct *dem;    /*structure to hold a DEM*/
  char demNamen[200];/*DEM filename*/
  char useDEM;       /*use a DEM switch*/
  float demTol;      /*DEM tolerance (if point this close to DEM, count as DEM*/
}control;


/*##################################################*/
/*data for a HDF5 file*/

typedef struct{
  float *gapArr;
  float *gapToArr;
  float *PAIbArr;
  float *PAIradArr;
}hdfVox;


/*##################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  tlsScan *scan=NULL;
  voxStruct *vox=NULL;
  tlsVoxMap map;
  float **waves=NULL;
  float **tlsBeamWave(tlsScan *,voxStruct *,tlsVoxMap *,control *);
  void calculateAllVoxCov(tlsScan *,voxStruct *,tlsVoxMap *,control *,int);
  void writeWaveform(float **,int,char *,int);
  void footprintBounds(control *);
  void simpleGapFrac(voxStruct *,control *);


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*set some values blank*/
  map.mapFile=NULL;
  map.mapPoint=NULL;
  map.nIn=NULL;

  /*read voxel bounds from data if needed. For beams, z bounds are read here*/
  fprintf(stdout,"Determining bounds\n");
  if(dimage->readBounds||(!dimage->doVoxCov))readBoundsFromTLS(&(dimage->bounds[0]),dimage->inList,dimage->nScans);
  if(!dimage->doVoxCov){  /*set x and y bounds from beam. Z bounds are set from TLS data above*/
    footprintBounds(dimage);
    beamVoxelBounds(&(dimage->origin[0]),&(dimage->grad[0]),dimage->fSigma,dimage->gaussFoot,&(dimage->bounds[0]));
    dimage->bounds[2]-=dimage->buffLen;
    dimage->bounds[5]+=dimage->buffLen;
    dimage->origin[2]=dimage->bounds[5];
  }
  fprintf(stdout,"Determined bounds\n");

  /*allocate voxel structure*/
  vox=voxAllocate(dimage->nScans,dimage->vRes,dimage->bounds,0);
  if(dimage->doVoxCov&&(!dimage->silouhetteVox))vox->savePts=0;
  vox->maxZen=dimage->maxZen;

  /*read DEM if needed*/
  if(dimage->useDEM){
    dimage->dem=readTifDEM(dimage->demNamen,dimage->bounds[0],dimage->bounds[1],dimage->bounds[3],dimage->bounds[4]);
    vox->dem=dimage->dem;
    vox->demTol=dimage->demTol;
  }

  /*loop over TLS files and add up gaps to voxel map*/
  for(i=0;i<dimage->nScans;i++){
    /*read TLS data within and populate voxel gaps. Save points if required*/
    scan=readOneTLS(dimage->inList[i],vox,0,&map,i,&dimage->lidPar);

    /*if single waveform simulation*/
    if(!dimage->doVoxCov){
      fprintf(stderr,"This needs replacing with a function that updates the siloughette image\n");
      exit(1);
      //calculateAllVoxCov(scan,vox,&map,dimage,i);
      /*make waveform*/
      //waves=tlsBeamWave(scan,vox,&map,dimage);
      /*output waveform*/
      //writeWaveform(waves,dimage->nBins,dimage->outNamen,dimage->nWaves);
    }/*all voxels or one beam switch*/

    /*tidy up*/
    scan=tidyTLScans(scan,1);
  }/*scan loop*/

  /*do simple gap fraction if needed*/
  if(dimage->doVoxCov&&(!dimage->silouhetteVox))simpleGapFrac(vox,dimage);

  /*tidy up*/
  TTIDY((void **)waves,dimage->nWaves);
  tidyVoxelMap(&map,vox->nVox);
  vox=tidyVox(vox);
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nScans);
    dimage->dem=tidyDEMstruct(dimage->dem);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##################################################*/
/*simple voxel gap fraction*/

void simpleGapFrac(voxStruct *vox,control *dimage)
{
  int i=0,j=0,k=0;
  int vPlace=0,n=0;
  int nPer=0;   /*number of parameters per voxel*/
  float gap=0,PAIb=0,gapTo=0;
  float thisGap=0,PAIrad=0;
  double origin[3],topZ=0;
  char checkFilled(float **,int,int);
  char writeCoord=0,underGround=0;
  char checkDTM(double,double,double,double,double,demStruct *);
  FILE *opoo=NULL;
  hdfVox voxData;  /*to hold data for HDF5*/
  void writeHDFvoxels(hdfVox *,voxStruct *,control *);


  /*ASCII or HDF5?*/
  if(dimage->writeHDF==0){  /*ASCII*/
    /*open file*/
    if((opoo=fopen(dimage->outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
      exit(1);
    }

    /*write the header*/
    nPer=13;   /*parameters per voxel*/
    if(!dimage->onlyPrintFilled){
      fprintf(opoo,"# 1 x, 2 y, 3 z");
      for(n=0;n<vox->nScans;n++)fprintf(opoo,", %d gap%d, %d inHit%d, %d inMiss%d, %d gapTo%d, %d hits%d, %d miss%d, %d PAIbel%d, %d volSamp%d, %d volTrans%d, %d PAIrad%d, %d meanRefl%d, %d meanZen%d, %d meanRange%d",nPer*n+4,n,nPer*n+5,n,nPer*n+6,n,nPer*n+7,n,nPer*n+8,n,nPer*n+9,n,nPer*n+10,n,nPer*n+11,n,nPer*n+12,n,nPer*n+13,n,nPer*n+14,n,nPer*n+15,n,nPer*n+16,n);
      fprintf(opoo,"\n");
    }
    voxData.gapArr=NULL;
    voxData.gapToArr=NULL;
    voxData.PAIbArr=NULL;
    voxData.PAIradArr=NULL;
  }else{   /*HDF5*/
    /*allocate space for data*/
    voxData.gapArr=falloc(vox->nScans*vox->nX*vox->nY*vox->nZ,"gapArr",0);
    voxData.gapToArr=falloc(vox->nScans*vox->nX*vox->nY*vox->nZ,"gapToArr",0);
    voxData.PAIbArr=falloc(vox->nScans*vox->nX*vox->nY*vox->nZ,"PAIbArr",0);
    voxData.PAIradArr=falloc(vox->nScans*vox->nX*vox->nY*vox->nZ,"PAIradArr",0);
  }/*format check*/


  /*loop over voxels*/
  for(i=0;i<vox->nX;i++){
    origin[0]=((double)i+0.5)*vox->res[0]+vox->bounds[0];
    for(j=0;j<vox->nY;j++){
      origin[1]=((double)j+0.5)*vox->res[1]+vox->bounds[1];
      for(k=0;k<vox->nZ;k++){  /*z axis is the right way up*/
        vPlace=k*vox->nX*vox->nY+j*vox->nX+i;
        origin[2]=vox->bounds[2]+(double)k*vox->res[2];

        /*determine if underground*/
        if(dimage->useDEM){
          topZ=origin[2]+(double)k*vox->res[2];
          underGround=checkDTM(topZ,origin[0]+(double)vox->res[0]/2.0,origin[1]+(double)vox->res[1]/2.0,(double)vox->res[0],(double)vox->res[1],dimage->dem);
        }else underGround=0;

        /*if we are writing all*/
        if(dimage->onlyPrintFilled&&(!dimage->writeHDF))writeCoord=checkFilled(vox->inHit,vox->nScans,vPlace);
        else                                            writeCoord=1;

        /*do we write this coord?*/
        if(writeCoord){
          if(dimage->writeHDF==0)fprintf(opoo,"%.3f %.3f %.3f",origin[0],origin[1],origin[2]);

          if(!dimage->onlyPrintFilled){
            /*loop over scans in this voxel*/
            for(n=0;n<vox->nScans;n++){
              /*are we underGround?*/
              if(!underGround){
                /*gap fraction to voxel*/
                if((vox->miss[n][vPlace]+vox->hits[n][vPlace])>0.0)gapTo=vox->hits[n][vPlace]/(vox->hits[n][vPlace]+vox->miss[n][vPlace]);
                else                                               gapTo=-1.0;

                /*gap fraction within voxel*/
                if((vox->inMiss[n][vPlace]+vox->inHit[n][vPlace])>0.0)gap=vox->inMiss[n][vPlace]/(vox->inMiss[n][vPlace]+vox->inHit[n][vPlace]);
                else                                                  gap=-1.0;

                /*pai within voxel*/
                if(vox->sampVol[n][vPlace]>0.0)PAIb=((float)vox->inHit[n][vPlace]/(vox->totVol[n][vPlace]))/0.6;  /*0.6 here is the Ross-G function for random orientation*/
                else                           PAIb=-1.0;
                if(gapTo>0.0){
                  thisGap=(gapTo<dimage->lidPar.minGap)?dimage->lidPar.minGap:gapTo;
                  PAIrad=(4.0*M_PI*(vox->sumRsq[n][vPlace]/(thisGap*thisGap))/2.0)/vox->volume;
                }else PAIrad=-1.0;
              }else{   /*we are underground*/
                gapTo=gap=PAIb=PAIrad=-2.0;
                vox->miss[n][vPlace]=vox->hits[n][vPlace]=vox->inMiss[n][vPlace]=vox->inHit[n][vPlace]=vox->sampVol[n][vPlace]=-2.0;
              }

              /*LAI from fitting to LAD, if we have sufficient data*/

              if(dimage->writeHDF==0){  /*write to ASCII file*/
                /*write results*/
                fprintf(opoo," %f %f %f %f %f %f %f",gap,vox->inHit[n][vPlace],vox->inMiss[n][vPlace],\
                                                 gapTo,vox->hits[n][vPlace],vox->miss[n][vPlace],PAIb);
                fprintf(opoo," %f %f %f %f %f %f",vox->sampVol[n][vPlace],vox->totVol[n][vPlace],PAIrad,\
                                    vox->meanRefl[n][vPlace],vox->meanZen[n][vPlace],vox->meanRange[n][vPlace]);
              }else{   /*save to HDF5 arrays*/
                voxData.gapArr[n*vox->nVox+vPlace]=gap;
                voxData.gapToArr[n*vox->nVox+vPlace]=gapTo;
                voxData.PAIbArr[n*vox->nVox+vPlace]=PAIb;
                voxData.PAIradArr[n*vox->nVox+vPlace]=PAIrad;
              }
            }/*scan loop*/
          }/*write all scan info check*/
          if(dimage->writeHDF==0)fprintf(opoo,"\n");
        }/*write point check*/
      }/*z loop*/
    }/*y loop*/
  }/*x loop*/

  /*close up*/
 if(dimage->writeHDF==0){  /*ASCII*/
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",dimage->outNamen);
  }else{   /*HDF5*/
    /*write arrays*/
    writeHDFvoxels(&voxData,vox,dimage);

    /*tidy up*/
    TIDY(voxData.gapArr);
    TIDY(voxData.gapToArr);
    TIDY(voxData.PAIbArr);
    TIDY(voxData.PAIradArr);
  }

  return;
}/*simpleGapFrac*/



/*##################################################*/
/*check if voxel is underground*/

char checkDTM(double topZ,double x,double y,double xRes,double yRes,demStruct *dem)
{
  int i=0,j=0,place=0;
  int i0=0,i1=0;
  int j0=0,j1=0;
  double minT=0;
  char underGround=0,tested=0;

  /*determine DTM pixel*/
  i0=(int)(((x-xRes/2.0)-dem->minX)/(double)dem->res);
  i1=(int)(((x+xRes/2.0)-dem->minX)/(double)dem->res);
  j0=(int)(((y-yRes/2.0)-dem->minY)/(double)dem->res);
  j1=(int)(((y+yRes/2.0)-dem->minY)/(double)dem->res);

  if(i0<0)            i0=0;
  else if(i0>=dem->nX)i0=dem->nX-1;
  if(i1<0)            i1=0;
  else if(i1>=dem->nX)i1=dem->nX-1;
  if(j0<0)            j0=0;
  else if(j0>=dem->nY)j0=dem->nY-1;
  if(j1<0)            j1=0;
  else if(j1>=dem->nY)j1=dem->nY-1;

  /*loop over intersecting DTM pixels*/
  tested=0;
  minT=100000000.0;
  if(topZ<dem->maxZ){
    for(i=i0;i<=i1;i++){
      for(j=j0;j<=j1;j++){
        place=j*dem->nX+i;
        if(dem->z[place]>(dem->noData+0.0000001)){
          if(dem->z[place]<minT)minT=dem->z[place];
          tested=1;
        }
      }
    }
  }
  if(tested&&(topZ<minT))underGround=1;
  else                   underGround=0;

  return(underGround);
}/*checkDTM*/


/*##################################################*/
/*write voxels in HDF5 format*/

void writeHDFvoxels(hdfVox *voxData,voxStruct *vox,control *dimage)
{
  int nVox;
  hid_t file;         /* Handles */
  float *temp=NULL;
  float *flattenArr(float **,int,int);

  /*to save calculations*/
  nVox=vox->nX*vox->nY*vox->nZ;

  /*open file*/
  file=H5Fcreate(dimage->outNamen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  /*write the header*/
  write1dIntHDF5(file,"NSCANS",&vox->nScans,1);
  write1dIntHDF5(file,"NX",&vox->nX,1);
  write1dIntHDF5(file,"NY",&vox->nY,1);
  write1dIntHDF5(file,"NZ",&vox->nZ,1);
  write1dDoubleHDF5(file,"BOUNDS",&(vox->bounds[0]),6);
  write1dDoubleHDF5(file,"RES",&(vox->res[0]),3);

  /*write the data, with gzip compression*/
  writeComp2dFloatHDF5(file,"GAP",voxData->gapArr,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"GAPTO",voxData->gapToArr,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"PAIB",voxData->PAIbArr,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"PAIRAD",voxData->PAIradArr,vox->nScans,nVox);

  /*flatten 2D arrays*/
  temp=flattenArr(vox->totVol,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"TOTVOL",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->sampVol,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"SAMPVOL",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->hits,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"HITS",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->miss,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"MISS",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->inHit,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"INHIT",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->inMiss,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"INMISS",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->meanRefl,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"MEANREFL",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->meanZen,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"MEANZEN",temp,vox->nScans,nVox);
  TIDY(temp);
  temp=flattenArr(vox->meanRange,vox->nScans,nVox);
  writeComp2dFloatHDF5(file,"MEANRANGE",temp,vox->nScans,nVox);
  TIDY(temp);


  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }
  fprintf(stdout,"Voxels written to %s\n",dimage->outNamen);
  return;
}/*writeHDFvoxels*/


/*##################################################*/
/*flatten an array for HDF format*/

float *flattenArr(float **jimlad,int nX,int nY)
{
  int i=0,j=0;
  float *temp=NULL;

  temp=falloc(nX*nY,"temporary",0);
  for(i=0;i<nX;i++){
    for(j=0;j<nY;j++)temp[i*nY+j]=jimlad[i][j];
  }

  return(temp);
}/*flattenArr*/


/*##################################################*/
/*see if a voxel is filled*/

char checkFilled(float **inHit,int nScans,int vPlace)
{
  int i=0;
  char isFilled=0;

  /*loop over scans*/
  for(i=0;i<nScans;i++){
    if(inHit[i][vPlace]>0.5){  /*inHit is a float. Use 0.5 for rounding*/
      isFilled=1;
      break;
    }
  }/*scan loop*/

  return(isFilled);
}/*checkFilled*/


/*##################################################*/
/*horizontal voxel bounds from footprint*/

void footprintBounds(control *dimage)
{
  double dz=0;
  double x0=0,y0=0;
  double minZ=0,maxZ=0;
  float rad=0,t=0,zen=0;

  minZ=dimage->bounds[2];
  maxZ=dimage->bounds[5];
  dimage->bounds[0]=dimage->bounds[1]=100000000.0;
  dimage->bounds[3]=dimage->bounds[4]=-100000000.0;

  /*determine radius*/
  zen=atan2(sqrt(dimage->grad[0]*dimage->grad[0]+dimage->grad[1]*dimage->grad[1]),dimage->grad[2]);
  if(dimage->gaussFoot)rad=determineGaussSep(dimage->fSigma,0.01)*cos(zen);
  else                 rad=dimage->fSigma*cos(zen);

  /*top*/
  dz=dimage->origin[2]-maxZ;
  t=(float)dz/dimage->grad[2];
  x0=(double)(t*dimage->grad[0])+dimage->origin[0];
  y0=(double)(t*dimage->grad[1])+dimage->origin[1];
  if((x0-rad)<dimage->bounds[0])dimage->bounds[0]=x0-rad;
  if((x0+rad)>dimage->bounds[3])dimage->bounds[3]=x0+rad;
  if((y0-rad)<dimage->bounds[1])dimage->bounds[1]=y0-rad;
  if((y0+rad)>dimage->bounds[4])dimage->bounds[4]=y0+rad;

  /*bottom*/
  dz=dimage->origin[2]-minZ;
  t=(float)dz/dimage->grad[2];
  x0=(double)(t*dimage->grad[0])+dimage->origin[0];
  y0=(double)(t*dimage->grad[1])+dimage->origin[1];
  if((x0-rad)<dimage->bounds[0])dimage->bounds[0]=x0-rad;
  if((x0+rad)>dimage->bounds[3])dimage->bounds[3]=x0+rad;
  if((y0-rad)<dimage->bounds[1])dimage->bounds[1]=y0-rad;
  if((y0+rad)>dimage->bounds[4])dimage->bounds[4]=y0+rad;

  return;
}/*footprintBounds*/


/*##################################################*/
/*write out waveforms*/

void writeWaveform(float **waves,int nBins,char *outNamen,int nWaves)
{
  int i=0,j=0;
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }

  fprintf(opoo,"# 1 elevation, 2 waveform, 3 true profile, 4 pulse smoothed profile, 5 mean gap to voxel, 6 max gap to voxel from a given scan\n");
  for(i=0;i<nBins;i++){
    for(j=0;j<nWaves;j++)fprintf(opoo,"%f ",waves[j][i]);
    fprintf(opoo,"\n");
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);
  return;
}/*writeWaveform*/


/*##################################################*/
/*calculate cover in all voxels*/

void calculateAllVoxCov(tlsScan *scans,voxStruct *vox,tlsVoxMap *map,control *dimage, int fInd)
{
  int i=0,j=0,k=0;
  int vPlace=0,n=0;
  int *voxList=NULL;
  rImageStruct *rImage=NULL;
  float **tlsWave=NULL;
  float gap=0;
  double origin[3],bounds[6];
  FILE *opoo=NULL;

  /*open file*/
  if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  }

  /*write header*/
  if(!dimage->silouhetteVox){
    fprintf(opoo,"# 1 x, 2 y, 3 z");
    for(n=0;n<vox->nScans;n++){
      fprintf(opoo,", %d gap%d, %d inHit%d, %d inMiss%d, %d hits%d, %d miss%d",n*5+4,n+1,n*5+5,n+1,n*5+6,n+1,n*5+7,n+1,n*5+8,n+1);
    }
    fprintf(opoo,"\n");
  }else{

  }

  /*intersecting voxel list. Dummy for here*/
  voxList=ialloc(vox->nZ,"voxList",0);
  tlsWave=fFalloc(dimage->nWaves,"tls waveform",0);   /*this is multiple waves now*/
  for(i=0;i<dimage->nWaves;i++)tlsWave[i]=falloc((uint64_t)vox->nZ,"tls waveform",i+1);
  for(i=0;i<vox->nZ;i++)tlsWave[0][i]=vox->bounds[5]-(float)i*vox->res[2];

  origin[2]=vox->bounds[5];
  bounds[2]=vox->bounds[2]-(double)dimage->buffLen;
  bounds[5]=vox->bounds[5]+(double)dimage->buffLen;
  for(i=0;i<vox->nX;i++){
    origin[0]=((double)i+0.5)*vox->res[0]+vox->bounds[0];
    bounds[0]=(double)i*vox->res[0]+vox->bounds[0];
    bounds[3]=bounds[0]+vox->res[0];
    for(j=0;j<vox->nY;j++){
      origin[1]=((double)j+0.5)*vox->res[1]+vox->bounds[1];
      bounds[1]=(double)j*vox->res[1]+vox->bounds[1];
      bounds[4]=bounds[1]+vox->res[1];


      /*check mode*/
      if(dimage->silouhetteVox){  /*trace through point cloud*/
        rImage=allocateRangeImage(vox->res[0]/2.0,vox->res[2],0.005,NULL,&(origin[0]),&(bounds[0]));
fprintf(stdout,"bins %d %d\n",rImage->nBins,vox->nZ);
exit(1);

        /*make gap images*/
        for(k=0;k<vox->nZ;k++)voxList[k]=k*vox->nX*vox->nY+j*vox->nX+i;
        silhouetteImage(dimage->nScans,NULL,scans,rImage,&dimage->lidPar,voxList,vox->nZ,map);

        /*file in the last, as we're looking down*/
        fillInRimageGround(rImage);

        /*calculate waveform*/
        waveFromImage(rImage,&(tlsWave[1]),0,100000.0);

        for(k=0;k<vox->nZ;k++){   /*note that the z axis is reflected in this mode*/
          if((tlsWave[1][k]+tlsWave[2][k])>0)fprintf(opoo,"%f %f %f %f %f\n",origin[0],origin[1],origin[2]-(double)k*rImage->rRes,tlsWave[1][k],tlsWave[2][k]);
        }
      }else{   /*output simple gap fraction*/
        for(k=0;k<vox->nZ;k++){  /*z axis is the right way up*/
          fprintf(opoo,"%f %f %f",origin[0],origin[1],vox->bounds[2]+(double)k*vox->res[2]);
          vPlace=k*vox->nX*vox->nY+j*vox->nX+i;
          for(n=0;n<vox->nScans;n++){
            if((vox->inMiss[n][vPlace]+vox->inHit[n][vPlace])>0.0)gap=vox->inMiss[n][vPlace]/(vox->inMiss[n][vPlace]+vox->inHit[n][vPlace]);
            else                                                  gap=-1.0;
            fprintf(opoo," %f %f %f %f %f",gap,vox->inHit[n][vPlace],vox->inMiss[n][vPlace],vox->hits[n][vPlace],vox->miss[n][vPlace]);
          }/*scan loop*/
          fprintf(opoo,"\n");
        }
      }

      /*tidy up*/
      if(rImage){
        TTIDY((void **)rImage->image,rImage->nBins);
        rImage->image=NULL;
        TIDY(rImage);
      }
    }/*voxel y loop*/
  }/*voxel x loop*/


  fprintf(stdout,"Output written to %s\n",dimage->outNamen);

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  TTIDY((void **)tlsWave,dimage->nWaves);
  TIDY(voxList);
  return;
}/*calculateAllVoxCov*/


/*##################################################*/
/*make waveforms from TLS data*/

float **tlsBeamWave(tlsScan *scans,voxStruct *vox,tlsVoxMap *map,control *dimage)
{
  int i=0;
  int nIn=0,*voxList=NULL;
  float **waves=NULL;
  float iRes=0;
  double *rangeList=NULL;
  double beamRad=0;
  rImageStruct *rImage=NULL;
  void meanGapProfile(float **,int,voxStruct *,double *,float);

  iRes=0.005;

  /*determine bounds*/
  if(dimage->gaussFoot)beamRad=determineGaussSep(dimage->fSigma,0.01);
  else                 beamRad=(double)dimage->fSigma;        /*circular footprint*/


  /*determine which voxels intersect*/
  fprintf(stdout,"Intersecting voxels\n");
  voxList=beamVoxels(&(dimage->grad[0]),dimage->origin[0],dimage->origin[1],dimage->origin[2],\
                     vox->bounds,vox->res,vox->nX,vox->nY,vox->nZ,&nIn,beamRad,&rangeList,vox->res[0]);
  fprintf(stdout,"Intersected %d\n",nIn);


  /*allocate image stack*/
  fprintf(stdout,"Allocating range image\n");
  rImage=allocateRangeImage(beamRad,dimage->rRes,iRes,&(dimage->grad[0]),&(dimage->origin[0]),dimage->bounds);
  fprintf(stdout,"Allocated %lu\n",(uint64_t)rImage->nX*(uint64_t)rImage->nY*(uint64_t)rImage->nBins);
  dimage->nBins=rImage->nBins;

  /*condense data into image stack*/
  fprintf(stdout,"Doing silhouette images\n");
  silhouetteImage(dimage->nScans,NULL,scans,rImage,&dimage->lidPar,voxList,nIn,map);
  fprintf(stdout,"Done silhouette images\n");

  /*file in the last, as we're looking down*/
  fillInRimageGround(rImage);

  /*calculate waveform*/
  waves=fFalloc(dimage->nWaves,"waveforms",0);
  for(i=0;i<dimage->nWaves;i++)waves[i]=falloc((uint64_t)rImage->nBins,"waves",0);
  setWaveformRange(&(waves[0][0]),rImage->z0,&(dimage->grad[0]),rImage->nBins,rImage->rRes);
  waveFromImage(rImage,&(waves[1]),dimage->gaussFoot,dimage->fSigma);

  /*smooth with pulse shape*/
  TIDY(waves[3]);
  waves[3]=smooth(dimage->pSigma,rImage->nBins,waves[1],rImage->rRes);

  /*mean gap fraction with height*/
  meanGapProfile(waves,rImage->nBins,vox,&(dimage->origin[0]),dimage->fSigma);

  /*tidy up*/
  TIDY(voxList);
  TIDY(rangeList);
  if(rImage){
    TTIDY((void **)rImage->image,rImage->nBins);
    rImage->image=NULL;
    TIDY(rImage);
  }
  return(waves);
}/*tlsBeamWave*/


/*##################################################*/
/*mean gap fraction with height*/

void meanGapProfile(float **waves,int nBins,voxStruct *vox,double *origin,float fSigma)
{
  int i=0,j=0,place=0;
  int xInd=0,yInd=0,zInd=0;
  float contN=0,weight=0;
  float scanWeight=0;
  float meanMaxGap=0;
  float gap=0;
  double x=0,y=0;

  /*loop over waveform*/
  for(i=0;i<nBins;i++){
    waves[4][i]=waves[5][i]=0.0;
    contN=0.0;

    zInd=(int)(((double)waves[0][i]-vox->bounds[2])/vox->res[2]);
    /*loop over voxels in this layer*/
    scanWeight=0.0;
    for(xInd=0;xInd<vox->nX;xInd++){
      for(yInd=0;yInd<vox->nY;yInd++){
        x=vox->bounds[0]+(double)xInd*vox->res[0];
        y=vox->bounds[1]+(double)yInd*vox->res[1];

        weight=(float)(gaussian(x,(double)fSigma,origin[0])*gaussian(y,(double)fSigma,origin[1]));
        place=zInd*vox->nX*vox->nY+yInd*vox->nX+xInd;
        meanMaxGap=0.0;
        for(j=0;j<vox->nScans;j++){
          if((vox->hits[j][place]+vox->miss[j][place])>0.0){
            gap=vox->hits[j][place]/(vox->hits[j][place]+vox->miss[j][place]);
          }else gap=0.0;
          waves[4][i]+=gap*weight;
          if(gap>meanMaxGap)meanMaxGap=gap;
          contN+=weight;
        }
        waves[5][i]+=weight*meanMaxGap;
        scanWeight+=weight;
      }/*y voxel loop*/
    }/*x voxel loop*/
    if(scanWeight>0.0)waves[5][i]/=scanWeight;
    if(contN>0.0)waves[4][i]/=contN;
  }/*bin loop*/

  return;
}/*meanGapProfile*/


/*##################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  float zen=0,az=0,tempAng=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control structure allocation.\n");
    exit(1);
  }

  /*defaults*/
  strcpy(dimage->outNamen,"teast.wave");
  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=10000000000000.0;
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=-10000000000000.0;
  dimage->vRes[0]=dimage->vRes[1]=dimage->vRes[2]=1.0;
  dimage->grad[0]=0.0;;            /*beam vector*/
  dimage->grad[1]=0.0;             /*beam vector*/
  dimage->grad[2]=-1.0;            /*beam vector*/
  dimage->doVoxCov=1;              /*do all voxels by default*/
  dimage->readBounds=0;            /*do not read bounds from data*/
  dimage->gaussFoot=0;             /*not Gaussian footprint by default*/
  dimage->fSigma=dimage->vRes[0];  /*same size as voxels by default*/
  dimage->rRes=dimage->vRes[2];    /*range resolution*/
  dimage->origin[0]=dimage->origin[1]=dimage->origin[2]=0.0;  /*origin of beam*/
  dimage->nWaves=6;
  dimage->silouhetteVox=0;
  dimage->maxZen=100000.0;         /*accept everything*/
  dimage->writeHDF=0;
  /*lidar parameters for voxelisation*/
  dimage->lidPar.minRefl=0;
  dimage->lidPar.maxRefl=1.0;
  dimage->lidPar.beamTanDiv=sin(M_PI*0.1/180.0)/cos(M_PI*0.1/180.0);
  dimage->lidPar.beamRad=0.001;
  dimage->lidPar.appRefl=0.5;
  dimage->lidPar.minGap=0.1;
  dimage->lidPar.correctR=0;
  /*switches*/
  dimage->onlyPrintFilled=0;
  /*dem*/
  dimage->useDEM=0;
  dimage->dem=NULL;
  dimage->demTol=0.2;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nScans);
        dimage->inList=readInList(&dimage->nScans,argv[++i]);
      }else if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->nScans);
        dimage->nScans=1;
        dimage->inList=chChalloc(dimage->nScans,"input name list",0);
        dimage->inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        dimage->readBounds=0;
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-vRes",5)){
        checkArguments(3,i,argc,"-vRes");
        for(j=0;j<3;j++)dimage->vRes[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->vRes[0]=dimage->vRes[1]=dimage->vRes[2]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rRes",5)){
        checkArguments(1,i,argc,"-rRes");
        dimage->rRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-coord",6)){
        checkArguments(3,i,argc,"coord");
        for(j=0;j<3;j++)dimage->origin[j]=atof(argv[++i]);
        dimage->doVoxCov=0;              /*do a single beam*/
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gauss",6)){
        dimage->gaussFoot=1;
      }else if(!strncasecmp(argv[i],"-vect",5)){
        checkArguments(3,i,argc,"-vect");
        for(j=0;j<3;j++)dimage->grad[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-voxGap",7)){
        dimage->silouhetteVox=0;
      }else if(!strncasecmp(argv[i],"-readBounds",11)){
        dimage->readBounds=1;
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        dimage->maxZen=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-ang",4)){
        checkArguments(2,i,argc,"-ang");
        zen=atof(argv[++i])*M_PI/180.0;
        az=atof(argv[++i])*M_PI/180.0;
        dimage->grad[0]=sin(zen)*sin(az);
        dimage->grad[1]=sin(zen)*cos(az);
        dimage->grad[2]=-1.0*cos(zen);
      }else if(!strncasecmp(argv[i],"-maxRefl",8)){
        checkArguments(1,i,argc,"-maxRefl");
        dimage->lidPar.maxRefl=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-appRefl",8)){
        checkArguments(1,i,argc,"-appRefl");
        dimage->lidPar.appRefl=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minGap",7)){
        checkArguments(1,i,argc,"-minGap");
        dimage->lidPar.minGap=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-beamDiv",8)){
        checkArguments(1,i,argc,"-beamDiv");
        tempAng=M_PI*atof(argv[++i])/180.0;
        dimage->lidPar.beamTanDiv=sin(tempAng)/cos(tempAng);
      }else if(!strncasecmp(argv[i],"-correctR",9)){
        dimage->lidPar.correctR=1;
      }else if(!strncasecmp(argv[i],"-printFilled",12)){
        dimage->onlyPrintFilled=1;
      }else if(!strncasecmp(argv[i],"-dtm",4)){
        checkArguments(1,i,argc,"-dtm");
        dimage->useDEM=1;
        strcpy(dimage->demNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-demTol",7)){
        checkArguments(1,i,argc,"-demTol");
        dimage->demTol=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hdf",4)){
        dimage->writeHDF=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n################################\nProgram to create GEDI waveforms from TLS las files\nor voxelise TLS data\n################################\n\n-input name;     input TLS filename. rxp or ptx\n-inList list;    input file list for multiple files\n-output name;    output filename\n-hdf;    write output as a compressed HDF5 file\n-bounds minX minY minZ maxX maxY maxZ:    voxel bounds\n-res res;        voxel resolution in metrs\n-vRes x y z;     voxel resolution in each dimension\n-coord x y z;    do a single beam at this origin\n-grad x y z;     beam vector\n-ang zen az;     beam angle in degrees\n-pSigma sigma;      pulse width\n-fSigma sigma;      footprint width\n-gauss;             Gaussian footprint\n-voxGap;            just do the voxel gap fractions\n-readBounds;        read bounds from data\n-maxRefl x;      DN for reflectance=1\n-beamDiv ang;    beam divergence, in degrees\n-correctR;       correct reflectance for range in gap fraction\n-minGap gap;     minimum trustable gaop fraction when scaling\n-appRefl rho;    scalar between TLS and ALS reflectance\n-printFilled;    only write out filled voxel coordinates\n\n# DTM\n-dtm name;       use a DTM to remove ground hits\n-demTol tol;     tolerance for counting points as ground from DTM\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  dimage->buffLen=4.0*dimage->pSigma+5.0;


  return(dimage);
}/*readCommands*/

/*the end*/
/*##################################################*/

