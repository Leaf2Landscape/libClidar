#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "hdf5.h"
#include "libDEMhandle.h"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidVoxel.h"
#include "tiffWrite.h"

/*#######################*/
/*# Voxelises NERC-ARSF #*/
/*# las 1.3 data        #*/
/*# S Hancock, 2014     #*/
/*#######################*/


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



#define TOL 0.0000000001
int dynRange;    /*dynamic range*/

/*REMOVE BEFORE USE*/
/*int counter;*/

/*############################################*/
/*control structure*/

typedef struct{
  /*input/output names*/
  char listNamen[200];    /*input file list*/
  char outRoot[100];
  /*options*/
  char doImage;    /*reflectance image switch*/
  char doBins;     /*individual bin output*/
  char writeGauss; /*write Gaussian parameter switch*/
  char discOut;    /*output discrete data or not*/
  char histoStats; /*work out histogram and statistics*/
  char varNoise;   /*variable noise*/
  char getGBIC;    /*gather data for GBIC table*/
  char appGBIC;    /*apply GBIC correction*/
  char balFlights; /*balance intensity between flights*/
  char agcImage;   /*AGC variable image switch*/
  char runCheck;       /*write a mark if we're running*/
  char grabDir[200];   /*directory to write mark to*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  char correctAtten;   /*correct attenuation switch*/

  /*lidar variables*/
  float res;       /*range resolution*/

  /*voxels*/
  char doVoxel;    /*voxel switch*/
  float vRes[3];   /*voxel resolution*/
  float beamRad;   /*beam radius*/
  char asciiVox;   /*ASCII output swicth*/
  char tiffVox;    /*geotiff output switch*/
  char doVoxRH;    /*geotiff of RH metric switch*/
  char voxRHbyte;  /*byte/float RH output switch*/
  char fromGr;     /*output height from ground*/

  /*canopy bounds*/
  char useCanBound;    /*switch*/
  char useCanTop;      /*top switch*/
  char canBounds;      /*determine canopy bound switch*/
  char canNamen[200];  /*bound input file*/
  canBstruct canB;     /*canopy bounds structure*/

  /*image*/
  int iNx;              /*image dimension*/
  int iNy;              /*image dimension*/
  float iRes;           /*image resolution*/
  unsigned char *image; /*image array*/
  int geoI[2];          /*geolocation coordinate*/
  double geoL[2];       /*geolocation coordinate*/
  float maxI;           /*maximum intensity*/

  /*global lidar parameters*/
  float bounds[6];      /*file bounds, minX, minY, minZ, maxX, maxY, maxZ*/
  double uBounds[6];     /*to use bounds, minX, minY, minZ, maxX, maxY, maxZ*/

  denPar denoise;      /*structure of denoising variables*/
  float nodata;        /*nodata flag*/
}control;  /*control structure*/


/*############################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  listLas *lasList=NULL;
  lasFile **lasIn=NULL;          /*las files*/
  void readHeader(lasFile *,char *);
  control *dimage=NULL;
  control *readCommands(int,char **);
  void setUpResults(control *,lasFile **,listLas *);
  void drawPicture(control *,lasFile **,listLas *);
  void canopyBounds(control *,lasFile **,listLas *);
  void makeVoxels(control *,lasFile **,listLas *);
  void writeBins(control *,lasFile **,listLas *);
  void histogramStats(control *,lasFile **,listLas *);
  void getGBICdata(control *,lasFile **,listLas *);
  void writeRunCheck(control *);
  void readCanBounds(canBstruct *,char *,double *);
  void readPulse(denPar *);


  /*defaults*/
  dimage=readCommands(argc,argv);

  /*if needed, mark that we're running*/
  if(dimage->runCheck)writeRunCheck(dimage);

  /*read input filename list*/
  lasList=readLasList(dimage->listNamen);

  /*allocate space for las files*/
  lasIn=lfalloc(lasList->nFiles);

  /*open and read file headers*/
  for(i=0;i<lasList->nFiles;i++)lasIn[i]=readLasHead(lasList->nameList[i],dimage->pBuffSize);

  /*read ALS pulse*/
  if(dimage->denoise.deconMeth>=0)readPulse(&dimage->denoise);

  /*read GBIC table*/
  readGBIC(dimage->appGBIC,dimage->balFlights,lasIn,lasList);

  /*determine results arrays*/
  setUpResults(dimage,lasIn,lasList);

  /*these should be combined for efficiency*/
  /*canopy bounds*/
  if(dimage->canBounds)       canopyBounds(dimage,lasIn,lasList);
  else if(dimage->useCanBound)readCanBounds(&dimage->canB,dimage->canNamen,dimage->uBounds);

  /*intensity image*/
  if(dimage->doImage)drawPicture(dimage,lasIn,lasList);

  /*voxel map*/
  if(dimage->doVoxel)makeVoxels(dimage,lasIn,lasList);

  /*individual bins*/
  if(dimage->doBins||dimage->writeGauss)writeBins(dimage,lasIn,lasList);

  /*histogram and stats*/
  if(dimage->histoStats)histogramStats(dimage,lasIn,lasList);

  /*GBIC data*/
  if(dimage->getGBIC)getGBICdata(dimage,lasIn,lasList);

  /*tidy all arrays*/
  if(lasIn){
    for(i=0;i<lasList->nFiles;i++){
      if(lasIn[i]->ipoo){
        fclose(lasIn[i]->ipoo);
        lasIn[i]->ipoo=NULL;
      }
      TIDY(lasIn[i]->gbic);
      TIDY(lasIn[i]);
    }
    TTIDY((void **)lasIn,lasList->nFiles);
  }
  tidyListLas(lasList);
  if(dimage){
    TIDY(dimage->image);
    TIDY(dimage->canB.canMin);
    TIDY(dimage->canB.canMax);
    TTIDY((void **)dimage->denoise.pulse,2);
    TIDY(dimage->denoise.matchPulse);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*############################################*/
/*get data for GBIC table*/

void getGBICdata(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0,j=0,k=0;
  float waveRefl=0;
  double xCent=0,yCent=0,zCent=0;
  /*uint64_t *hist=NULL;*/
  char boundsCheck(double,double,double,double *);
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.gbic",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"# 1 AGC, 2 disc, 3 wave, 4 waveLen\n");

  /*if(!(hist=(uint64_t *)calloc(257,sizeof(uint64_t)))){
    fprintf(stderr,"error contN allocation.\n");
    exit(1);
  }
  for(i=0;i<257;i++)hist[i]=0;*/

  for(i=0;i<lasList->nFiles;i++){
    fprintf(stdout,"Processing file %d of %d\n",i+1,lasList->nFiles);
    if(checkFileBounds(lasIn[i],dimage->uBounds[0],dimage->uBounds[3],dimage->uBounds[1],dimage->uBounds[4])){
      for(j=0;j<lasIn[i]->nPoints;j++){   /*point loop*/
        /*read point*/
        readLasPoint(lasIn[i],j);
        setCoords(&xCent,&yCent,&zCent,lasIn[i]);
        /*hist[lasIn[i]->psID]++;*/

        if(boundsCheck(xCent,yCent,zCent,dimage->uBounds)){
          if((lasIn[i]->packetDes)&&(lasIn[i]->field.nRet==1)){  /*read wave if this is the first point*/
            lasIn[i]->wave=readLasWave(lasIn[i]->waveMap,lasIn[i]->waveLen,lasIn[i]->ipoo,lasIn[i]->waveStart);
            waveRefl=0.0;
            for(k=0;k<lasIn[i]->waveLen;k++){   /*GBIC is applied below for blancing flight lines*/
              if((float)lasIn[i]->wave[k]>dimage->denoise.thresh)waveRefl+=((float)lasIn[i]->wave[k]-dimage->denoise.meanN)/lasIn[i]->gbic[lasIn[i]->psID];
            }/*wave loop*/
            TIDY(lasIn[i]->wave);
          }else if(lasIn[i]->packetDes==0){/*read waveform*/
            waveRefl=-1.0;
          }
          fprintf(opoo,"%d %d %f %d %f %f %f\n",lasIn[i]->psID,lasIn[i]->refl,waveRefl,lasIn[i]->waveLen,xCent,yCent,zCent);
        }/*check point is within*/
      }/*point loop*/
    }/*check file b(unds*/
  }/*file loop*/

  /*for(i=0;i<257;i++)fprintf(opoo,"%d %ld\n",i,hist[i]);
  TIDY(hist);*/


  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",namen);
  return;
}/*getGBICdata*/


/*############################################*/
/*write out grid of min/max*/

void canopyBounds(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0,j=0,k=0;                    /*loop controllers*/
  int xBin=0,yBin=0,place=0;          /*image indices*/
  double x=0,y=0,z=0;                 /*bin coordinates*/
  double xCent=0,yCent=0,zCent=0;     /*point coordinates*/
  double buffer=0;                    /*buffer for beam bounds*/
  char namen[200];
  FILE *opoo=NULL;
  /*denoise and deconvolve*/
  float *decon=NULL;

  buffer=6.0;

  for(i=0;i<lasList->nFiles;i++){ /*file loop*/
    fprintf(stdout,"Processing file %d of %d\n",i+1,lasList->nFiles);

    /*check the file bounds*/
    if(checkFileBounds(lasIn[i],dimage->uBounds[0],dimage->uBounds[3],dimage->uBounds[1],dimage->uBounds[4])){
      for(j=0;j<lasIn[i]->nPoints;j++){   /*point loop*/
        /*read point*/
        readLasPoint(lasIn[i],j);
        if(checkOneWave(lasIn[i])){

          /*find a point on the line. Last return*/
          setCoords(&xCent,&yCent,&zCent,lasIn[i]);

          /*we should check over a buffer here, say 100m?*/
          if((xCent>=(dimage->uBounds[0]-buffer))&&(yCent>=(dimage->uBounds[1]-buffer))&&\
             (zCent>=(dimage->uBounds[2]-buffer))&&(xCent<=(dimage->uBounds[3]+buffer))&&\
             (yCent<=(dimage->uBounds[4]+buffer))&&(zCent<=(dimage->uBounds[5]+buffer))){

            /*read waveform*/
            lasIn[i]->wave=readLasWave(lasIn[i]->waveMap,lasIn[i]->waveLen,lasIn[i]->ipoo,lasIn[i]->waveStart);

            /*denoise and deconvolve*/
            decon=processWave(lasIn[i]->wave,(int)lasIn[i]->waveLen,&(dimage->denoise),lasIn[i]->gbic[lasIn[i]->psID]);

            /*determine highest and lowest point*/
            for(k=0;k<(int)lasIn[i]->waveLen;k++){
              if(decon[k]>TOL){ /*brightness check*/
                binPosition(&x,&y,&z,k,xCent,yCent,zCent,lasIn[i]->time,lasIn[i]->grad);
                if((x>=dimage->uBounds[0])&&(y>=dimage->uBounds[1])&&(z>=dimage->uBounds[2])&&\
                   (x<=dimage->uBounds[3])&&(y<=dimage->uBounds[4])&&(z<=dimage->uBounds[5])){
                  xBin=(int)((x-dimage->uBounds[0])/(double)dimage->canB.cRes+0.5);
                  yBin=(int)((y-dimage->uBounds[1])/(double)dimage->canB.cRes+0.5);

                  if((xBin>=0)&&(xBin<dimage->canB.cNx)&&(yBin>=0)&&(yBin<dimage->canB.cNy)){
                    place=yBin*dimage->canB.cNx+xBin;
                    if(z<dimage->canB.canMin[place])dimage->canB.canMin[place]=z;
                    if(z>dimage->canB.canMax[place])dimage->canB.canMax[place]=z;
                  }/*image index check*/
                }/*bin bounds check*/
              }/*threshold check*/
            }/*bin loop*/
            TIDY(decon);
            TIDY(dimage->denoise.gPar);
            TIDY(lasIn[i]->wave);
          }/*point bound check*/
        }/*point has wave check*/
      }/*point loop*/
    }/*file bounds check*/
  }/*file loop*/

  dimage->canB.cUbound[0]=dimage->uBounds[0];
  dimage->canB.cUbound[1]=dimage->uBounds[1];
  dimage->canB.cUbound[2]=dimage->uBounds[3];
  dimage->canB.cUbound[3]=dimage->uBounds[4];

  /*output results*/
  sprintf(namen,"%s.canBounds",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  for(i=0;i<dimage->canB.cNx;i++){
    for(j=0;j<dimage->canB.cNy;j++){
      place=j*dimage->canB.cNx+i;
      fprintf(opoo,"%f %f %f %f\n",((float)i+0.5)*dimage->canB.cRes+dimage->uBounds[0],\
             ((float)j+0.5)*dimage->canB.cRes+dimage->uBounds[1],dimage->canB.canMin[place],dimage->canB.canMax[place]);
    }/*y loop*/
  }/*x loop*/
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Canopy bounds to %s\n",namen);

  return;
}/*canopyBounds*/


/*############################################*/
/*write out individual bins*/

void writeBins(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0,j=0,k=0;                    /*loop controllers*/
  int wStart=0,wEnd=0;
  double x=0,y=0,z=0;
  double xCent=0,yCent=0,zCent=0;
  double buffer=0;                    /*buffer for beam bounds*/
  double setCanGround(double,double,canBstruct *);
  void setCanInd(int *,int *,double,double,double,lasFile *,canBstruct *);
  void imposeCanBound(float *,int,int,int);
  char namen[200],discNamen[200],gaussNamen[200];
  char useIt=0;
  FILE *opoo=NULL;
  FILE *discFile=NULL;
  FILE *gaussPoo=NULL;
  /*denoise and deconvolve*/
  float *decon=NULL;
  float *correctAttenuation(float *,int);
  float *trueArea=NULL;



  buffer=6.0;

  /*individual bins*/
  if(dimage->doBins){
    sprintf(namen,"%s.wave",dimage->outRoot);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
    fprintf(opoo,"# 1 lon, 2 lat, 3 height, 4 intensity, 5 file numb, 6 wave numb, 7 AGC\n");
  }

  /*Gaussian parameters*/
  if(dimage->writeGauss){
    sprintf(gaussNamen,"%s.gauss",dimage->outRoot);
    if((gaussPoo=fopen(gaussNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",gaussNamen);
      exit(1);
    }
    fprintf(gaussPoo,"# 1 lon, 2 lat, 3 height, 4 amplitude, 5 width, 6 energy\n");
  }


  if(dimage->discOut){
    sprintf(discNamen,"%s.disc",dimage->outRoot);
    if((discFile=fopen(discNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",discNamen);
      exit(1);
    }
    fprintf(discFile,"# 1 lon, 2 lat, 3 height, 4 intensity, 5 file numb, 6 wave numb, 7 projected area\n");
  }/*open discrete output file*/


  for(i=0;i<lasList->nFiles;i++){ /*file loop*/
    fprintf(stdout,"Processing file %d of %d\n",i+1,lasList->nFiles);

    /*check the file bounds*/
    if(checkFileBounds(lasIn[i],dimage->uBounds[0],dimage->uBounds[3],dimage->uBounds[1],dimage->uBounds[4])){
      for(j=0;j<lasIn[i]->nPoints;j++){   /*point loop*/
        /*read point*/
        readLasPoint(lasIn[i],j);
        if(checkOneWave(lasIn[i])){

          /*find a point on the line. Last return*/
          setCoords(&xCent,&yCent,&zCent,lasIn[i]);

          /*we should check over a buffer here, say 100m?*/
          if((xCent>=(dimage->uBounds[0]-buffer))&&(yCent>=(dimage->uBounds[1]-buffer))&&\
             (zCent>=(dimage->uBounds[2]-buffer))&&(xCent<=(dimage->uBounds[3]+buffer))&&\
             (yCent<=(dimage->uBounds[4]+buffer))&&(zCent<=(dimage->uBounds[5]+buffer))){

            /*read waveform*/
            lasIn[i]->wave=readLasWave(lasIn[i]->waveMap,lasIn[i]->waveLen,lasIn[i]->ipoo,lasIn[i]->waveStart);

            /*denoise and deconvolve*/
            decon=processWave(lasIn[i]->wave,(int)lasIn[i]->waveLen,&(dimage->denoise),lasIn[i]->gbic[lasIn[i]->psID]);

            if(dimage->doBins){  /*write out individual bins*/
              /*set canopy bounds if used*/
              if(dimage->useCanBound){
                setCanInd(&wStart,&wEnd,xCent,yCent,zCent,lasIn[i],&dimage->canB);
                imposeCanBound(decon,wStart,wEnd,(int)lasIn[i]->waveLen);
              }else{
                wStart=0;
                wEnd=lasIn[i]->waveLen;
              }

              /*correct for attenuation*/
              if(dimage->correctAtten)trueArea=correctAttenuation(decon,lasIn[i]->waveLen);
              else                    trueArea=decon;

              /*write out accepted bins*/
              for(k=wStart;k<wEnd;k++){/*bin loop*/
                if((decon[k]>TOL)||(trueArea[k]>TOL)){ /*brightness check*/
                  binPosition(&x,&y,&z,k,xCent,yCent,zCent,lasIn[i]->time,lasIn[i]->grad);
                  if((x>=dimage->uBounds[0])&&(y>=dimage->uBounds[1])&&(z>=dimage->uBounds[2])&&\
                     (x<=dimage->uBounds[3])&&(y<=dimage->uBounds[4])&&(z<=dimage->uBounds[5])){
                    fprintf(opoo,"%.2f %.2f %.2f %f %d %d %d %f\n",x,y,z,decon[k],i,j,(int)lasIn[i]->psID,trueArea[k]);
                  }
                }/*threshold check*/
              }/*bin loop*/
            }/*write out individual bins*/

            if(dimage->writeGauss){  /*write out Gaussians*/
              for(k=0;k<dimage->denoise.nGauss;k++){ /*Gaussian loop*/
                x=xCent+(double)dimage->denoise.gPar[3*k]*(double)lasIn[i]->grad[0];
                y=yCent+(double)dimage->denoise.gPar[3*k]*(double)lasIn[i]->grad[1];
                z=zCent+(double)dimage->denoise.gPar[3*k]*(double)lasIn[i]->grad[2];
                if((x>=dimage->uBounds[0])&&(y>=dimage->uBounds[1])&&(z>=dimage->uBounds[2])&&\
                   (x<=dimage->uBounds[3])&&(y<=dimage->uBounds[4])&&(z<=dimage->uBounds[5])){
                  /*check canopy bounds*/
                  if(dimage->useCanBound){
                    if(z>=setCanGround(x,y,&dimage->canB))useIt=1;
                    else                                  useIt=0;
                  }else useIt=1;
                  if(useIt)fprintf(gaussPoo,"%.2f %.2f %.2f %f %f %f\n",x,y,z,dimage->denoise.gPar[3*k+1],\
                    dimage->denoise.gPar[3*k+2],sqrt(2.0*M_PI)*dimage->denoise.gPar[3*k+1]*dimage->denoise.gPar[3*k+2]);
                }/*bounds check*/
              }/*Gaussian loop*/
            }/*write out Gaussians*/

            TIDY(decon);
            TIDY(lasIn[i]->wave);
            TIDY(dimage->denoise.gPar);
            if(dimage->correctAtten){
              TIDY(trueArea);
            }else{
              trueArea=NULL;
            }
          }/*point bounds check*/
        }/*check the point has a waveform*/

        if(dimage->discOut){ /*write discrete*/
          setCoords(&xCent,&yCent,&zCent,lasIn[i]);
          if((xCent>=dimage->uBounds[0])&&(yCent>=dimage->uBounds[1])&&(zCent>=dimage->uBounds[2])&&\
             (xCent<=dimage->uBounds[3])&&(yCent<=dimage->uBounds[4])&&(zCent<=dimage->uBounds[5])){
            fprintf(discFile,"%f %f %f %d %d %d\n",xCent,yCent,zCent,lasIn[i]->refl,i,j);
          }/*bounds check*/
        }/*write discrete*/
      }/*point loop*/
    }/*file bound check*/
  }/*file loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
    fprintf(stdout,"Bins written to %s\n",namen);
  }
  if(gaussPoo){
    fclose(gaussPoo);
    gaussPoo=NULL;
    fprintf(stdout,"Gaussians written to %s\n",gaussNamen);
  }
  if(discFile){
    fclose(discFile);
    discFile=NULL;
    fprintf(stdout,"Discrete written to %s\n",discNamen);
  }
  return;
}/*writeBins*/


/*############################################*/
/*make voxels*/

void makeVoxels(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0;        /*loop controllers*/
  uint32_t j=0;
  int nRH=0;                  /*number of RH metrics to output*/
  int wStart=0,wEnd=0;        /*waveform bounds*/
  double xCent=0,yCent=0,zCent=0;     /*a point on the beam*/
  double buffer=0;                    /*buffer for beam bounds*/
  float *correctAttenuation(float *,int);
  float *decon=NULL,*trueArea=NULL;
  float **voxelRHmetrics(voxStruct *,canBstruct *,int *,control *);
  float **voxMetrics=NULL;
  voxStruct *vox=NULL;
  voxStruct *voxAllocate(int,float *,double *,char);
  voxStruct *tidyVox(voxStruct *);
  void writeTiffVox(voxStruct *,char *,canBstruct *,char);
  void setCanInd(int *,int *,double,double,double,lasFile *,canBstruct *);
  void imposeCanBound(float *,int,int,int);
  void voxelate(voxStruct *,float *,lasFile *,double,double,double,float);
  void writeAsciiVox(voxStruct *,char *);
  void writeVoxMetrics(voxStruct *,float **,int,char *,char);
  char usable=0;

  /*buffer around plot edges to check waveforms*/
  buffer=30.0;

  vox=voxAllocate(1,&(dimage->vRes[0]),dimage->uBounds,0);

  for(i=0;i<lasList->nFiles;i++){ /*file loop*/
    fprintf(stdout,"Voxelising file %d of %d\n",i+1,lasList->nFiles);
    /*check the file bounds*/
    if(checkFileBounds(lasIn[i],dimage->uBounds[0],dimage->uBounds[3],dimage->uBounds[1],dimage->uBounds[4])){
      /*point loop*/
      for(j=0;j<lasIn[i]->nPoints;j++){
        /*read point*/
        readLasPoint(lasIn[i],j);
        if((lasIn[i]->packetDes)&&(lasIn[i]->field.nRet==lasIn[i]->field.retNumb)){ /*check there is a waveform and just one per beam*/
          /*find a point on the line. Last return*/
          setCoords(&xCent,&yCent,&zCent,lasIn[i]);

          /*we should check over a buffer here, say 100m?*/
          if((xCent>=(dimage->uBounds[0]-buffer))&&(yCent>=(dimage->uBounds[1]-buffer))&&\
             (zCent>=(dimage->uBounds[2]-buffer))&&(xCent<=(dimage->uBounds[3]+buffer))&&\
             (yCent<=(dimage->uBounds[4]+buffer))&&(zCent<=(dimage->uBounds[5]+buffer))){

            /*read waveform*/
            lasIn[i]->wave=readLasWave(lasIn[i]->waveMap,lasIn[i]->waveLen,lasIn[i]->ipoo,lasIn[i]->waveStart);
            /*denoise and deconvolve*/
            decon=processWave(lasIn[i]->wave,(int)lasIn[i]->waveLen,&(dimage->denoise),lasIn[i]->gbic[lasIn[i]->psID]);
            /*impose canopy bounds*/
            if(dimage->useCanBound){
              setCanInd(&wStart,&wEnd,xCent,yCent,zCent,lasIn[i],&dimage->canB);
              imposeCanBound(decon,wStart,wEnd,(int)lasIn[i]->waveLen);
            }
            /*correct for attenuation*/
            if(dimage->correctAtten){
              trueArea=correctAttenuation(decon,lasIn[i]->waveLen);
              TIDY(decon);
            }else{
              trueArea=decon;
              decon=NULL;
            }

            /*check waveform is usable*/
            usable=checkWaveform(trueArea,lasIn[i]->waveLen);

            /*map to voxels*/
            if(usable)voxelate(vox,trueArea,lasIn[i],xCent,yCent,zCent,dimage->beamRad);
            TIDY(trueArea);
            TIDY(dimage->denoise.gPar);
          }/*bounds check*/
        }/*waveform check*/
      }/*point loop*/
    }/*file overlap check*/
  }/*file loop*/

  /*normalise*/
  for(i=vox->nX*vox->nY*vox->nZ-1;i>=0;i--){
    if(vox->contN[i]>0){
      vox->hits[0][i]/=(float)vox->contN[i];
      vox->miss[0][i]/=(float)vox->contN[i];
    }else{
      vox->hits[0][i]=vox->miss[0][i]=dimage->nodata;
    }
  }

  /*write out ascii*/
  if(dimage->asciiVox)writeAsciiVox(vox,dimage->outRoot);

  /*write out geotiffs*/
  if(dimage->tiffVox){
    if(dimage->useCanBound)writeTiffVox(vox,dimage->outRoot,&(dimage->canB),dimage->fromGr);
    else                   writeTiffVox(vox,dimage->outRoot,NULL,0);
  }

  /*calculate RH metrics*/
  if(dimage->doVoxRH){
    voxMetrics=voxelRHmetrics(vox,&dimage->canB,&nRH,dimage);
    writeVoxMetrics(vox,voxMetrics,nRH,dimage->outRoot,dimage->voxRHbyte);
    TTIDY((void **)voxMetrics,nRH);
  }
  

  vox=tidyVox(vox);
  return;
}/*makeVoxels*/


/*############################################*/
/*wrirte voxel metrics as geotiffs*/

void writeVoxMetrics(voxStruct *vox,float **voxMetrics,int nRH,char *outRoot,char voxRHbyte)
{
  int i=0,j=0;
  int geoI[2];
  float min=0,max=0;
  char **roots=NULL;
  char nameRoot[200];
  unsigned char *image=NULL;
  double geoL[2];

  /*set up name roots*/
  roots=chChalloc(nRH,"name roots",0);
  for(i=0;i<nRH;i++)roots[i]=challoc(50,"name roots",i+1);
  strcpy(roots[0],"RH25");
  strcpy(roots[1],"RH50");
  strcpy(roots[2],"RH75");
  strcpy(roots[3],"RH98");
  strcpy(roots[4],"FHD");
;
  /*geolocation information*/
  geoI[0]=geoI[1]=0;
  geoL[0]=vox->bounds[0]+0.5*vox->res[0];
  geoL[1]=vox->bounds[4]-0.5*vox->res[1];

  if(voxRHbyte)image=uchalloc((uint64_t)vox->nX*(uint64_t)vox->nY,"image",0);

  /*print out each RH metric file*/
  for(i=0;i<nRH;i++){
    if(voxRHbyte){  /*1 byte image*/
      sprintf(nameRoot,"%s.%s.tif",outRoot,roots[i]);
      /*determine min/max*/
      min=0.0;
      max=-1000000000.0;
      for(j=vox->nX*vox->nY-1;j>=0;j--){
        /*if(voxMetrics[i][j]<min)min=voxMetrics[i][j];*/
        if(voxMetrics[i][j]>max)max=voxMetrics[i][j];
      }
      /*scale image*/
      for(j=vox->nX*vox->nY-1;j>=0;j--){
        if((voxMetrics[i][j]>=min)&&(voxMetrics[i][j]))image[j]=(unsigned char)((voxMetrics[i][j]-min)*(255.0/(max-min)));
        else if(voxMetrics[i][j]<min)                 image[j]=0;
        else if(voxMetrics[i][j]>max)                 image[j]=255;
      }
      drawTiff(nameRoot,&(geoL[0]),&(geoI[0]),vox->res[0],image,vox->nX,vox->nY,(double)((max-min)/255.0),27700);
    }else{   /*floating point image*/
      sprintf(nameRoot,"%s.%s",outRoot,roots[i]);
      drawTiffFlo(nameRoot,&(geoL[0]),&(geoI[0]),vox->res[0],&(voxMetrics[i][0]),vox->nX,vox->nY,1.0,27700);
    }
  }

  TIDY(image);
  TTIDY((void **)roots,nRH);
  return;
}/*writeVoxMetrics*/


/*############################################*/
/*calculate RH metrics from voxels*/

float **voxelRHmetrics(voxStruct *vox,canBstruct *canB,int *nRH,control *dimage)
{
  int i=0,j=0,k=0;
  int vPlace=0,mPlace=0;
  int ind=0;
  float **voxMetrics=NULL;
  float zG=0,*wave=NULL;
  float total=0,cumul=0;
  float denom=0;
  double x=0,y=0;
  double setCanGround(double,double,canBstruct *);
  char *done=NULL;


  if(canB->canMin==NULL){
    fprintf(stderr,"-voxRH requires a DEM specifiying\n");
    exit(1);
  }

  wave=falloc((uint64_t)vox->nZ,"voxel waveform",0);
  *nRH=5;
  voxMetrics=fFalloc(*nRH,"voxMetrics",0);
  done=challoc((uint64_t)(*nRH),"done flag",0);
  for(i=0;i<(*nRH);i++)voxMetrics[i]=falloc((uint64_t)vox->nX*(uint64_t)vox->nY,"voxMetrics",(*nRH)+1);


  for(i=0;i<vox->nX;i++){
    x=(double)i*vox->res[0]+vox->bounds[0];
    for(j=0;j<vox->nY;j++){
      mPlace=(vox->nY-j-1)*vox->nX+i;   /*tiff starts from top right*/
      y=(double)j*vox->res[1]+vox->bounds[1];

      /*ground elevation*/
      zG=(float)setCanGround(x,y,canB)-vox->bounds[2];

      /*make voxel waveform*/
      total=0.0;
      for(k=0;k<vox->nZ;k++)wave[k]=0.0;  /*reset*/
      for(k=0;k<vox->nZ;k++){
        vPlace=k*vox->nX*vox->nY+j*vox->nX+i;
        denom=vox->hits[0][vPlace]+vox->miss[0][vPlace];
        if(denom>0.0)wave[k]=vox->hits[0][vPlace]/denom;
        total+=wave[k];
      }

      /*check that we have some information here*/
      if(total>0.0){
        /*calculate RH metrics*/
        for(k=0;k<(*nRH);k++)done[k]=0;
        cumul=0.0;
        for(k=0;k<vox->nZ;k++){
          cumul+=wave[k]/total;

          ind=0;
          if((cumul>=0.25)&&(done[ind]==0)){
            done[ind]=1;
            voxMetrics[ind][mPlace]=(float)((double)k*vox->res[2]-zG);
          }

          ind=1;
          if((cumul>=0.5)&&(done[ind]==0)){
            done[ind]=1;
            voxMetrics[ind][mPlace]=(float)((double)k*vox->res[2]-zG);
          }

          ind=2;
          if((cumul>=0.75)&&(done[ind]==0)){
            done[ind]=1;
            voxMetrics[ind][mPlace]=(float)((double)k*vox->res[2]-zG);
          }

          ind=3;
          if((cumul>=0.98)&&(done[ind]==0)){
            done[ind]=1;
            voxMetrics[ind][mPlace]=(float)((double)k*vox->res[2]-zG);
          }
        }

        /*prevent aliasing*/
        for(ind=0;ind<4;ind++)if(fabs(voxMetrics[ind][mPlace])<=vox->res[2])voxMetrics[ind][mPlace]=0.0;


        /*foliage height diversity*/
        ind=4;
        voxMetrics[ind][mPlace]=foliageHeightDiversity(wave,vox->nZ);
      }else{/*empty check*/
        for(ind=0;ind<(*nRH);ind++)voxMetrics[ind][mPlace]=dimage->nodata;
      }
    }/*y loop*/
  }/*x loop*/

  TIDY(wave);
  return(voxMetrics);
}/*voxelRHmetrics*/


/*############################################*/
/*write out geotiff of voxels*/

void writeTiffVox(voxStruct *vox,char *outRoot,canBstruct *canB,char fromGr)
{
  int i=0,j=0,k=0;
  int zBin=0,place=0;
  int iPlace=0;
  int nLayers=0;
  int geoI[2];
  float *image=NULL;
  double x=0,y=0,z=0;
  double geoL[2];
  double setCanGround(double,double,canBstruct *);
  char nameRoot[200];

  image=falloc((uint64_t)vox->nX*(uint64_t)vox->nY,"",0);
  nLayers=(int)((vox->bounds[5]-vox->bounds[2])/vox->res[2])+1;


  geoI[0]=geoI[1]=0;
  geoL[0]=vox->bounds[0]+0.5*vox->res[0];
  geoL[1]=vox->bounds[4]-0.5*vox->res[1];

  /*loop over layers*/
  for(k=0;k<nLayers;k++){
    sprintf(nameRoot,"%s.%d",outRoot,k);
    /*zero image*/
    for(i=vox->nX*vox->nY-1;i>=0;i--)image[i]=-0.1;

    /*build up image*/
    for(i=0;i<vox->nX;i++){
      x=((double)i+0.5)*vox->res[0]+vox->bounds[0];
      for(j=0;j<vox->nY;j++){
        y=((double)j+0.5)*vox->res[1]+vox->bounds[1];
        /*bin above ground to use*/
        if((fromGr)&&(canB)){  /*this is not working yet*/
          z=setCanGround(x,y,canB);
          if(z>=0.0)zBin=k+(int)((z-vox->bounds[2])/vox->res[2]+0.5);
          else      zBin=k;
        }else zBin=k;

        iPlace=(vox->nY-j-1)*vox->nX+i;   /*tiff starts from top row*/
        if((zBin>=0)&&(zBin<vox->nZ)){
          place=zBin*vox->nX*vox->nY+j*vox->nX+i;
          if(vox->contN[place]>0)image[iPlace]=vox->hits[0][place]/(vox->hits[0][place]+vox->miss[0][place]);
          else                   image[iPlace]=-0.1;
        }else image[iPlace]=-0.1;
      }/*y*/
    }/*x*/

    /*draw the image*/
    drawTiffFlo(nameRoot,&(geoL[0]),&(geoI[0]),vox->res[0],image,vox->nX,vox->nY,1.0,27700);
  }/*layer loop*/

  TIDY(image);
  return;
}/*writeTiffVox*/


/*############################################*/
/*draw intensity image*/

void drawPicture(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0,k=0;
  int x=0,y=0,place=0;   /*image coordinates*/
  int *tempAGC=NULL;
  uint32_t j=0;
  uint64_t *contN=NULL;
  float *temp=NULL;
  float max=0;
  double lat=0,lon=0,height=0;
  unsigned char *agcImage=NULL;
  char namen[200];


  temp=falloc((uint64_t)dimage->iNx*(uint64_t)dimage->iNy,"temp image",0);
  if(dimage->agcImage)tempAGC=ialloc(dimage->iNx*dimage->iNy,"temp AGC image",0);
  if(!(contN=(uint64_t *)calloc(dimage->iNx*dimage->iNy,sizeof(uint64_t)))){
    fprintf(stderr,"error contN allocation.\n");
    exit(1);
  }
  for(i=dimage->iNx*dimage->iNy-1;i>=0;i--){
    temp[i]=0.0;
    contN[i]=0;
    if(dimage->agcImage)tempAGC[i]=0;
  }

  dimage->geoI[0]=-1;
  for(i=0;i<lasList->nFiles;i++){ /*file loop*/
    fprintf(stdout,"File %d of %d\n",i+1,lasList->nFiles);
    for(j=0;j<lasIn[i]->nPoints;j++){  /*point loop*/
      readLasPoint(lasIn[i],j);
      if((lasIn[i]->packetDes)&&(lasIn[i]->field.nRet==lasIn[i]->field.retNumb)){ /*check there is a waveform and just one per beam*/

        /*determine point on image*/
        setCoords(&lon,&lat,&height,lasIn[i]);
        x=(int)((lon-dimage->bounds[0])/dimage->iRes);
        y=(int)((lat-dimage->bounds[1])/dimage->iRes);

        /*set geolocation*/
        if(dimage->geoI[0]<0){
          dimage->geoI[0]=x;
          dimage->geoI[1]=dimage->iNy-y-1;
          dimage->geoL[0]=lon;
          dimage->geoL[1]=lat;
        }

        if((x>=0)&&(x<dimage->iNx)&&(y>=0)&&(y<dimage->iNy)){
          place=(dimage->iNy-y-1)*dimage->iNx+x;
          lasIn[i]->wave=readLasWave(lasIn[i]->waveMap,lasIn[i]->waveLen,lasIn[i]->ipoo,lasIn[i]->waveStart);
          if(lasIn[i]->wave){
            for(k=0;k<lasIn[i]->waveLen;k++){
              if((float)lasIn[i]->wave[k]>dimage->denoise.thresh){
                temp[place]+=((float)lasIn[i]->wave[k]-dimage->denoise.meanN)/lasIn[i]->gbic[lasIn[i]->psID];
              }
            }   /*check positive so we don't wrap around in very dark signals, like roads*/
            if(dimage->agcImage)tempAGC[place]+=(int)lasIn[i]->psID;
            contN[place]++;
          }
        }
        TIDY(lasIn[i]->wave);
      }/*check there is a waveform and just one per beam*/
    }/*point loop*/
  }/*file loop*/

  /*normalise*/
  max=-10000.0;
  for(i=dimage->iNx*dimage->iNy-1;i>=0;i--){
    if(contN[i]>0){
      temp[i]/=(float)contN[i];
      if(temp[i]>max)max=temp[i];
    }
  }
  //max/=1.333; //4.0;
  if(max>dimage->maxI)max=dimage->maxI;   /*to avoid glint*/

  /*set dynamic range*/
  if(dimage->agcImage)agcImage=uchalloc((uint64_t)dimage->iNx*(uint64_t)dimage->iNy,"AGC image",0);
  for(i=dimage->iNx*dimage->iNy-1;i>=0;i--){
    if(temp[i]<=max)dimage->image[i]=(unsigned char)(temp[i]*255.0/max);
    else            dimage->image[i]=255;
    if(dimage->agcImage)if(contN[i]>0)agcImage[i]=(unsigned char)((int64_t)tempAGC[i]/contN[i]);
  }
  TIDY(temp);
  TIDY(contN);
  TIDY(tempAGC);

  /*now we must write the image*/
  sprintf(namen,"%s.tif",dimage->outRoot);
  drawTiff(namen,dimage->geoL,dimage->geoI,dimage->iRes,dimage->image,dimage->iNx,dimage->iNy,1.0,27700);

  if(dimage->agcImage){ /*write AGC image*/
    sprintf(namen,"%s.AGC.tif",dimage->outRoot);
    drawTiff(namen,dimage->geoL,dimage->geoI,dimage->iRes,agcImage,dimage->iNx,dimage->iNy,1.0,27700);
    TIDY(agcImage);
  }/*write AGC image*/

  return;
}/*drawPicture*/


/*############################################*/
/*set up results arrays*/

void setUpResults(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0,j=0;

  /*determine bounds*/
  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=100000000.0;
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=-100000000.0;
  for(i=0;i<lasList->nFiles;i++){
    for(j=0;j<3;j++){
      if(lasIn[i]->minB[j]<dimage->bounds[j])dimage->bounds[j]=lasIn[i]->minB[j];
      if(lasIn[i]->maxB[j]>dimage->bounds[j+3])dimage->bounds[j+3]=lasIn[i]->maxB[j];
    }
  }/*bound determining*/

  /*set up image*/
  if(dimage->doImage){
    dimage->iNx=(int)((dimage->bounds[3]-dimage->bounds[0])/dimage->iRes+0.5);
    dimage->iNy=(int)((dimage->bounds[4]-dimage->bounds[1])/dimage->iRes+0.5);
    dimage->image=uchalloc((uint64_t)dimage->iNx*(uint64_t)dimage->iNy,"image",0);
  }/*image set up*/

  /*bounds left blank, use full extent*/
  if(dimage->uBounds[0]<0.0)for(i=0;i<6;i++)dimage->uBounds[i]=dimage->bounds[i];

  /*set up canopy bounds*/
  if(dimage->canBounds){
    dimage->canB.cNx=(int)((dimage->uBounds[3]-dimage->uBounds[0])/dimage->canB.cRes)+1;
    dimage->canB.cNy=(int)((dimage->uBounds[4]-dimage->uBounds[1])/dimage->canB.cRes)+1;
    dimage->canB.canMax=falloc((uint64_t)dimage->canB.cNx*(uint64_t)dimage->canB.cNy,"canopy max bound",0);
    dimage->canB.canMin=falloc((uint64_t)dimage->canB.cNx*(uint64_t)dimage->canB.cNy,"canopy min bound",0);
    for(i=dimage->canB.cNx*dimage->canB.cNy-1;i>=0;i--){
      dimage->canB.canMax[i]=-100000.0;
      dimage->canB.canMin[i]=100000.0;
    }
  }/*canopy bounds*/

  return;
}/*setUpResults*/


/*############################################*/
/*calculate histogram and statistics*/

void histogramStats(control *dimage,lasFile **lasIn,listLas *lasList)
{
  int i=0,j=0;
  int *hist=NULL;
  int *waveHist(unsigned char *,uint32_t);
  float waveMean(unsigned char *,uint32_t);
  float mean=0;
  double xCent=0,yCent=0,zCent=0;
  double buffer=0;
  unsigned char getMode(int *);
  unsigned char modal=0;
  unsigned char findMedianWave(unsigned char*,uint32_t);
  unsigned char medWave=0;
  unsigned char modQuart=0; /*modal top quartile*/
  void modalDeviation(unsigned char,unsigned char *,uint32_t,unsigned char *);
  FILE *opoo=NULL;
  char namen[200];

  buffer=10.0;
  sprintf(namen,"%s.histStats",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"# 1 modal, 2 mean, 3 median, 4 modal quartile, 5 file, 6 wave\n");

  for(i=0;i<lasList->nFiles;i++){ /*file loop*/
    fprintf(stdout,"Processing file %d of %d\n",i+1,lasList->nFiles);
    /*check the file bounds*/
    if(checkFileBounds(lasIn[i],dimage->uBounds[0],dimage->uBounds[3],dimage->uBounds[1],dimage->uBounds[4])){
      for(j=0;j<lasIn[i]->nPoints;j++){   /*point loop*/
        /*read point*/
        readLasPoint(lasIn[i],j);
        if(checkOneWave(lasIn[i])){
          /*find a point on the line. Last return*/
          setCoords(&xCent,&yCent,&zCent,lasIn[i]);

          /*we should check over a buffer here, say 100m?*/
          if((xCent>=(dimage->uBounds[0]-buffer))&&(yCent>=(dimage->uBounds[1]-buffer))&&\
             (zCent>=(dimage->uBounds[2]-buffer))&&(xCent<=(dimage->uBounds[3]+buffer))&&\
             (yCent<=(dimage->uBounds[4]+buffer))&&(zCent<=(dimage->uBounds[5]+buffer))){

            /*read waveform*/
            lasIn[i]->wave=readLasWave(lasIn[i]->waveMap,lasIn[i]->waveLen,lasIn[i]->ipoo,lasIn[i]->waveStart);

            /*calculate histogram*/
            hist=waveHist(lasIn[i]->wave,lasIn[i]->waveLen);

            /*calculate mean, median and modal value*/
            modal=getMode(hist);
            mean=waveMean(lasIn[i]->wave,lasIn[i]->waveLen);
            medWave=findMedianWave(lasIn[i]->wave,lasIn[i]->waveLen);

            /*mode and top quartile of the deviation from the mode*/
            modalDeviation(modal,lasIn[i]->wave,lasIn[i]->waveLen,&modQuart);

            fprintf(opoo,"%d %f %d %d %d %d\n",modal,mean,medWave,modQuart,i,j);
            TIDY(hist);
          }/*point bounds check*/
        }/*waveform check*/
      }/*point loop*/
    }/*file bound check*/
  }/*file loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Stats written to %s\n",namen);
  return;
}/*histogramStats*/


/*############################################*/
/*find waveform median*/

unsigned char findMedianWave(unsigned char *wave,uint32_t waveLen)
{
  unsigned char median=0;
  unsigned char *temp=NULL;
  int compUchar(const void *x,const void *y);  /*function needed by qsort()*/

  temp=uchalloc((uint64_t)waveLen,"temp median",0);
  memcpy(temp,wave,waveLen);

  qsort(temp,waveLen,sizeof(unsigned char),compUchar);  /*put the contents of temp in order*/
  median=temp[waveLen/2];

  TIDY(temp);
  return(median);
}/*findMedianWave*/


/*############################################*/
/*find waveform mean*/

float waveMean(unsigned char *wave,uint32_t waveLen)
{
  uint32_t i=0;
  float mean=0;

  mean=0.0;
  for(i=0;i<waveLen;i++)mean+=(float)wave[i];
  mean/=(float)waveLen;

  return(mean);
}/*waveMean*/


/*############################################*/
/*find mode from histogram*/

unsigned char getMode(int *hist)
{
  int i=0,max=0;
  unsigned char modal=0;

  for(i=0;i<dynRange;i++){
    if(hist[i]>max){
      max=hist[i];
      modal=(unsigned char)i;
    }
  }
  return(modal);
}/*getMode*/


/*############################################*/
/*wave histogram*/

int *waveHist(unsigned char *wave,uint32_t waveLen)
{
  int *hist=NULL;
  uint32_t i=0;

  hist=ialloc(dynRange,"wave histogram",0);
  for(i=0;i<dynRange;i++)hist[i]=0;
  for(i=0;i<waveLen;i++)hist[wave[i]]++;

  return(hist);
}/*waveHist*/


/*##########################################*/
/*write a run check file*/

void writeRunCheck(control *dimage)
{
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s/%s.check",dimage->grabDir,dimage->outRoot);

  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening rucn check file %s\n",namen);
    exit(1);
  }

  fprintf(opoo,"Running in this file\n");

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  return;
}/*writeRunCheck*/


/*############################################*/
/*read the command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  dynRange=256;    /*dynamic range*/

  /*defaults*/
  dimage->doImage=0;
  dimage->iRes=2.0;
  dimage->image=NULL;
  strcpy(dimage->outRoot,"waveTeast");
  strcpy(dimage->listNamen,"/home/sh563/data/bess/lists/las1.3.BD.dat");
  dimage->doVoxel=0;
  dimage->asciiVox=0;
  dimage->tiffVox=0;
  dimage->doVoxRH=0;
  dimage->voxRHbyte=1;
  dimage->fromGr=1;
  dimage->vRes[0]=dimage->vRes[1]=dimage->vRes[2]=2.0;
  dimage->uBounds[0]=dimage->uBounds[1]=dimage->uBounds[2]=-10000000000.0;
  dimage->uBounds[3]=dimage->uBounds[4]=dimage->uBounds[5]=10000000000.0;
  dimage->doBins=0;
  dimage->writeGauss=0;
  dimage->res=0.15;
  dimage->discOut=0;
  dimage->histoStats=0;
  dimage->varNoise=0;
  dimage->getGBIC=0;
  dimage->appGBIC=0;
  dimage->balFlights=0;
  dimage->agcImage=0;
  dimage->maxI=1250.0;
  dimage->beamRad=0.165;   /*33 cm*/
  dimage->runCheck=0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->correctAtten=1;   /*correct attenuation by default*/
  dimage->nodata=-9999.0;

  /*canopy bounds*/
  dimage->canBounds=0;        /*don't calculate canopy bounds*/
  dimage->canB.cRes=10.0;
  dimage->canB.canMax=NULL;
  dimage->canB.canMin=NULL;
  dimage->canB.cNx=dimage->canB.cNy=0;
  dimage->useCanBound=0;
  strcpy(dimage->canNamen,"/home/sh563/data/bess/analysis/can_bounds/canopyBounds.filt");
  dimage->useCanTop=0;


  /*denoising*/
  setDenoiseDefault(&dimage->denoise);
  dimage->denoise.meanN=13.0;
  dimage->denoise.thresh=16.0;
  dimage->denoise.deconGauss=0;  /*for this, never use a Gaussian pulse*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->listNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-iRes",5)){
        checkArguments(1,i,argc,"-iRes");
        dimage->doImage=1;
        dimage->iRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-iMax",5)){
        checkArguments(1,i,argc,"-iMax");
        dimage->maxI=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-vRes",5)){
        checkArguments(3,i,argc,"-vRes");
        for(j=0;j<3;j++)dimage->vRes[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-voxAscii",9)){
        dimage->doVoxel=1;
        dimage->asciiVox=1;
      }else if(!strncasecmp(argv[i],"-voxTiff",8)){
        dimage->doVoxel=1;
        dimage->tiffVox=1;
      }else if(!strncasecmp(argv[i],"-voxRH",6)){
        dimage->doVoxel=1;
        dimage->doVoxRH=1;
      }else if(!strncasecmp(argv[i],"-voxFloatRH",11)){
        dimage->doVoxel=1;
        dimage->doVoxRH=1;
        dimage->voxRHbyte=0;
      }else if(!strncasecmp(argv[i],"-voxAbs",7)){
        dimage->fromGr=0;
      }else if(!strncasecmp(argv[i],"-vBounds",8)){
        checkArguments(4,i,argc,"-vBounds");
        for(j=0;j<6;j++)dimage->uBounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-doBins",7)){
        dimage->doBins=1;
      }else if(!strncasecmp(argv[i],"-writeGauss",11)){
        dimage->writeGauss=1;
      }else if(!strncasecmp(argv[i],"-thresh",7)){
        checkArguments(1,i,argc,"-thresh");
        dimage->denoise.thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-meanN",6)){
        checkArguments(1,i,argc,"-meanN");
        dimage->denoise.meanN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-tailThresh",11)){
        checkArguments(1,i,argc,"-tailThresh");
        dimage->denoise.tailThresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minW",5)){
        checkArguments(1,i,argc,"-minW");
        dimage->denoise.minWidth=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gold",5)){  /*Gold's method*/
        dimage->denoise.deconMeth=0;
      }else if(!strncasecmp(argv[i],"-RL",3)){    /*Richardson-Lucy*/
        dimage->denoise.deconMeth=1;
      }else if(!strncasecmp(argv[i],"-pScale",7)){    /*deconvolution pulse scale*/
        checkArguments(1,i,argc,"-pScale");
        dimage->denoise.pScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-smooth",7)){
        checkArguments(1,i,argc,"-smooth");
        dimage->denoise.sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-median",7)){
        checkArguments(1,i,argc,"-median");
        dimage->denoise.medLen=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-preSmooth",10)){
        checkArguments(1,i,argc,"-preSmooth");
        dimage->denoise.psWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-deThresh",9)){    /*deconvolution threshold*/
        checkArguments(1,i,argc,"-deThresh");
        dimage->denoise.deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pFile",6)){    /*deconvolution threshold*/
        checkArguments(1,i,argc,"-pFile");
        strcpy(dimage->denoise.pNamen,argv[++i]);
        dimage->denoise.deconGauss=0;                  /*do not use a Gaussian pulse*/
      }else if(!strncasecmp(argv[i],"-discOut",8)){    /*output discrete return*/
        dimage->discOut=1;
        dimage->doBins=1;
      }else if(!strncasecmp(argv[i],"-hist",5)){    /*output histogram statistics*/
        dimage->histoStats=1;
      }else if(!strncasecmp(argv[i],"-varNoise",9)){    /*variable noise switch*/
        dimage->denoise.varNoise=1;
        dimage->denoise.medStats=1;
      }else if(!strncasecmp(argv[i],"-varScale",9)){
        checkArguments(1,i,argc,"-varScale");
        dimage->denoise.threshScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noTrack",8)){
        dimage->denoise.noiseTrack=0;
      }else if(!strncasecmp(argv[i],"-getGBIC",8)){    /*gather data to calculate GBIC*/
        dimage->getGBIC=1;
      }else if(!strncasecmp(argv[i],"-GBIC",5)){    /*apply GBIC*/
        dimage->appGBIC=1;
      }else if(!strncasecmp(argv[i],"-balFlights",11)){    /*balance intensity between flights*/
        dimage->balFlights=1;
      }else if(!strncasecmp(argv[i],"-agcImage",9)){    /*balance intensity between flights*/
        dimage->agcImage=1;
      }else if(!strncasecmp(argv[i],"-fitGauss",9)){    /*fit Gaussians*/
        dimage->denoise.fitGauss=1;
      }else if(!strncasecmp(argv[i],"-gWidth",7)){    /*gather data to calculate GBIC*/
        checkArguments(1,i,argc,"-gWidth");
        dimage->denoise.gWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-canBound",9)){    /*get canopy bounds*/
        dimage->canBounds=1;
      }else if(!strncasecmp(argv[i],"-cRes",5)){    /*canopy bounds resolution*/
        checkArguments(1,i,argc,"-cRes");
        dimage->canB.cRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-useGround",10)){    /*use canopy bottom bound*/
        dimage->useCanBound=1;
      }else if(!strncasecmp(argv[i],"-useCanTop",10)){    /*use canopy top bound*/
        dimage->useCanTop=1;
      }else if(!strncasecmp(argv[i],"-canNamen",9)){    /*gather data to calculate GBIC*/
        checkArguments(1,i,argc,"-canNamen");
        strcpy(dimage->canNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-gaussPulse",11)){
        dimage->denoise.gaussPulse=1;
      }else if(!strncasecmp(argv[i],"-gaussHard",10)){
        dimage->denoise.gaussFilt=1;
      }else if(!strncasecmp(argv[i],"-findHard",9)){
        dimage->denoise.matchHard=1;
      }else if(!strncasecmp(argv[i],"-hardTol",8)){
        checkArguments(1,i,argc,"-hardTol");
        dimage->denoise.hardTol=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-matchHard",10)){
        dimage->denoise.matchHard=1;
      }else if(!strncasecmp(argv[i],"-runCheck",9)){
        checkArguments(1,i,argc,"-runCheck");
        dimage->runCheck=1;
        strcpy(dimage->grabDir,argv[++i]);
      }else if(!strncasecmp(argv[i],"-leaveAtten",11)){
        dimage->correctAtten=0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;     list of input lasfiles\n-outRoot root;   output filename root\n\nVoxels\n-voxAscii;       write ascii voxels\n-voxTiff;        write voxel geotiff\n-voxRH;          make geotiffs of voxel RH metrics. 1 byte output\n-voxFloatRH;     make geotiffs of voxel RH metrics. Float output\n-vRes x y z;     voxel resolution in axes\n-voxAbs;         output absolute height rather than from ground\n\n-doBins;         output individual bins above threshold\n-vBounds minX minY minZ maxX maxY maxZ;     limit voxel bounds\n\nImage\n-iRes res;       make an image of this resolution\n-iMax max;       maximum image intensity\n-runCheck dir;   write a run check file to a directory\n\nDenoising\n-thresh x;       waveform threshold\n-meanN n;        mean noise level\n-tailThresh t;   trailing edge threshold\n-minWidth w;     minimum feature width to accept\n-preSmooth s;    pre-smoothing width\n-smooth s;       smoothing width\n-median n;       median filter length\n-varNoise;       use variable noise thresholds\n-varScale n;     factor to scale variable noise by\n-noTrack;        don't use noise tracking\n-leaveAtten;     do not correct for attenuation\n\nDeconvolution\n-gold;           deconvolve with Gold's method\n-RL;             deconvolve with Richardson-Lucy\n-deThresh t;     deconvolution thresold\n-pScale s;       scale decon pulse length\n-pFile name;     pulse filename\n\nSwitches\n-discOut;        output discrete point cloud\n-hist;           output histogram and statistics\n-getGBIC;        write out data for GBIC table, needs vBounds\n-GBIC;           apply GBIC\n-balFlights;     balance intensity between flights, read from input list\n-agcImage;       make an image of the AGC. Needs -iRes\n\nGaussian\n-fitGauss;       fit Gaussians\n-gWidth sig;     smooth before setting initial Gaussian estimates\n-gaussPulse;     convert pulse to a Gaussian before decon\n-writeGauss;     write Gaussian output\n\nCanopy bounds\n-canBounds;      output canopy bounds\n-cRes res;       canopy bound resolution\n-useGround;      use ground bound\n-useCanTop;      use canopy top bound\n-canNamen n;     canopy bound filename\n\n-matchHard;      correlation hard target finding\n-gaussHard;      Gaussian fitting to identify hard targets\n-hardTol tol;    tolerance to scale width by\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }
  if(dimage->denoise.tailThresh<0)dimage->denoise.tailThresh=dimage->denoise.thresh;
  return(dimage);
}/*readCommands*/

/*the end*/
/*############################################*/

