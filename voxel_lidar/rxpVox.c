#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLidVoxel.h"
#include "libLasProcess.h"
#include "libOptWaveform.h"

#define TLSVOX


/*#####################*/
/*# Reads binary dump #*/
/*# of Riegl rxp and  #*/
/*# voxelises         #*/
/*# S Hancock 03 2015 #*/
/*#####################*/

/*#######################################*/
/*# Copyright 2014-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of the lidar voxel library, compareWaves.          #*/
/*#                                                                      #*/
/*# compareWaves is free software: you can redistribute it and/or modify #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# compareWaves is distributed in the hope that it will be useful,      #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#  You should have received a copy of the GNU General Public License   #*/
/*#  along with compareWaves.  If not, see <http://www.gnu.org/licenses/>#*/
/*########################################################################*/



/*####################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  char alsNamen[200];   /*list of ALS files*/
  char tlsNamen[200];   /*list of TLS files*/
  char outRoot[200];    /*output filename root*/
  char pNamen[200];     /*pulse file*/

  /*switches*/
  char waveOut;       /*output waveforms*/
  char appGBIC;       /*apply GBIC correction*/
  char balFlights;    /*balance intensity between flights*/
  char waveBitsOut;   /*flag to output bits to make TLS waveform*/
  char compBeams;     /*compare individual beams*/
  char compVox;       /*compare ALS ad TLS voxels*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/

  /*TLS wave bits*/
  int nTLSwaves;        /*number of TLS waves per beam*/
  int binOff;           /*number of prebins to look for geolocation*/
  double beamRad;       /*beam radius*/
  uint64_t counter;     /*output number of waves output*/
  char pointOut;        /*output individual points*/
  lidVoxPar tlsPar;     /*tls voxelising parameters*/

  /*voxel area*/
  double bounds[6];     /*area of interest, minX minY minZ, maxX, maxY maxZ*/
  float vRes[3];        /*voxel resolution*/

  /*error stat controls*/
  float hErrRes;        /*height error resolution*/
  float gRes;           /*cumulative gap resolution*/
  char errOut[200];     /*ouput filename*/

  /*denoising parameters*/
  denPar denoise;      /*structure of denoising variables*/
  optControl ofimage;  /*optimisation control structure*/

  /*canopy bounds*/
  char canNamen[200]; /*ground file*/
  char useCanBound;   /*switch*/
  canBstruct canB;    /*bound structure*/
}control;


/*####################################*/
/*important las point data*/

typedef struct{
  int32_t x;           /*index to calculate coordinate*/
  int32_t y;           /*index to calculate coordinate*/
  int32_t z;           /*index to calculate coordinate*/
  uint16_t refl;       /*point intensity*/
  b_field field;        /*bit field. Defined in lasHead.h*/
  unsigned char scanAng;  /*scan angle*/
  uint16_t psID;          /*point source ID, used for Leica's AGC*/
  unsigned char packetDes;/*waveform packed description*/
  uint64_t waveMap;       /*pointer to waveform in file*/
  uint32_t waveLen;       /*length of waveform in bins*/
  float time;             /*time in picoseconds of this wave*/
  float grad[3];          /*waveform gradient: Vector*/
  unsigned char *wave;    /*waveform*/
}lasPoint;


/*####################################*/
/*voxel gap structure*/

typedef struct{
  /*the voxel*/
  voxStruct *vox;   /*TLS voxels. Gap within voxel*/
  voxStruct *toTLS; /*gap to TLS voxels*/
  float **meanGap;  /*mean minimum gap for voxels*/
  float **meanRefl; /*mean reflectance for voxels*/
  int **contN;

  /*TLS map*/
  int **mapFile;        /*file per voxel*/
  uint32_t **mapPoint;  /*point per voxel*/
  int *nIn;             /*number of points per voxel*/

  /*TLS data for optimisation*/
  overStruct *tlsDat;
}tlsVoxStruct;
  

/*####################################*/
/*ALS data*/

typedef struct{
  lasPoint **point;    /*AlS point arrays*/
  lasFile **lasIn;     /*las files*/
  uint32_t *nPoints;   /*number of point records used per file*/
  int nFiles;          /*number of las files*/
}alsData;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  listLas *alsList=NULL,*tlsList=NULL;
  alsData *als=NULL;
  alsData *readALSwithin(listLas *,control *);
  tlsScan *tls=NULL;
  tlsScan *readTLSwithinn(listLas *,control *,tlsVoxStruct *);
  tlsVoxStruct *vox=NULL;
  tlsVoxStruct *tlsVoxAlloc(int,control *);
  tlsVoxStruct *tidyTLSvox(tlsVoxStruct *,int);
  void tidyALS(alsData *);
  void intersectVoxels(tlsScan **,listLas *,listLas *,control *);
  void matchWaves(alsData *,tlsScan *,int,tlsVoxStruct *,control *);
  void readCanBounds(canBstruct *,char *,double *);
  void readPulse(denPar *);


  /*read commands*/
  dimage=readCommands(argc,argv);

  /*read input filename list*/
  alsList=readLasList(dimage->alsNamen);

  /*if required read ALS pulse*/
  if(dimage->denoise.deconMeth>=0)readPulse(&dimage->denoise);

  /*if required read canopy bounds*/
  if(dimage->useCanBound)readCanBounds(&dimage->canB,dimage->canNamen,dimage->bounds);

  /*store relevent ALS points*/
  als=readALSwithin(alsList,dimage);

  /*read TLS input list*/
  tlsList=readLasList(dimage->tlsNamen);

  /*allocate voxel structure*/
  vox=tlsVoxAlloc(tlsList->nFiles,dimage);

  /*read TLS data*/
  tls=readTLSwithinn(tlsList,dimage,vox);
 
  /*now match ALS to TLS and create waveforms*/
  matchWaves(als,tls,tlsList->nFiles,vox,dimage);

  /*tidy up arrays*/
  tls=tidyTLScans(tls,tlsList->nFiles);
  vox=tidyTLSvox(vox,tlsList->nFiles);
  tidyListLas(alsList);
  tidyListLas(tlsList);
  tidyALS(als);
  if(dimage){
    TTIDY((void **)dimage->denoise.pulse,2);
    TIDY(dimage->canB.canMin);
    TIDY(dimage->canB.canMax);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################*/
/*match TLS to ALS and make waveforms*/

void matchWaves(alsData *als,tlsScan *tls,int nTLS,tlsVoxStruct *vox,control *dimage)
{
  int i=0;
  int nVox=0,*voxList=NULL;
  int wStart=0,wEnd=0;
  uint32_t j=0;
  float *denoised=NULL;
  float **tlsWave=NULL;
  float **tlsWithinBeam(int,int *,lasFile *,tlsScan *,int,double,double,double,control *,tlsVoxStruct *,double,double,double,float *denoised);
  float tlsRMSE=0,*minGap=NULL,*appRefl=NULL;   /*variables to calculate per beam and pass to voxel building*/
  double x0=0,y0=0,z0=0;
  double xCent=0,yCent=0,zCent=0;
  double discRange=0;                 /*discrete return range*/
  void writeALSpoint(lasPoint *,lasFile *);
  void outputWaves(unsigned char *,float **,int32_t,char *,uint64_t,control *,double,\
                   double,double,lasFile *,float *,int,double,uint16_t,float);
  void compareBeams(tlsVoxStruct *,unsigned char *,float *,int,double,double,double,lasFile *,control *,float *,float *,float *);
  overStruct *tidyOverStruct(overStruct *);
  overStruct *setupOverStruct(control *);
  voxStruct *alsVox=NULL;
  void setupTLSstruct(overStruct *);
  void tidyDataStruct(overStruct *);
  void buildVoxel(voxStruct *,float *,double,double,double,lasFile *,float,float,float *,float *,tlsVoxStruct *);
  void compareVox(voxStruct *,voxStruct *,control *,tlsVoxStruct *);
  void setCanInd(int *,int *,double,double,double,lasFile *,canBstruct *);
  void imposeCanBound(float *,int,int,int);
  void setTopVoxBlank(voxStruct *);
  float *correctAttenuation(float *,int);
  float *area=NULL;


  /*set up data structure*/
  if(dimage->compBeams){
    vox->tlsDat=setupOverStruct(dimage);  
    minGap=falloc(vox->vox->nScans,"minGap",0);
    appRefl=falloc(vox->vox->nScans,"appRefl",0);
  }
  if(dimage->compVox)alsVox=voxAllocate(1,&(dimage->vRes[0]),dimage->bounds,0);
  


  /*ALS file loop*/
  for(i=0;i<als->nFiles;i++){
    for(j=0;j<als->nPoints[i];j++){  /*ALS point loop*/
      als->point[i][j].wave=readLasWave(als->point[i][j].waveMap,als->point[i][j].waveLen,\
                                                als->lasIn[i]->ipoo,als->lasIn[i]->waveStart);

      /*denoise and deconvolve*/
      denoised=processWave(als->point[i][j].wave,(int)als->point[i][j].waveLen,&(dimage->denoise),1.0);

      /*work out which voxels it intersects*/
      writeALSpoint(&(als->point[i][j]),als->lasIn[i]);
      setCoords(&xCent,&yCent,&zCent,als->lasIn[i]);
      binPosition(&x0,&y0,&z0,0,xCent,yCent,zCent,als->lasIn[i]->time,als->lasIn[i]->grad);
      discRange=sqrt((x0-xCent)*(x0-xCent)+(y0-yCent)*(y0-yCent)+(z0-zCent)*(z0-zCent));
      voxList=beamVoxels(als->point[i][j].grad,x0,y0,z0,&(vox->vox->bounds[0]),&(vox->vox->res[0])\
                                   ,vox->vox->nX,vox->vox->nY,vox->vox->nZ,&nVox,dimage->beamRad,NULL,-1.0);
      /*impose canopy bounds*/
      if(dimage->useCanBound){ 
        setCanInd(&wStart,&wEnd,xCent,yCent,zCent,als->lasIn[i],&dimage->canB);
        imposeCanBound(denoised,wStart,wEnd,(int)als->point[i][j].waveLen);
      }
      /*correct attenuation*/
      area=correctAttenuation(denoised,(int)als->point[i][j].waveLen);


      if(nVox>0){
        /*set up arrays to hold data*/
        if(dimage->compBeams)setupTLSstruct(vox->tlsDat);

        /*calculate TLS waves*/
        als->lasIn[i]->wave=als->point[i][j].wave;  /*pass point bits to lasFIle structure*/
        als->lasIn[i]->psID=als->point[i][j].psID;  /*pass point bits to lasFIle structure*/
        tlsWave=tlsWithinBeam(nVox,voxList,als->lasIn[i],tls,nTLS,x0,y0,z0,dimage,vox,xCent,yCent,zCent,denoised);
        als->lasIn[i]->wave=NULL;

        /*output the waves*/
        if(dimage->waveOut){
          outputWaves(als->point[i][j].wave,tlsWave,als->point[i][j].waveLen,dimage->outRoot,\
                      dimage->counter,dimage,xCent,yCent,zCent,als->lasIn[i],denoised,dimage->binOff,\
                      discRange,als->point[i][j].refl,als->lasIn[i]->gbic[als->point[i][j].psID]);
        }

        /*compare beams*/
        if(dimage->compBeams)compareBeams(vox,als->point[i][j].wave,area,als->point[i][j].waveLen,\
                                          xCent,yCent,zCent,als->lasIn[i],dimage,&tlsRMSE,minGap,appRefl);

        /*build voxels for comparison*/
        if(dimage->compVox)buildVoxel(alsVox,area,xCent,yCent,zCent,als->lasIn[i],dimage->beamRad,tlsRMSE,minGap,appRefl,vox);
    

        dimage->counter++; /*this is also used within tlsWithinBeam()*/
        tidyDataStruct(vox->tlsDat);
        TTIDY((void **)tlsWave,dimage->nTLSwaves);
      }/*TLS data check*/

      TIDY(area);
      TIDY(voxList);
      TIDY(denoised);
      TIDY(als->point[i][j].wave);
    }/*ALS point loop*/
  }/*ALS file loop*/

  /*tidy tls optimisation structure*/
  if(dimage->compVox){
    setTopVoxBlank(alsVox);
    compareVox(alsVox,vox->vox,dimage,vox);
    alsVox=tidyVox(alsVox);
  }

  if(dimage->compBeams){
    vox->tlsDat=tidyOverStruct(vox->tlsDat);
    fprintf(stdout,"Beam errors to %s\n",dimage->errOut);
  }
  TIDY(minGap);
  TIDY(appRefl);
  return;
}/*matchWaves*/


/*############################################*/
/*compare voxel models*/

void compareVox(voxStruct *alsVox,voxStruct *tlsVox,control *dimage,tlsVoxStruct *vox)
{
  int i=0,j=0,k=0;
  int place=0,numb=0;
  float alsCov=0,tlsCov=0;
  float totALS=0,totTLS=0;
  float contN=0;
  float *tlsProj=NULL;
  float *projectVox(tlsVoxStruct *,int,int);
  double x=0,y=0,z=0;
  char namen[200];
  FILE *opoo=NULL;

  /*check that both are the same*/
  if((alsVox->nX!=tlsVox->nX)||(alsVox->nY!=tlsVox->nY)||(alsVox->nZ!=tlsVox->nZ)){
    fprintf(stderr,"Voxel number mismatch\n");
    exit(1);
  }

  sprintf(namen,"%s.vErr",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"# 1 tlsCov, 2 alsCov, 3 totTLS , 4 totALS, 5 x, 6 y, 7 z, 8 tlsRMSE\n");

  for(i=0;i<alsVox->nX;i++){
    x=(double)(i+0.5)*(double)alsVox->res[0]+alsVox->bounds[0];
    for(j=0;j<alsVox->nY;j++){
      y=(double)(j+0.5)*(double)alsVox->res[1]+alsVox->bounds[1];

      /*condense the TLS data to calcuate projected area*/
      //tlsProj=projectVox(vox,i,j);

      for(k=0;k<alsVox->nZ;k++){
        z=(double)(k+0.5)*(double)alsVox->res[2]+alsVox->bounds[2];

        place=k*alsVox->nY*alsVox->nX+j*alsVox->nX+i;

        /*tls area*/
        tlsCov=contN=0.0;
        for(numb=0;numb<tlsVox->nScans;numb++){
          totTLS=tlsVox->hits[numb][place]+tlsVox->miss[numb][place];
          if(totTLS>0.0){
            tlsCov+=tlsVox->hits[numb][place];
            contN+=totTLS;
          }
        }
        if(contN>0.0)tlsCov/=contN;
        else         tlsCov=-1.0;

        /*als area*/
        if(alsVox->contN[place]>0){
          alsVox->hits[0][place]/=(float)alsVox->contN[place];
          alsVox->miss[0][place]/=(float)alsVox->contN[place];
          alsVox->rmse[place]/=(float)alsVox->contN[place];
          totALS=alsVox->hits[0][place]+alsVox->miss[0][place];
        }else totALS=0.0;
        if(totALS>0.0)alsCov=alsVox->hits[0][place]/totALS;
        else          alsCov=-1.0;

        fprintf(opoo,"%f %f %f %d %f %f %f %f\n",tlsCov,alsCov,contN,alsVox->contN[place],x,y,z,alsVox->rmse[place]);
      }/*z loop*/
      TIDY(tlsProj);
    }/*y loop*/
  }/*x loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Voxel errors to %s\n",namen);

  return;
}/*compareVox*/


/*############################################*/
/*project TLS voxels upwards*/

float *projectVox(tlsVoxStruct *vox,int i,int j)
{
  int k=0,n=0,place=0;
  int contN=0;
  float *tlsProj=NULL;
  void makeTLSwaves(overStruct *);
  overStruct *tempDat=NULL;
  overStruct *tidyOverStruct(overStruct *);
  void setPars(overStruct *);



  /*copy data into arrays for making TLS waveform*/
  if(!(tempDat=(overStruct *)calloc(1,sizeof(overStruct)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  if(!(tempDat->dimage=(optControl *)calloc(1,sizeof(optControl)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  tempDat->dimage->nFiles=vox->vox->nScans;
  tempDat->dimage->res=(float)vox->vox->res[2];
  tempDat->dimage->errOut=NULL;
  setPars(tempDat);

  /*normalise meanRefl and meanGap and load into structure*/
  for(n=0;n<vox->vox->nScans;n++){
    tempDat->dArr[n].tlsPar.appRefl=0.0;
    tempDat->dArr[n].tlsPar.minGap=0.0;
    contN=0;
    for(k=0;k<vox->vox->nZ;k++){
      place=k*vox->vox->nY*vox->vox->nX+j*vox->vox->nX+i;
      if(vox->contN[n][place]>0){
        tempDat->dArr[n].tlsPar.appRefl+=vox->meanRefl[n][place]/(float)vox->contN[n][place];
        tempDat->dArr[n].tlsPar.minGap+=vox->meanGap[n][place]/(float)vox->contN[n][place];
        contN++;
      }
    }
    if(contN>0){
      tempDat->dArr[n].tlsPar.appRefl/=(float)contN;
      tempDat->dArr[n].tlsPar.minGap/=(float)contN;
    }
  }

  /*load points into structure*/
  if(!(tempDat->dArr=(dataStruct *)calloc(tempDat->dimage->nFiles,sizeof(dataStruct)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }


  /*make TLS waveform*/
  makeTLSwaves(tempDat);

  /*copy relevant bit*/

  /*tidy up*/
  tempDat=tidyOverStruct(tempDat);

  return(tlsProj);
}/*projectVox*/


/*############################################*/
/*build up ALS voxels*/

void buildVoxel(voxStruct *alsVox,float *area,double xCent,double yCent,double zCent,lasFile *lasIn,float beamRad,float tlsRMSE,float *minGap,float *appRefl,tlsVoxStruct *vox)
{
  int i=0,j=0,zInd=0;
  int lastBin=0;
  int *totList=NULL,nTot=0;
  char *passList=NULL;
  float *thisGap=NULL;
  double zTop=0,zBot=0;
  double x0=0,y0=0,z0=0;


  /*determine when the beam is blocked*/
  lastBin=lasIn->waveLen+1;
  for(i=0;i<lasIn->waveLen;i++){
    if(area[i]>0.0)lastBin=i;
  }

  /*map for all beams*/
  binPosition(&x0,&y0,&z0,0,xCent,yCent,zCent,lasIn->time,lasIn->grad);
  totList=beamVoxels(lasIn->grad,x0,y0,z0,&(alsVox->bounds[0]),&(alsVox->res[0]),alsVox->nX,alsVox->nY,alsVox->nZ,&nTot,beamRad,NULL,-1.0);
  passList=challoc(nTot,"pass list",0);
  for(j=0;j<nTot;j++)passList[j]=0;

  /*gap fraction for each voxel for this beam*/
  thisGap=falloc(nTot,"thisGap",0);
  for(i=0;i<nTot;i++)thisGap[i]=1.0;

  /*loop along until last beam*/
  for(i=0;i<lastBin;i++){
    binPosition(&x0,&y0,&z0,i,xCent,yCent,zCent,lasIn->time,lasIn->grad);
    for(j=0;j<nTot;j++){
      zInd=(int)((float)totList[j]/(float)(alsVox->nX*alsVox->nY));
      zBot=(double)(zInd)*(double)alsVox->res[2]+alsVox->bounds[2];
      zTop=(double)(zInd+1)*(double)alsVox->res[2]+alsVox->bounds[2];
      if((z0<zTop)&&(z0>zBot)){  /*check that this voxel is within bin*/
        if(area[i]>0.0000001){
          /*fprintf(stdout,"Found one\n");*/
          alsVox->hits[0][totList[j]]+=area[i];
          thisGap[j]-=area[i];
        }
        passList[j]=1;
      }/*else fprintf(stdout,"WIthout %f %f %f %d %d %f\n",z0,zBot,zTop,zInd,alsVox->nZ,alsVox->bounds[4]);*/
    }/*voxel loop*/
  }/*bin loop*/

  /*we need to normalise by number of beams*/
  for(i=0;i<nTot;i++){
    if(passList[i]&&totList[i]){    /*if this voxel was passed through*/
      alsVox->miss[0][totList[i]]+=thisGap[i];
      alsVox->rmse[totList[i]]+=tlsRMSE;
      alsVox->contN[totList[i]]++;
      /*TLS parameters for projected area*/
      for(j=0;j<vox->vox->nScans;j++){
        vox->meanGap[j][totList[i]]+=minGap[j];
        vox->meanRefl[j][totList[i]]+=appRefl[j];
        vox->contN[j][totList[i]]++;
      }
    }
  }

  TIDY(thisGap);
  TIDY(totList);
  TIDY(passList);
  return;
}/*buildVoxel*/


/*############################################*/
/*compare ALS and TLS beams*/

void compareBeams(tlsVoxStruct *vox,unsigned char *rawALS,float *area,int numb,double xCent,double yCent,double zCent,lasFile *lasIn,control *dimage,float *tlsPass,float *minGap,float *appRefl)
{
  int i=0,hBin=0,gBin=0;
  int nGbins=0,nHbins=0;
  int minHbin=0,maxHbin=0;
  int cBin=0;
  float tls=0,als=0;            /*projected areas*/
  float err=0,gap=0;
  float aE=0,raE=0,rtE=0;       /*waveform energies*/
  float tlsRMSE=0,alsRMSE=0;    /*raw ALS and process ALS RMSEs*/
  float *hErr=NULL,*gErr=NULL;  /*error with height and gap arrays*/
  float *smoothALS=NULL;        /*smoothed ALS to compare*/
  float falseNeg=0,falsePos=0;  /*false negative and positive errors*/
  float *coarseTLS=NULL,*coarseALS=NULL;
  double x=0,y=0,z=0;      /*coordinates for matching to ground map*/
  double gZ=0;
  double setCanGround(double,double,canBstruct *);
  void fitTLSparams(overStruct *);
  void makeTLSwaves(overStruct *);
  void setPars(overStruct *);
  void setDenoiseDefault(denPar *);
  char *gIn=NULL,*hIn=NULL;


  if(vox->tlsDat->dArr->nTLS==0)return;  /*no TLS data*/

  /*make TLS waveform*/
  vox->tlsDat->dArr->nALS=numb;
  setPars(vox->tlsDat);
  if(vox->tlsDat->dArr->nALS==0)return;
  vox->tlsDat->dArr->als=processWave(rawALS,numb,&(vox->tlsDat->den),1.0);
  fitTLSparams(vox->tlsDat);  /*tune TLS*/
  makeTLSwaves(vox->tlsDat);  /*make waveform. Will be true area*/

  /*copy TLS parameters*/
  for(i=0;i<vox->vox->nScans;i++){
    appRefl[i]=vox->tlsDat->dArr->tlsPar.appRefl;
    minGap[i]=vox->tlsDat->dArr->tlsPar.minGap;
  }

  /*smooth ALS for error*/
  smoothALS=smooth(vox->tlsDat->dimage->denWidth,numb,area,vox->tlsDat->dimage->res);
  aE=raE=rtE=0.0;
  for(i=0;i<numb;i++){
    aE+=smoothALS[i];
    raE+=vox->tlsDat->dArr->als[i];
    rtE+=vox->tlsDat->dArr->tlsWave[0][i];
  }

  /*TLS to raw ALS fit error*/
  tlsRMSE=0.0;
  for(i=0;i<numb;i++){
    err=vox->tlsDat->dArr->als[i]/raE-vox->tlsDat->dArr->tlsWave[0][i]/rtE;
    tlsRMSE+=err*err;
  }
  tlsRMSE=sqrt(tlsRMSE/(float)numb);

  /*processed ALS to projected TLS area error*/
  alsRMSE=0.0;
  /*set up arrays for height error*/
  nHbins=(int)((float)numb*vox->tlsDat->dimage->res/dimage->hErrRes+0.5)+1;
  hErr=falloc(nHbins,"herror",0);
  hIn=challoc(nHbins,"height in",0);
  for(i=0;i<nHbins;i++){
    hErr[i]=0.0;
    hIn[i]=0;
  }
  minHbin=1000;
  maxHbin=-1000;
  /*set up arrays for gap error*/
  nGbins=(int)(1.0/dimage->gRes)+2;
  gErr=falloc(nGbins,"gap error",0);
  gIn=challoc(nGbins,"gap in",0);
  for(i=0;i<nGbins;i++){
    gErr[i]=0.0;
    gIn[i]=0;
  }
  gap=1.0;
  /*coarsened waves for ommission/commission errors*/
  coarseALS=falloc(nHbins,"coarsened ALS",0);
  coarseTLS=falloc(nHbins,"coarsened TLS",0);
  for(i=0;i<nHbins;i++)coarseALS[i]=coarseTLS[i]=0.0;
  for(i=0;i<numb;i++){
    /*total RMSE*/
    als=smoothALS[i]/aE;
    tls=vox->tlsDat->dArr->smoothTLS[i]/vox->tlsDat->dArr->tE;
    err=als-tls;
    alsRMSE+=err*err;

    /*error as a function of height above ground*/
    binPosition(&x,&y,&z,i,xCent,yCent,zCent,lasIn->time,lasIn->grad);
    if(dimage->useCanBound)gZ=setCanGround(x,y,&dimage->canB);
    else                   gZ=0.0;
    hBin=(int)((z-gZ)/(double)dimage->hErrRes+0.5);
    if((hBin>=0)&&(hBin<nHbins)){
      hErr[hBin]+=err;
      if((als>0.0)||(tls>0.0))hIn[hBin]=1;
      if(hBin<minHbin)minHbin=hBin;
      if(hBin>maxHbin)maxHbin=hBin;
    }

    /*error as a function of gap fraction*/
    gBin=nGbins-((int)(gap/dimage->gRes+0.5)+1);
    if((gBin<0)||(gBin>=nGbins)){
      fprintf(stderr,"Bin error from gap %f bin %d of %d\n",gap,gBin,nGbins);
      exit(1);
    }
    gErr[gBin]+=err;
    if((als>0.0)||(tls>0.0))gIn[gBin]=1;
    gap-=tls;

    /*coarse waveforms for ommission/commission*/
    cBin=(int)((float)(i*nHbins)/(float)numb);
    if((cBin>=0)&&(cBin<nHbins)){
      coarseALS[cBin]+=als;
      coarseTLS[cBin]+=tls;
    }
  }
  alsRMSE=sqrt(alsRMSE/(float)numb);

  /*ommission/commission error*/
  falsePos=falseNeg=0.0;
  for(i=0;i<nHbins;i++){
    if((coarseALS[i]>0.0)&&(coarseTLS[i]==0.0))falsePos+=coarseALS[i];
    if((coarseALS[i]==0.0)&&(coarseTLS[i]>0.0))falseNeg+=coarseTLS[i];
  }
  TIDY(coarseALS);
  TIDY(coarseTLS);

  /*write out results*/
  fprintf(vox->tlsDat->dimage->errOut,"%f %f %f %f %d %f %f %f %f",tlsRMSE,alsRMSE,\
                           falsePos,falseNeg,numb,aE,raE,rtE,vox->tlsDat->dArr->tE);
  for(i=0;i<nGbins;i++){
    if(gIn[i])fprintf(vox->tlsDat->dimage->errOut," %f",gErr[i]);
    else      fprintf(vox->tlsDat->dimage->errOut," ?");
  }
  for(i=maxHbin;i>=minHbin;i--){
    if(hIn[i])fprintf(vox->tlsDat->dimage->errOut," %f",hErr[i]);
    else      fprintf(vox->tlsDat->dimage->errOut," ?");
  }
  /*for(i=wEnd-wStart;i<numb;i++)fprintf(vox->tlsDat->dimage->errOut," ?");*/
  fprintf(vox->tlsDat->dimage->errOut,"\n");

  *tlsPass=tlsRMSE;

  /*tidy arrays*/
  TIDY(hIn);
  TIDY(gIn);
  TIDY(hErr);
  TIDY(gErr);
  TIDY(smoothALS);
  TIDY(vox->tlsDat->dArr->als);
  return;
}/*compareBeams*/


/*############################################*/
/*tidy TLS structure*/

void tidyDataStruct(overStruct *data)
{
  int i=0;

  if(data){
    if(data->dArr){
      for(i=0;i<data->dimage->nFiles;i++){
        TIDY(data->dArr[i].tls);
        TIDY(data->dArr[i].als);
        TIDY(data->dArr[i].rawALS);
        TTIDY((void **)data->dArr[i].tlsWave,3);
        TIDY(data->dArr[i].smoothTLS);
      }
    }
    TIDY(data->dArr);
  }
  return;
}/*tidyDataStruct*/


/*############################################*/
/*set up TLS optimisation structure*/

overStruct *setupOverStruct(control *dimage)
{
  int i=0,nGbins=0;
  overStruct *data=NULL;
  void pulseReadForTLS(char *,overStruct *);
  void setDenoiseDefault(denPar *);

  if(!(data=(overStruct *)calloc(1,sizeof(overStruct)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }

  /*read pulse for convolution*/
  pulseReadForTLS(dimage->pNamen,data);

  data->dimage=&dimage->ofimage;
  data->dimage->nFiles=1;
  data->dimage->res=0.15;  /*this may also be set in pulseReadForTLS()*/

  /*use defaults for TLS optimisation*/
  setDenoiseDefault(&data->den);
  data->den.medStats=1;

  /*open output file*/
  sprintf(dimage->errOut,"%s.err",dimage->outRoot);
  if((dimage->ofimage.errOut=fopen(dimage->errOut,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->errOut);
    exit(1);
  }

  fprintf(dimage->ofimage.errOut,"# 1 tlsRMSE, 2 alsRMSE, 3 falsePos, 4 falseNeg, 5 numb, 6 aE, 7 raE, 8 rtE, 9 tE");
  nGbins=(int)(1.0/dimage->gRes)+1;
  for(i=0;i<nGbins;i++)fprintf(dimage->ofimage.errOut,", %d gErr %d",i+10,i);
  for(i=0;i<256;i++)fprintf(dimage->ofimage.errOut,", %d hErr %d",i+10+nGbins,i);
  fprintf(dimage->ofimage.errOut,"\n");

  return(data);
}/*setupOverStruct*/


/*############################################*/
/*set up TLS structure*/

void setupTLSstruct(overStruct *data)
{
  data->thisTLS=0;        /*only one TLS structure at a time*/

  if(!(data->dArr=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  data->dArr->nTLS=0;
  data->dArr->nALS=0;
  data->dArr->rawALS=NULL;
  data->dArr->als=NULL;

  return;
}/*setupTLSstruct*/


/*############################################*/
/*write waveforms from different sources*/

void outputWaves(unsigned char *wave,float **tlsWave,int32_t waveLen,char *outRoot,uint64_t nOut,control *dimage,double lon,double lat,double elev,lasFile *lasIn,float *denoised,int binOff,double discRange,uint16_t discRefl,float gbic)
{
  int i=0,j=0;
  char namen[200];
  double xWave=0,yWave=0,zWave=0;
  FILE *opoo=NULL;

  sprintf(namen,"%s.%lld.dat",outRoot,(long long)nOut);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  fprintf(opoo,"# 1 range, 2 ALS, 3 TLS count, 4 TLS visible cover, 5 TLS bin cover, 6 TLS summit, 7 denoised, 8 lon, 9 lat, 10 elev, 11 gbic\n");
  fprintf(opoo,"# lon %f lat %f elev %f\n",lon,lat,elev);
  fprintf(opoo,"# discRange %f %d\n",discRange,(int)discRefl);
  for(i=0-binOff;i<0;i++){
    fprintf(opoo,"%f ?",(float)i*0.15);
    for(j=0;j<dimage->nTLSwaves;j++)fprintf(opoo," %f",tlsWave[j][i+binOff]);
    binPosition(&xWave,&yWave,&zWave,i,lon,lat,elev,lasIn->time,lasIn->grad);
    fprintf(opoo," ? %f %f %f ?",xWave,yWave,zWave);
    fprintf(opoo,"\n");
  }
  for(i=0;i<(int)waveLen;i++){
    fprintf(opoo,"%f %d",(float)i*0.15,wave[i]);
    for(j=0;j<dimage->nTLSwaves;j++)fprintf(opoo," %f",tlsWave[j][i+binOff]);
    binPosition(&xWave,&yWave,&zWave,i,lon,lat,elev,lasIn->time,lasIn->grad);
    fprintf(opoo," %f %f %f %f %f",denoised[i],xWave,yWave,zWave,gbic);
    fprintf(opoo,"\n");
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",namen);
  return;
}/*outputWaves*/


/*####################################*/
/*make TLS waves within ALS*/

float **tlsWithinBeam(int nVox,int *voxList,lasFile *tempLas,tlsScan *tls,int nTLS,double x0,double y0,double z0,control *dimage,tlsVoxStruct *vox,double xDisc,double yDisc,double zDisc,float *denoised)
{
  int i=0,j=0;
  int tF=0,tP=0;                /*tls file and point indices*/
  int rNx=0,rNy=0;              /*range image dimensions*/
  int bin=0;
  float **tlsWave=NULL;        /*multiple TLS waves. 0 point count, 1 projected cover*/
  float sepSq=0,minSepSq=0;
  float zen=0,az=0;
  char **rImage=NULL;    /*range image, a stack nBins long*/
  float rimRes=0;         /*range image resolution*/
  float gap=0;            /*point radius and gap fraction to that voxel*/
  double *vect=NULL;
  double x=0,y=0,z=0;
  void rotateX(double *,double);
  void rotateZ(double *,double);
  void saveTLSpoint(overStruct *,double *,float,float,float,int);
  void makeBinImage(double *,float *,float *,char *,int,float,double,float,float,int,int,float,lidVoxPar *);
  void waveFromImageOld(char **,float **,int,int,int);;
  FILE *waveBits=NULL;     /*file to write bits to make TLS waveform from*/
  char bitsNamen[200];     /*filename for the above*/
  char found=0;


  /*For outputting individual points*/
  double xWave=0,yWave=0,zWave=0;
  FILE *opoo=NULL;
  char namen[200];
  if(dimage->pointOut){ /*if outputting open up a file*/
    sprintf(namen,"%s.%lld.pts",dimage->outRoot,(long long)dimage->counter);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
  }/*output points or not*/

  if(dimage->waveBitsOut){  /*if we are to output these, open the file*/
    sprintf(bitsNamen,"%s.%lld.wBits",dimage->outRoot,(long long)dimage->counter);
    if((waveBits=fopen(bitsNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",bitsNamen);
      exit(1);
    }
  }/*TLS waveform bits file*/


  /*allocate arrays*/
  vect=dalloc(3,"vector",0);
  tlsWave=fFalloc(dimage->nTLSwaves,"tls waveform",0);   /*this is multiple waves now*/
  for(i=0;i<dimage->nTLSwaves;i++)tlsWave[i]=falloc(tempLas->waveLen+dimage->binOff,"tls waveform",i+1);

  /*range image bits*/
  rimRes=0.001;   /*1 mm*/
  rNx=rNy=(int)(2.0*dimage->beamRad/rimRes);
  rImage=chChalloc(tempLas->waveLen+dimage->binOff,"",0);
  for(j=0;j<tempLas->waveLen+dimage->binOff;j++){
    rImage[j]=challoc(rNx*rNy,"range image",j+1);
    for(i=rNx*rNy-1;i>=0;i--)rImage[j][i]=0;
  }

  minSepSq=dimage->beamRad*dimage->beamRad;   /*beam radius squared*/
  zen=(float)atan2(sqrt((double)tempLas->grad[0]*(double)tempLas->grad[0]+\
      (double)tempLas->grad[1]*(double)tempLas->grad[1]),(double)tempLas->grad[2]);
  az=(float)atan2((double)tempLas->grad[0],(double)tempLas->grad[1]);

  for(i=0;i<nVox;i++){   /*TLS voxel loop*/
    for(j=0;j<vox->nIn[voxList[i]];j++){ /*points within voxel loop*/
      /*read pixel and load into array. Translate to put beam at origin*/
      tF=vox->mapFile[voxList[i]][j];
      tP=vox->mapPoint[voxList[i]][j];

      x=(double)tls[tF].point[tP].x+tls[tF].xOff;
      y=(double)tls[tF].point[tP].y+tls[tF].yOff;
      z=(double)tls[tF].point[tP].z+tls[tF].zOff;
      /*fprintf(stdout,"tlsAfter %f %f %f %d\n",x,y,z,tF);*/

      vect[0]=(double)x-x0;
      vect[1]=(double)y-y0;
      vect[2]=(double)z-z0;

      /*rotate to x-y plane*/
      rotateZ(vect,(double)(-1.0*az));
      rotateX(vect,(double)(-1.0*zen));

      /*calculate distance from beam*/
      sepSq=(float)(vect[0]*vect[0]+vect[1]*vect[1]);
      if(sepSq<=minSepSq){

        bin=(int)(vect[2]/0.15+0.5)+dimage->binOff;
        if((bin>=0)&&(bin<(tempLas->waveLen+dimage->binOff))){  /*within beam*/
          if((vox->vox->miss[tF][voxList[i]]+vox->vox->hits[tF][voxList[i]])>0.0){
            gap=vox->toTLS->miss[tF][voxList[i]]/(vox->toTLS->miss[tF][voxList[i]]+vox->toTLS->hits[tF][voxList[i]]);
          }else gap=-1.0;
          makeBinImage(&(vect[0]),&(tlsWave[0][bin]),&(tlsWave[3][bin]),&(rImage[bin][0]),\
                       tempLas->waveLen+dimage->binOff,gap,(double)tls[tF].point[tP].r,\
                       tls[tF].point[tP].refl,dimage->beamRad,rNx,rNy,rimRes,&dimage->tlsPar);

          /*write point cloud*/
          if(dimage->pointOut)fprintf(opoo,"%f %f %f ? ? ?\n",x,y,z);

          /*write out data to make waveform for external optimisation*/
          if(dimage->waveBitsOut)fprintf(waveBits,"%d %f %f %f %f %f %d %f\n",bin,vect[0],vect[1],vect[2],\
             gap,tls[tF].point[tP].r,(int)tls[tF].point[tP].refl,vox->toTLS->miss[tF][voxList[i]]+vox->toTLS->hits[tF][voxList[i]]);
          found=1;

          /*save to array for optimisation*/
          if(dimage->compBeams)saveTLSpoint(vox->tlsDat,vect,gap,tls[tF].point[tP].r,tls[tF].point[tP].refl,bin);
        }/*within beam*/
      }/*separation check*/
    }/*tls point loop*/
  }/*map voxel loop*/

  /*condense image into waveform by distance along z axis*/
  waveFromImageOld(rImage,&(tlsWave[1]),tempLas->waveLen+dimage->binOff,rNx,rNy);


  /*output ALS waveform as individual points*/
  if(dimage->pointOut){
    for(i=0;i<tempLas->waveLen;i++){
      if(tempLas->wave[i]>16){
        binPosition(&xWave,&yWave,&zWave,i,xDisc,yDisc,zDisc,tempLas->time,tempLas->grad);
        fprintf(opoo,"? ? ? %f %f %f\n",xWave,yWave,zWave);
      }
    }
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Points written to %s\n",namen);
  }/*output individual points*/

  /*waveform bits*/
  if(dimage->waveBitsOut){
    /*write ALS waveform*/
    if(found){
      for(bin=0;bin<tempLas->waveLen;bin++)fprintf(waveBits,"ALS %d %d %f %f\n",\
                       bin+dimage->binOff,tempLas->wave[bin],denoised[bin],tempLas->gbic[tempLas->psID]);
      fprintf(waveBits,"# %.2f %.2f %.2f\n",x0,y0,z0);
      fprintf(stdout,"Waveform bits written to %s\n",bitsNamen);
    }
    if(waveBits){
      fclose(waveBits);
      waveBits=NULL;
    }
  }/*waveform bits*/

  TIDY(vect);
  TTIDY((void **)rImage,tempLas->waveLen+dimage->binOff);
  return(tlsWave);
}/*tlsWithinBeam*/


/*#######################################*/
/*save TLS point for later optimisation*/

void saveTLSpoint(overStruct *tlsDat,double *vect,float gap,float r,float refl,int bin)
{

  if(!(tlsDat->dArr->tls=(tlsPoint *)realloc(tlsDat->dArr->tls,(tlsDat->dArr->nTLS+1)*sizeof(tlsPoint)))){
    fprintf(stderr,"Balls\n");
    exit(1);
  }

  tlsDat->dArr->tls[tlsDat->dArr->nTLS].bin=bin;       /*bin number*/
  tlsDat->dArr->tls[tlsDat->dArr->nTLS].x=(float)vect[0];       /*coordinate*/
  tlsDat->dArr->tls[tlsDat->dArr->nTLS].y=(float)vect[1];       /*coordinate*/
  tlsDat->dArr->tls[tlsDat->dArr->nTLS].z=(float)vect[2];       /*coordinate*/
  tlsDat->dArr->tls[tlsDat->dArr->nTLS].gap=gap;     /*voxel gap fraction*/
  tlsDat->dArr->tls[tlsDat->dArr->nTLS].r=r;       /*range*/
  tlsDat->dArr->tls[tlsDat->dArr->nTLS].refl=(uint16_t)(refl); /*reflectance*/
  tlsDat->dArr->nTLS++;
  return;
}/*saveTLSpoint*/


/*####################################*/
/*read TLS within plot*/

tlsScan *readTLSwithinn(listLas *tlsList,control *dimage,tlsVoxStruct *vox)
{
  int i=0,k=0,n=0;
  int xBin=0,yBin=0,zBin=0;
  uint32_t j=0;
  uint32_t nUsed=0;
  uint32_t nBuff=0;        /*a buffer to deal with multiple point hits*/
  tlsScan *tls=NULL,*tempTLS=NULL;
  double xCent=0,yCent=0,zCent=0;   /*beam origin coordinates*/
  double x=0,y=0,z=0;               /*hit coordinates*/
  double minX=0,minY=0,minZ=0;
  double maxX=0,maxY=0,maxZ=0;
  double maxR=0,lastR=0;
  double grad[3];
  double *rangeList=NULL;
  char pointingTowards(double,double,double,tlsBeam *,double *);
  void setTLSpoint(float,float,float,float,float,uint8_t,uint8_t,tlsPoint *);
  void mapVox(double,double,double,int,uint32_t,tlsVoxStruct *);
  int nPix=0;
  int *tempList=NULL;
  char found=0;

  if(!(tls=(tlsScan *)calloc(tlsList->nFiles,sizeof(tlsScan)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }

  /*Riegl maximum range for bounds check*/
  maxR=300.0;

  /*loop over TLS files*/
  for(i=0;i<tlsList->nFiles;i++){
    /*read a scan*/
    tempTLS=readTLSpolarBinary(tlsList->nameList[i]);
    nBuff=3*tempTLS->nBeams;
    if(!(tls[i].point=(tlsPoint *)calloc(nBuff,sizeof(tlsPoint)))){
      fprintf(stderr,"error tls allocation.\n");
      exit(1);
    }
    tls[i].beam=NULL;
    tls[i].xOff=tempTLS->xOff;
    tls[i].yOff=tempTLS->yOff;       
    tls[i].zOff=tempTLS->zOff;       
    nUsed=0;

    /*are we within 300 m of the bounds?*/
    if(((dimage->bounds[0]-tempTLS->xOff)>maxR)||((dimage->bounds[1]-tempTLS->yOff)>maxR)||
       ((dimage->bounds[2]-tempTLS->zOff)>maxR)||((dimage->bounds[3]-tempTLS->xOff)<(-1.0*maxR))||
       ((dimage->bounds[4]-tempTLS->yOff)<(-1.0*maxR))||((dimage->bounds[5]-tempTLS->zOff)<(-1.0*maxR)))continue;


    /*loop through beams and see if they intersect the plot*/
    for(j=0;j<tempTLS->nBeams;j++){
      xCent=(double)tempTLS->beam[j].x+tempTLS->xOff;
      yCent=(double)tempTLS->beam[j].y+tempTLS->yOff;
      zCent=(double)tempTLS->beam[j].z+tempTLS->zOff;

      /*##### VOXELS #####*/
      /*map voxels*/
      grad[0]=tempTLS->beam[j].zen;
      grad[1]=tempTLS->beam[j].az;
      grad[2]=-99999.0;
      nPix=0;
      tempList=findVoxels(&(grad[0]),xCent,yCent,zCent,vox->vox->bounds,\
               &(vox->vox->res[0]),&nPix,vox->vox->nX,vox->vox->nY,vox->vox->nZ,&rangeList);

      if(tempTLS->beam[j].nHits>0)lastR=(double)tempTLS->beam[j].r[tempTLS->beam[j].nHits-1];
      else                        lastR=maxR;

      /*count up hits and misses*/
      for(k=0;k<nPix;k++){ /*loop over hits in beam*/
        if((tempList[k]<0)||(tempList[k]>vox->vox->nVox)){
          fprintf(stdout,"Index error %d\n",tempList[k]);
          exit(1);
        }

        /*gap fraction within voxel*/
        if(rangeList[k]<=lastR){   /*check last hit range*/
          if(tempTLS->beam[j].nHits==0)vox->vox->miss[i][tempList[k]]+=1.0; /*no returns*/
          else{      /*see if there is a hit within*/
            zBin=tempList[k]/(vox->vox->nX*vox->vox->nY);
            yBin=(tempList[k]-zBin*vox->vox->nX*vox->vox->nY)/vox->vox->nX;
            xBin=tempList[k]-(zBin*vox->vox->nX*vox->vox->nY+yBin*vox->vox->nX);
            minX=(float)xBin*vox->vox->res[0]+vox->vox->bounds[0];
            minY=(float)yBin*vox->vox->res[1]+vox->vox->bounds[1];
            minZ=(float)zBin*vox->vox->res[2]+vox->vox->bounds[2];
            maxX=(float)(xBin+1)*vox->vox->res[0]+vox->vox->bounds[0];
            maxY=(float)(yBin+1)*vox->vox->res[1]+vox->vox->bounds[1];
            maxZ=(float)(zBin+1)*vox->vox->res[2]+vox->vox->bounds[2];
            found=0;
            for(n=0;n<tempTLS->beam[j].nHits;n++){
              x=xCent+tempTLS->beam[j].r[n]*sin(tempTLS->beam[j].az)*sin(tempTLS->beam[j].zen);
              y=yCent+tempTLS->beam[j].r[n]*cos(tempTLS->beam[j].az)*sin(tempTLS->beam[j].zen);
              z=zCent+tempTLS->beam[j].r[n]*cos(tempTLS->beam[j].zen);
              if((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)&&(z>=minZ)&&(z<=maxZ)){
                vox->vox->hits[i][tempList[k]]+=1.0;
                found=1;
              }
            }
            if(found)vox->vox->miss[i][tempList[k]]+=1.0; /*if none within look to see if hit after*/
          }
        }/*last hit range check*/

        /*gap fraction to voxel*/
        if(rangeList[k]>maxR)break;  /*gone beyond maximum range. No data*/
        if(tempTLS->beam[j].nHits==0)                                              vox->toTLS->miss[i][tempList[k]]+=1.0; /*no returns*/
        else if((double)tempTLS->beam[j].r[tempTLS->beam[j].nHits-1]<=rangeList[k])vox->toTLS->hits[i][tempList[k]]+=1.0; /*hit before*/
        else                                                                       vox->toTLS->miss[i][tempList[k]]+=1.0; /*hit after*/
      }/*voxel intersect loop*/


      /*##### RETURN POINTS #####*/
      /*copy relevant Cartesian hits*/
      for(k=0;k<tempTLS->beam[j].nHits;k++){
        x=xCent+tempTLS->beam[j].r[k]*sin(tempTLS->beam[j].az)*sin(tempTLS->beam[j].zen);
        y=yCent+tempTLS->beam[j].r[k]*cos(tempTLS->beam[j].az)*sin(tempTLS->beam[j].zen);
        z=zCent+tempTLS->beam[j].r[k]*cos(tempTLS->beam[j].zen);

        /*check bounds and copy point if within*/
        if((x>=dimage->bounds[0])&&(y>=dimage->bounds[1])&&(z>=dimage->bounds[2])&&\
           (x<=dimage->bounds[3])&&(y<=dimage->bounds[4])&&(z<=dimage->bounds[5])){

          /*fprintf(stdout,"tls %f %f %f %d\n",x,y,z,i);*/
          if(nUsed>nBuff){
            fprintf(stderr,"Buffer not big enough %d %d\n",nUsed,nBuff);
            exit(1);
          }
          setTLSpoint((float)(x-tempTLS->xOff),(float)(y-tempTLS->yOff),(float)(z-tempTLS->zOff),\
             tempTLS->beam[j].r[k],tempTLS->beam[j].refl[k],k,tempTLS->beam[j].nHits,&(tls[i].point[nUsed]));

          /*map to voxels to ease later searches*/
          mapVox(x,y,z,i,nUsed,vox);
          nUsed++;
        }/*is hit within area*/
      }/*hit loop*/
      TIDY(rangeList);
      TIDY(tempList);
    }/*beam loop*/

    /*trim arrays*/
    if(nUsed>0){
      if(!(tls[i].point=(tlsPoint *)realloc(tls[i].point,nUsed*sizeof(tlsPoint)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
    }else TIDY(tls[i].point);
    tls[i].nPoints=nUsed;
    tls[i].nBeams=0;
    tempTLS=tidyTLScan(tempTLS);
  }/*file loop*/

  return(tls);
}/*readTLSwithinn*/


/*###################################################################*/
/*map voxels to TLS for later*/

void mapVox(double x,double y,double z,int nFile,uint32_t nBeam,tlsVoxStruct *vox)
{
  int xBin=0,yBin=0,zBin=0;
  int vPlace=0;
  int *markInt(int,int *,int);
  uint32_t *markUint32(int,uint32_t *,uint32_t);

  xBin=(int)((x-vox->vox->bounds[0])/vox->vox->res[0]);
  yBin=(int)((y-vox->vox->bounds[1])/vox->vox->res[1]);
  zBin=(int)((z-vox->vox->bounds[2])/vox->vox->res[2]);

  if((xBin>=0)&&(xBin<vox->vox->nX)&&(yBin>=0)&&(yBin<vox->vox->nY)&&(zBin>=0)&&(zBin<vox->vox->nZ)){
    vPlace=xBin+yBin*vox->vox->nX+zBin*vox->vox->nX*vox->vox->nY;
    vox->mapFile[vPlace]=markInt(vox->nIn[vPlace],&(vox->mapFile[vPlace][0]),nFile);
    vox->mapPoint[vPlace]=markUint32(vox->nIn[vPlace],&(vox->mapPoint[vPlace][0]),nBeam);
    vox->nIn[vPlace]++;
  }
  return;
}/*mapVox*/


/*####################################*/
/*copy a beam over*/

void copyBeam(tlsBeam *in,tlsBeam *out)
{
  uint8_t i=0;

  if(in->nHits>0){
    out->r=falloc(in->nHits,"range",0);
    out->refl=falloc(in->nHits,"refl",0);
  }else out->r=out->refl=NULL;

  out->zen=in->zen;
  out->az=in->az;
  out->x=in->x;
  out->y=in->y;
  out->z=in->z;
  out->shotN=in->shotN;
  out->nHits=in->nHits;
  for(i=0;i<out->nHits;i++){
    out->r[i]=in->r[i];
    out->refl[i]=in->refl[i];
  }

  return;
}/*copyBeam*/


/*####################################*/
/*copy a beam over*/

void setTLSpoint(float x,float y,float z,float r,float refl,uint8_t hitN,uint8_t nHits,tlsPoint *out)
{
  out->x=x;
  out->y=y;
  out->z=z;
  out->hitN=hitN;
  out->nHits=nHits;
  out->r=r;
  out->refl=refl;

  return;
}/*setTLSpoint*/


/*####################################*/
/*tidy ALS structure*/

void tidyALS(alsData *als)
{
  int i=0;

  if(als){
    TTIDY((void **)als->point,als->nFiles);
    TIDY(als->nPoints);
    for(i=0;i<als->nFiles;i++){
      tidyLasFile(als->lasIn[i]);
      als->lasIn[i]=NULL;
    }
    TIDY(als->lasIn);      /*las files*/
    TIDY(als);
  }
  return;
}/*tidyALS*/


/*####################################*/
/*intersect ALS with TLS*/

alsData *readALSwithin(listLas *alsList,control *dimage)
{
  int i=0,j=0;
  uint32_t nUsed=0;
  alsData *als=NULL;          /*las files*/
  double xCent=0,yCent=0,zCent=0;
  void copyALSpoint(lasPoint *,lasFile *);

  /*allocate space*/
  if(!(als=(alsData *)calloc(alsList->nFiles,sizeof(alsData)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  als->lasIn=lfalloc(alsList->nFiles);
  als->nFiles=alsList->nFiles;
  if(!(als->point=(lasPoint **)calloc(als->nFiles,sizeof(lasPoint *)))){
    fprintf(stderr,"error in als allocation.\n");
    exit(1);
  }
  if(!(als->nPoints=(uint32_t *)calloc(als->nFiles,sizeof(uint32_t)))){
    fprintf(stderr,"error in als allocation.\n");
    exit(1);
  }


  /*read file headers*/
  for(i=0;i<alsList->nFiles;i++)als->lasIn[i]=readLasHead(alsList->nameList[i],dimage->pBuffSize);

  /*read GBIC table*/
  readGBIC(dimage->appGBIC,dimage->balFlights,als->lasIn,alsList);


  /*file loop*/
  for(i=0;i<alsList->nFiles;i++){
    nUsed=0;
    if(!(als->point[i]=(lasPoint *)calloc(als->lasIn[i]->nPoints,sizeof(lasPoint)))){
      fprintf(stderr,"error tls allocation.\n");
      exit(1);
    }

    if(checkFileBounds(als->lasIn[i],dimage->bounds[0],dimage->bounds[3],dimage->bounds[1],dimage->bounds[4])){
      for(j=0;j<als->lasIn[i]->nPoints;j++){  /*point loop*/
        readLasPoint(als->lasIn[i],j);
        /*check whether it's the first point*/
        if(als->lasIn[i]->field.retNumb==1){
          /*check bounds*/
          setCoords(&xCent,&yCent,&zCent,als->lasIn[i]);
          if((xCent>=dimage->bounds[0])&&(xCent<=dimage->bounds[3])&&\
             (yCent>=dimage->bounds[1])&&(yCent<=dimage->bounds[4])){
            /*if satisfies then save*/
            /*fprintf(stdout,"als %f %f %f %d\n",xCent,yCent,zCent,i);*/
            copyALSpoint(&(als->point[i][nUsed]),als->lasIn[i]);
            nUsed++;
          }/*point bounds check*/
        }/*first return check*/
      }/*point loop*/
    }/*check file bounds*/
    /*trim array*/
    als->nPoints[i]=nUsed;
    if(nUsed>0){
      if(!(als->point[i]=(lasPoint *)realloc(als->point[i],nUsed*sizeof(lasPoint)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
    }else TIDY(als->point[i]);
    fprintf(stdout,"Using %d from ALS file %d\n",als->nPoints[i],i);
  }/*file loop*/
  return(als);
}/*readALSwithin*/


/*####################################*/
/*copy a las point*/

void copyALSpoint(lasPoint *point,lasFile *lasIn)
{
  point->x=lasIn->x;
  point->y=lasIn->y;
  point->z=lasIn->z;
  point->refl=lasIn->refl;
  point->field=lasIn->field;
  point->scanAng=lasIn->scanAng;
  point->psID=lasIn->psID;
  point->packetDes=lasIn->packetDes;
  point->waveMap=lasIn->waveMap;
  point->waveLen=lasIn->waveLen;
  point->time=lasIn->time;
  point->grad[0]=lasIn->grad[0];
  point->grad[1]=lasIn->grad[1];
  point->grad[2]=lasIn->grad[2];

  return;
}/*copyALSpoint*/


/*####################################*/
/*copy a las point to a lasFile structure*/

void writeALSpoint(lasPoint *point,lasFile *lasIn)
{
  lasIn->x=point->x;
  lasIn->y=point->y;
  lasIn->z=point->z;
  lasIn->refl=point->refl;
  lasIn->field=point->field;
  lasIn->scanAng=point->scanAng;
  lasIn->psID=point->psID;
  lasIn->packetDes=point->packetDes;
  lasIn->waveMap=point->waveMap;
  lasIn->waveLen=point->waveLen;
  lasIn->time=point->time;
  lasIn->grad[0]=point->grad[0];
  lasIn->grad[1]=point->grad[1];
  lasIn->grad[2]=point->grad[2];

  return;
}/*copyALSpoint*/


/*####################################*/
/*tidy voxel structure*/

tlsVoxStruct *tidyTLSvox(tlsVoxStruct *vox,int nFiles)
{
  overStruct *tidyOverStruct(overStruct *);

  if(vox){
    TTIDY((void **)vox->mapFile,vox->vox->nVox);
    TTIDY((void **)vox->mapPoint,vox->vox->nVox);
    TTIDY((void **)vox->meanGap,vox->vox->nScans);
    TTIDY((void **)vox->meanRefl,vox->vox->nScans);
    TTIDY((void **)vox->contN,vox->vox->nScans);
    vox->tlsDat=tidyOverStruct(vox->tlsDat);
    vox->vox=tidyVox(vox->vox);
    vox->toTLS=tidyVox(vox->toTLS);
    TIDY(vox->nIn);
    TIDY(vox);
  }
  return(NULL);
}/*tidyTLSvox*/


/*####################################*/
/*allocate tls voxel structure*/

tlsVoxStruct *tlsVoxAlloc(int nFiles,control *dimage)
{
  int i=0,j=0;
  tlsVoxStruct *vox=NULL;

  if(!(vox=(tlsVoxStruct *)calloc(1,sizeof(tlsVoxStruct)))){
    fprintf(stderr,"error voxel structure allocation.\n");
    exit(1);
  }

  /*voxel structures*/
  vox->vox=voxAllocate(nFiles,&(dimage->vRes[0]),dimage->bounds,0);
  vox->toTLS=voxAllocate(nFiles,&(dimage->vRes[0]),dimage->bounds,0);

  /*mean factors for projecting area upwards*/
  vox->meanRefl=fFalloc(vox->vox->nScans,"meanRefl",0);
  vox->meanGap=fFalloc(vox->vox->nScans,"meanGap",0);
  vox->contN=iIalloc(vox->vox->nScans,"cal count",0);
  for(i=0;i<vox->vox->nScans;i++){
    vox->meanRefl[i]=falloc(vox->vox->nVox,"meanRefl",i+1);
    vox->meanGap[i]=falloc(vox->vox->nVox,"meanGap",i+1);
    vox->contN[i]=ialloc(vox->vox->nVox,"cal count",i+1);
    for(j=vox->vox->nVox-1;j>=0;j--){
      vox->meanRefl[i][j]=vox->meanGap[i][j]=0.0;
      vox->contN[i][j]=0;
    }
  }

  /*TLS specific things to map back to individual scans*/
  vox->mapFile=iIalloc(vox->vox->nVox,"voxel file map",0);
  vox->nIn=ialloc(vox->vox->nVox,"number in",0);
  if(!(vox->mapPoint=(uint32_t **)calloc(vox->vox->nVox,sizeof(uint32_t *)))){
    fprintf(stderr,"error voxel map allocation.\n");
    exit(1);
  }

  return(vox);
}/*tlsVoxAlloc*/


/*####################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);


  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control structure allocation.\n");
    exit(1);
  }

  strcpy(dimage->alsNamen,"/home/sh563/data/teastlas1.3.teast.dat");
  strcpy(dimage->tlsNamen,"/home/sh563/data/teast/binary.list");
  strcpy(dimage->pNamen,"/Users/stevenhancock/data/bess/leica_shape/leicaPulse.dat");
  strcpy(dimage->outRoot,"teast");
  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=-100000000.0;
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=100000000.0;

  /*comparison*/
  dimage->compBeams=0;   /*do not compare individual beams*/
  dimage->compVox=0;     /*do not compare voxels*/

  /*TLS parameters for passing around*/
  dimage->nTLSwaves=4;  /*0 target profile, 1 visible projected area, 2 projected area, 3 point count*/
  dimage->beamRad=0.165;    /*33 cm footprints*/
  dimage->binOff=0;         /*start TLS wave from ALS wave*/
  dimage->counter=0;        /*reset counter*/
  dimage->pointOut=0;       /*don't output points*/
  dimage->waveOut=1;
  dimage->waveBitsOut=0;    /*don't output bits to make TLS waveforms*/
  dimage->pBuffSize=(uint64_t)200000000;
  /*TLS voxel parameters*/
  dimage->tlsPar.minGap=0.5;
  dimage->tlsPar.appRefl=0.2;  /*apparent reflectance, at 1545 nm actually*/
  dimage->tlsPar.beamTanDiv=sin(0.35/(2.0*1000.0))/cos(0.35/(2.0*1000.0));  /*tan half of 0.35 mrad*/
  dimage->tlsPar.beamRad=0.035;     /*7 mm diameter beam start*/

  /*GBIC*/
  dimage->appGBIC=0;
  dimage->balFlights=0;

  /*centre on plot 1 for the tests*/
  dimage->bounds[0]=507201.6157387340-10.0;
  dimage->bounds[1]=225391.1395131220-10.0;
  dimage->bounds[2]=127.0679805710-10.0;
  dimage->bounds[3]=507201.6157387340+10.0;
  dimage->bounds[4]=225391.1395131220+10.0;
  dimage->bounds[5]=127.0679805710+20.0;
  dimage->vRes[0]=dimage->vRes[1]=dimage->vRes[2]=1.0;

  /*optimsation things. Mainly for optimiseTLS.c and ignored here*/
  dimage->ofimage.weight=0;       /*weight to avoid subterranean*/
  dimage->ofimage.optSmoo=0;      /*optimise smoothing*/
  dimage->ofimage.optPsmoo=0;     /*optimise pre-smoothing*/
  dimage->ofimage.fixPS=0;        /*fix pulse length scaling*/
  dimage->ofimage.denWidth=0.0;   /*smoothing width when testing denoising accuracy*/
  dimage->ofimage.testAtten=1;    /*optimise to the attenuation corrected waveform*/
  dimage->ofimage.flWeight=0.0;   /*optimisation weight*/
  dimage->ofimage.doGBIC=0;       /*apply GBIC*/
  dimage->ofimage.errOut=NULL;
  dimage->ofimage.fixGap=0;       /*don't fix gap fraction*/

  /*canopy bounds*/
  strcpy(dimage->canNamen,"/home/sh563/data/bess/analysis/can_bounds/gBounds.dat");
  dimage->useCanBound=0;

  /*error settings*/
  dimage->hErrRes=0.15;
  dimage->gRes=0.05;  /*5% steps*/

  /*denoising defaults*/
  setDenoiseDefault(&dimage->denoise);
  strcpy(dimage->denoise.pNamen,dimage->pNamen);

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-tls",4)){
        checkArguments(1,i,argc,"-tls");
        strcpy(dimage->tlsNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-als",4)){
        checkArguments(1,i,argc,"-als");
        strcpy(dimage->alsNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-pFile",6)){
        checkArguments(1,i,argc,"-pFile");
        strcpy(dimage->pNamen,argv[++i]);
        strcpy(dimage->denoise.pNamen,dimage->pNamen);
      }else if(!strncasecmp(argv[i],"-noWaves",8)){
        dimage->waveOut=0;
      }else if(!strncasecmp(argv[i],"-waveBits",9)){
        dimage->waveBitsOut=1;
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-vRes",5)){
        checkArguments(3,i,argc,"-vRes");
        for(j=0;j<3;j++)dimage->vRes[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-smooth",7)){
        checkArguments(1,i,argc,"-smooth");
        dimage->denoise.sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-median",7)){
        checkArguments(1,i,argc,"-median");
        dimage->denoise.medLen=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-preSmooth",10)){
        checkArguments(1,i,argc,"-preSmooth");
        dimage->denoise.psWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-meanN",6)){
        checkArguments(1,i,argc,"-meanN");
        dimage->denoise.meanN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-tailThresh",11)){
        checkArguments(1,i,argc,"-tailThresh");
        dimage->denoise.tailThresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-thresh",7)){
        checkArguments(1,i,argc,"-thresh");
        dimage->denoise.thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minWidth",9)){
        checkArguments(1,i,argc,"-thresh");
        dimage->denoise.minWidth=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pScale",7)){
        checkArguments(1,i,argc,"-pScale");
        dimage->denoise.pScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxIter",8)){
        checkArguments(1,i,argc,"-maxIter");
        dimage->denoise.maxIter=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-tol",4)){
        checkArguments(1,i,argc,"-tol");
        dimage->denoise.deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gold",5)){
        dimage->denoise.deconMeth=0;
      }else if(!strncasecmp(argv[i],"-RL",3)){    /*Richardson-Lucy*/
        dimage->denoise.deconMeth=1;
      }else if(!strncasecmp(argv[i],"-varNoise",9)){    /*variable noise switch*/
        dimage->denoise.varNoise=1;
        dimage->denoise.medStats=1;
      }else if(!strncasecmp(argv[i],"-varScale",9)){
        checkArguments(1,i,argc,"-varScale");
        dimage->denoise.threshScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noTrack",8)){
        dimage->denoise.noiseTrack=0;
      }else if(!strncasecmp(argv[i],"-minGap",7)){
        checkArguments(1,i,argc,"-minGap");
        dimage->tlsPar.minGap=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-appRefl",8)){
        checkArguments(1,i,argc,"-appRefl");
        dimage->tlsPar.appRefl=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-GBIC",5)){    /*apply GBIC*/
        dimage->appGBIC=1;
      }else if(!strncasecmp(argv[i],"-pointOut",9)){
        dimage->pointOut=1;
      }else if(!strncasecmp(argv[i],"-balFlights",11)){    /*balance intensity between flights*/
        dimage->balFlights=1;
      }else if(!strncasecmp(argv[i],"-compBeams",10)){
        dimage->compBeams=1;
      }else if(!strncasecmp(argv[i],"-compVox",8)){
        dimage->compVox=1;
      }else if(!strncasecmp(argv[i],"-canBound",9)){
        dimage->useCanBound=1;
      }else if(!strncasecmp(argv[i],"-canFile",8)){
        checkArguments(1,i,argc,"-canFile");
        strcpy(dimage->canNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-gRes",5)){
        checkArguments(1,i,argc,"-gRes");
        dimage->gRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hErrRes",8)){
        checkArguments(1,i,argc,"-hErrRes");
        dimage->hErrRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-matchHard",10)){
        dimage->denoise.matchHard=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-tls name;       tls file list\n-als name;       als file list\n-vRes x y z;    voxel resolution in each axis\n-bounds minX minY minZ maxX maxY maxZ;   voxel bounds\n\nDenoising:\n-smooth w;     Gaussian smoothing\n-median w;     median filter width\n-preSmooth w;  pre-smoothing width\n-meanN n;      mean noise level\n-thresh n;     noise threshold\n-tailThresh n; tail noise threshold\n-minWidth w;   minimum feature width\n-varNoise;     use variable noise thresholds\n-varScale n;   factor to scale variable noise by\n-noTrack;      don't use noise tracking\n\nDeconvolution:\n-gold;         Gold's method\n-RL;           Richardson-Lucy deconvolution\n-pScale x;     pulse scale factor\n-maxIter n;    maximum number of iterations\n-tol tol;      deconvolution tolerance\n-pFile name;   pulse shape filename\n-matchHard;    look for hard targets\n\nTLS scaling\n-minGap gap;   minimum gap correction factor\n-appRefl n;    scale between 1064 and 1545 reflectance\n\n-waveBits;     output bits to make TLS waveforms\n-GBIC;         apply GBIC\n-balFlights;   balance intensity between flights, read from input list\n-pointOut;     output points\n-compBeams;    compare individual TLS and ALS beams\n-compVox;      compare TLS and ALS voxels\n-noWaves;      don't output waves\n-canBound;     use canopy bounds\n-canFile n;    canopy bound filename\n\nError statistics\n-gRes res;     cumulative gap fraction resolution\n-hErrRes res;  height error resolution\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*arg loop*/

  /*set to meanN if left blank*/
  if(dimage->denoise.tailThresh<0.0)dimage->denoise.tailThresh=dimage->denoise.meanN;

  return(dimage);
}/*readCommands*/

/*the end*/
/*####################################*/


