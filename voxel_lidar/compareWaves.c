#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"


/*#######################*/
/*# Matched NERC-ARSF   #*/
/*# las 1.3 data to TLS #*/
/*# S Hancock, 2014     #*/
/*#######################*/


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



/*############################################*/
/*control structure*/

typedef struct{
  char listNamen[200];  /*ALS input file list*/
  char tlsNamen[200];   /*TLS input file list*/
  char alignNamen[200]; /*TLS alignment data*/
  char outRoot[100];    /*output filename root*/
  double bounds[6];     /*lon and lat bounds to use*/
  uint64_t maxWaves;    /*maximum number to output*/
  char pointOut;        /*output points or not*/
  uint64_t counter;     /*output number of waves output*/
  float beamRad;        /*ALS beam radius*/
  float tlsTanDiv;      /*tan() of tls beam divergence*/
  float tlsRad;         /*tls beam start diameter*/
  int nTLSwaves;        /*number of TLS waves*/
  int binOff;           /*number of prebins to look for geolocation*/
  /*scan centres*/
  char useCent;         /*use scan centres switch*/
  char centNamen[200];  /*scan centre filename*/
  /*tls parameters*/
  float appRefl;        /*apparent reflectance scaling*/
  float minGap;         /*minimum trustable gap fraction*/
  char calTLS;          /*calibrate TLS as we go*/

  /*options*/
  char waveOut;         /*output waveforms*/
  char heightOut;       /*output top range*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/

  /*denoising parameters*/
  denPar denoise;      /*structure of denoising variables*/
}control;


/*############################################*/
/*TLS data*/

typedef struct{
  uint64_t numb;
  double *x;
  double *y;
  double *z;
  float *rad;      /*beam radius*/
  /*map TLS points*/
  float mapRes;   /*map resolution*/
  int mNx;        /*map x dimension*/
  int mNy;        /*map y dimension*/
  int mNz;        /*map z dimension*/
  int **map;      /*map from pixels to point coords. May need to be 64 bit at some point*/
  int *nIn;       /*number of points in this pixel*/
}tlsData;


/*############################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  listLas *lasList=NULL;
  tlsData *tls=NULL;
  tlsData *readTLS(control *);
  void readPulse(denPar *);
  void compareWaves(listLas *,tlsData *,control *);

  /*read the command line*/
  dimage=readCommands(argc,argv);

  /*read ALS input filename list*/
  fprintf(stdout,"Reading ALS headers\n");
  lasList=readLasList(dimage->listNamen);

  /*read TLS data*/
  fprintf(stdout,"Reading TLS data\n");
  tls=readTLS(dimage);

  /*read ALS pulse*/
  if(dimage->denoise.deconMeth>=0)readPulse(&dimage->denoise);

  /*match up ALS and TLS points*/
  fprintf(stdout,"Matching waves\n");
  compareWaves(lasList,tls,dimage);

  /*tidy up arrays*/
  if(tls){/*tls*/
    TIDY(tls->x);
    TIDY(tls->y);
    TIDY(tls->z);
    TIDY(tls->rad);
    TTIDY((void **)tls->map,tls->mNx*tls->mNy*tls->mNz);
    TIDY(tls->nIn);
    TIDY(tls);
  }/*tls*/
  if(dimage){/*control*/
    TTIDY((void **)dimage->denoise.pulse,2);
    TIDY(dimage);
  }/*control*/
  tidyListLas(lasList);
  return(0);
}/*main*/


/*############################################*/
/*match up waveforms*/

void compareWaves(listLas *lasList,tlsData *tls,control *dimage)
{
  int i=0,j=0;
  int nVox=0,*voxList=NULL;
  uint64_t lastWave=0;
  double xCent=0,yCent=0,zCent=0;
  double x0=0,y0=0,z0=0;
  double discRange=0;
  double tempRes[3];
  float **tlsWave=NULL;
  float **intersectBeam(int,int *,tlsData *,lasFile *,double,double,double,double,double,double,control *);
  lasFile *lasIn=NULL;
  void outputWaves(unsigned char *,float **,int32_t,char *,uint64_t,control *,double,double,double,lasFile *,float *,int,double,uint16_t);
  /*denoise and deconvolve*/
  float *denoised=NULL;
  /*height stuff*/
  void writeHeights(float *,float **,int,float,float *,uint64_t,FILE *,int);
  char hNamen[200];
  FILE *heightFile=NULL;

  for(i=0;i<3;i++)tempRes[i]=(double)tls->mapRes;


  if(dimage->heightOut){  /*open height file if needed*/
    sprintf(hNamen,"%s.height",dimage->outRoot);
    if((heightFile=fopen(hNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",hNamen);
      exit(1);
    }
    fprintf(heightFile,"# 1 tlsRange, 2 waveRange, 3 discRange, 4 counter\n");
  }/*open height file if needed*/

  dimage->counter=0;  /*reset counter*/
  for(i=0;i<lasList->nFiles;i++){  /*file loop*/
    lasIn=readLasHead(lasList->nameList[i],dimage->pBuffSize);

    /*check bounds*/
    if((lasIn->minB[0]<=dimage->bounds[1])&&(lasIn->maxB[0]>=dimage->bounds[0])\
      &&(lasIn->minB[1]<=dimage->bounds[3])&&(lasIn->maxB[1]>=dimage->bounds[2])){
      lastWave=3999;
      for(j=0;j<lasIn->nPoints;j++){  /*point loop*/
        readLasPoint(lasIn,j);
        if(lasIn->waveMap!=lastWave){
          setCoords(&xCent,&yCent,&zCent,lasIn);
          if((xCent>=dimage->bounds[0])&&(xCent<=dimage->bounds[1])&&\
             (yCent>=dimage->bounds[2])&&(yCent<=dimage->bounds[3])){
            lasIn->wave=readLasWave(lasIn->waveMap,lasIn->waveLen,lasIn->ipoo,lasIn->waveStart);

            /*deconvolve*/
            denoised=processWave(lasIn->wave,(int)lasIn->waveLen,&(dimage->denoise),lasIn->gbic[lasIn->psID]);

            /*work out what map pixels this intersects*/
            /*REMOVE BEFORE USE - need to include the beam width at start and end*/
            binPosition(&x0,&y0,&z0,0,xCent,yCent,zCent,lasIn->time,lasIn->grad);
            discRange=sqrt((x0-xCent)*(x0-xCent)+(y0-yCent)*(y0-yCent)+(z0-zCent)*(z0-zCent));
            voxList=beamVoxels(lasIn->grad,x0,y0,z0,&(dimage->bounds[0]),\
                          &(tempRes[0]),tls->mNx,tls->mNy,tls->mNz,&nVox,0.165,NULL,-1.0);


            /*now loop through the TLS data to see what point lies within this cone*/
            tlsWave=intersectBeam(nVox,voxList,tls,lasIn,x0,y0,z0,xCent,yCent,zCent,dimage);

            /*output the two waves*/
            if(dimage->waveOut)outputWaves(lasIn->wave,tlsWave,lasIn->waveLen,dimage->outRoot,dimage->counter,dimage,xCent,yCent,zCent,lasIn,denoised,dimage->binOff,discRange,lasIn->refl);
            if(dimage->heightOut)writeHeights(denoised,tlsWave,lasIn->waveLen,lasIn->time,lasIn->grad,dimage->counter,heightFile,dimage->binOff);
            dimage->counter++;

            TIDY(lasIn->wave);
            TTIDY((void **)tlsWave,dimage->nTLSwaves);
            TIDY(voxList)
            TIDY(denoised);
            lastWave=lasIn->waveMap;
          }/*bounds check*/
        }/*once per waveform*/
      }/*point loop*/
    }/*file extent check*/
    tidyLasFile(lasIn);
  }/*file loop*/

  if(heightFile){
    fprintf(stdout,"Written to %s\n",hNamen);
    fclose(heightFile);
    heightFile=NULL;
  }

  return;
}/*compareWaves*/


/*############################################*/
/*write out top ranges*/

void writeHeights(float *denoised,float **tlsWave,int waveLen,float time,float *grad,uint64_t counter,FILE *heightFile,int binOff)
{
  int i=0;
  float tlsTol=0,waveTol=0;
  float discRange=-1000.0;
  float waveRange=-1000.0;
  float tlsRange=-1000.0; 
  char tlsOut[20],waveOut[20],discOut[20];

  tlsTol=0.00001;  /*a floating point tolerance*/
  waveTol=1.0;

  /*signal is already denoised and noise tracked. First return is range*/
  for(i=0-binOff;i<0;i++){
    if((tlsRange<0.0)&&(tlsWave[1][i+binOff]>tlsTol)){
      tlsRange=(float)i*0.15;
      break;
    }
  }

  /*first bin of either is start*/
  for(i=0;i<waveLen;i++){
    if((waveRange<0.0)&&(denoised[i]>waveTol))waveRange=(float)i*0.15;
    if((tlsRange<0.0)&&(tlsWave[1][i+binOff]>tlsTol))tlsRange=(float)i*0.15;
  }
  discRange=time*sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);

  /*put in blanks for gnuplot*/
  if(tlsRange>-100.0)sprintf(tlsOut,"%f",tlsRange);
  else               sprintf(tlsOut,"?");
  if(waveRange>-100.0)sprintf(waveOut,"%f",waveRange);
  else               sprintf(waveOut,"?");
  if(discRange>-100.0)sprintf(discOut,"%f",discRange);
  else                sprintf(discOut,"?");
 
  fprintf(heightFile,"%s %s %s %ld\n",tlsOut,waveOut,discOut,(long int)counter);
  /*fflush(heightFile);*/

  return;
}/*writeHeights*/


/*############################################*/
/*write waveforms from different sources*/

void outputWaves(unsigned char *wave,float **tlsWave,int32_t waveLen,char *outRoot,uint64_t nOut,control *dimage,double lon,double lat,double elev,lasFile *lasIn,float *denoised,int binOff,double discRange,uint16_t discRefl)
{
  int i=0,j=0;
  char namen[200];
  double xWave=0,yWave=0,zWave=0;
  FILE *opoo=NULL;

  sprintf(namen,"%s.%llu.dat",outRoot,nOut);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  fprintf(opoo,"# 1 range, 2 ALS, 3 TLS count, 4 TLS visible cover, 5 TLS bin cover, 6 TLS summit, 7 denoised, 8 lon, 9 lat, 10 elev\n");
  fprintf(opoo,"# lon %f lat %f elev %f\n",lon,lat,elev);
  fprintf(opoo,"# discRange %f %d\n",discRange,(int)discRefl);
  for(i=0-binOff;i<0;i++){
    fprintf(opoo,"%f ?",(float)i*0.15);
    for(j=0;j<dimage->nTLSwaves;j++)fprintf(opoo," %f",tlsWave[j][i+binOff]);
    binPosition(&xWave,&yWave,&zWave,i,lon,lat,elev,lasIn->time,lasIn->grad);
    fprintf(opoo," ? %f %f %f",xWave,yWave,zWave);
    fprintf(opoo,"\n");
  }
  for(i=0;i<(int)waveLen;i++){
    fprintf(opoo,"%f %d",(float)i*0.15,wave[i]);
    for(j=0;j<dimage->nTLSwaves;j++)fprintf(opoo," %f",tlsWave[j][i+binOff]);
    binPosition(&xWave,&yWave,&zWave,i,lon,lat,elev,lasIn->time,lasIn->grad);
    fprintf(opoo," %f %f %f %f",denoised[i],xWave,yWave,zWave);
    fprintf(opoo,"\n");
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",namen);
  return;
}/*outputWaves*/


/*############################################*/
/*find TLS points within ALS beam*/

float **intersectBeam(int nPix,int *pixList,tlsData *tls,lasFile *lasIn,double x0,double y0,double z0,double xCent,double yCent,double zCent,control *dimage)
{
  int i=0,j=0;
  int pixInd=0;
  int xInd=0,yInd=0,rPlace=0;   /*range image indices*/
  int xIcent=0,xStart=0,xEnd=0; /*range image indices*/
  int yIcent=0,yStart=0,yEnd=0; /*range image indices*/
  int rNx=0,rNy=0;              /*range image dimensions*/
  int bin=0;
  float **tlsWave=NULL;        /*multiple TLS waves. 0 point count, 1 projected cover*/
  float sepSq=0,minSepSq=0;
  float zen=0,az=0;
  float alsRadSq=0;       /*ALS beam radius squared*/
  float **rImage=NULL;    /*range image, a stack nBins long*/
  float rimRes=0;         /*range image resolution*/
  float rSepSq=0,maxRsepSq=0;   /*range image separation*/
  double *vect=NULL;
  void rotateX(double *,double);
  void rotateY(double *,double);
  void rotateZ(double *,double);
  char doneIt=0;    /*a marker*/

  /*For outputting points*/
  double xWave=0,yWave=0,zWave=0;
  FILE *opoo=NULL;
  char namen[200];
  if(dimage->pointOut){
    sprintf(namen,"%s.%llu.pts",dimage->outRoot,dimage->counter);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
  }/*output points or not*/

  /*allocate arrays*/
  alsRadSq=dimage->beamRad*dimage->beamRad;
  vect=dalloc(3,"vector",0);
  tlsWave=fFalloc(dimage->nTLSwaves,"tls waveform",0);   /*this is multiple waves now*/
  for(i=0;i<dimage->nTLSwaves;i++)tlsWave[i]=falloc(lasIn->waveLen+dimage->binOff,"tls waveform",i+1);
  /*range image bits*/
  rimRes=0.001;   /*1 mm*/
  rNx=rNy=(int)(2.0*dimage->beamRad/rimRes);
  rImage=fFalloc(lasIn->waveLen+dimage->binOff,"",0);
  for(j=0;j<lasIn->waveLen+dimage->binOff;j++){
    rImage[j]=falloc(rNx*rNy,"range image",j+1);
    for(i=rNx*rNy-1;i>=0;i--)rImage[j][i]=-1000.0;
  }

  minSepSq=dimage->beamRad*dimage->beamRad;   /*beam radius squared*/
  zen=(float)atan2(sqrt((double)lasIn->grad[0]*(double)lasIn->grad[0]+\
      (double)lasIn->grad[1]*(double)lasIn->grad[1]),(double)lasIn->grad[2]);
  az=(float)atan2((double)lasIn->grad[0],(double)lasIn->grad[1]);

  for(i=0;i<nPix;i++){   /*TLS voxel loop*/
    for(j=0;j<tls->nIn[pixList[i]];j++){ /*points within voxel loop*/
      /*read pixel and load into array. Translate to put beam at origin*/
      pixInd=tls->map[pixList[i]][j];
      vect[0]=tls->x[pixInd]-x0;
      vect[1]=tls->y[pixInd]-y0;
      vect[2]=tls->z[pixInd]-z0;

      /*rotate to x-y plane*/
      rotateZ(vect,(double)(-1.0*az));
      rotateX(vect,(double)(-1.0*zen));

      /*calculate distance from beam*/
      sepSq=(float)(vect[0]*vect[0]+vect[1]*vect[1]);
      if(sepSq<=minSepSq){
        bin=(int)(vect[2]/0.15+0.5)+dimage->binOff;
        if((bin>=0)&&(bin<(lasIn->waveLen+dimage->binOff))){  /*within beam*/
          tlsWave[0][bin]+=tls->rad[pixInd]*tls->rad[pixInd]*dimage->appRefl/alsRadSq;
          tlsWave[3][bin]+=1.0;
          /*write point cloud*/
          if(dimage->pointOut)fprintf(opoo,"%f %f %f ? ? ?\n",tls->x[pixInd],tls->y[pixInd],tls->z[pixInd]);

          /*range image*/
          xIcent=(int)((vect[0]+dimage->beamRad)/rimRes);
          yIcent=(int)((vect[1]+dimage->beamRad)/rimRes);
          xStart=xIcent-(int)(tls->rad[pixInd]/rimRes);
          xEnd=xIcent+(int)(tls->rad[pixInd]/rimRes);
          yStart=yIcent-(int)(tls->rad[pixInd]/rimRes);
          yEnd=yIcent+(int)(tls->rad[pixInd]/rimRes);
          maxRsepSq=tls->rad[pixInd]*tls->rad[pixInd];

          if(xStart<0)xStart=0;      /*enforce bounds*/
          if(xEnd>=rNx)xEnd=rNx-1;
          if(yStart<0)yStart=0;
          if(yEnd>=rNy)yEnd=rNy-1;   /*enforce bounds*/

          for(xInd=xStart;xInd<=xEnd;xInd++){
            for(yInd=yStart;yInd<=yEnd;yInd++){
              rSepSq=(float)((xInd-xStart)*(xInd-xStart)+(yInd-yStart)*(yInd-yStart))*rimRes*rimRes;
              if(rSepSq<=maxRsepSq){
                rPlace=yInd*rNx+xInd;
                if((rImage[bin][rPlace]<0.0)||(vect[2]<rImage[bin][rPlace]))rImage[bin][rPlace]=vect[2];
              }/*check within point*/
            }/*loop around point*/
          }/*loop around point*/
        }/*within beam*/
      }/*separation check*/
    }/*tls point loop*/
  }/*map voxel loop*/

  /*turn range images into a waveforms*/
  for(i=rNx*rNy-1;i>=0;i--){
    doneIt=0;
    for(j=0;j<lasIn->waveLen+dimage->binOff;j++){
      if(rImage[j][i]>=0.0){
        if(doneIt==0){
          tlsWave[1][j]+=dimage->appRefl*(rimRes*rimRes)/(M_PI*dimage->beamRad*dimage->beamRad);
          doneIt=1;
        }/*first hit only*/
        tlsWave[2][j]+=dimage->appRefl*(rimRes*rimRes)/(M_PI*dimage->beamRad*dimage->beamRad);  /*hits per bin*/
      }
    }/*range image loop*/
  }/*range image bin loop*/

  /*output ALS waveform as individual points*/
  if(dimage->pointOut){
    for(i=0;i<lasIn->waveLen;i++){
      if(lasIn->wave[i]>16){
        binPosition(&xWave,&yWave,&zWave,i,xCent,yCent,zCent,lasIn->time,lasIn->grad);
        fprintf(opoo,"? ? ? %f %f %f\n",xWave,yWave,zWave);
      }
    }
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Points written to %s\n",namen);
  }/*output individual points*/

  TIDY(vect);
  TTIDY((void **)rImage,lasIn->waveLen+dimage->binOff);
  return(tlsWave);
}/*intersectBeam*/


/*############################################*/
/*read TLS file ASCII (.pts), whichever way*/

tlsData *readTLS(control *dimage)
{
  int i=0,j=0;
  int xBin=0,yBin=0,zBin=0,mapPlace=0;    /*map array coordinates*/
  double xCent=0,yCent=0,zCent=0;
  double range=0;
  double *markDo(int,double *,double);
  double setRange(double,double,double,double *);
  float pointSize(double,uint16_t,float,float,float,float,float,float);
  tlsData *tls=NULL;
  listLas *tlsList=NULL;
  lasFile *lasIn=NULL;

  /*allocate space for structures*/
  if(!(tls=(tlsData *)calloc(1,sizeof(tlsData)))){
    fprintf(stderr,"error in tls structure.\n");
    exit(1);
  }
  tls->mapRes=1.0;     /*hard wired in for now*/

  /*read TLS input filename list*/
  tlsList=readLasList(dimage->tlsNamen);

  /*read alignment and scan centres*/
  TTIDY((void **)tlsList->align,tlsList->nFiles);  /*as they're allocated in readLasList*/
  tlsList->align=readCoordList(tlsList->nFiles,tlsList->nameList,dimage->alignNamen);
  if(dimage->useCent){
    TTIDY((void **)tlsList->scanCent,tlsList->nFiles);  /*as they're allocated in readLasList*/
    tlsList->scanCent=readCoordList(tlsList->nFiles,tlsList->nameList,dimage->centNamen);
  }

  /*allocate space for map*/
  tls->mNx=(int)((dimage->bounds[1]-dimage->bounds[0])/tls->mapRes)+1;
  tls->mNy=(int)((dimage->bounds[3]-dimage->bounds[2])/tls->mapRes)+1;
  tls->mNz=(int)((dimage->bounds[5]-dimage->bounds[4])/tls->mapRes)+1;
  tls->nIn=ialloc(tls->mNx*tls->mNy*tls->mNz,"nIn",0);
  tls->map=iIalloc(tls->mNx*tls->mNy*tls->mNz,"map",0);

  /*loop throgh files*/
  tls->numb=0;
  for(i=0;i<tlsList->nFiles;i++){
    lasIn=readLasHead(tlsList->nameList[i],dimage->pBuffSize);
    /*translate alignment to global coords*/
    for(j=0;j<2;j++){
      lasIn->minB[j]+=tlsList->align[i][j];
      lasIn->maxB[j]+=tlsList->align[i][j];
    }
    /*check file extent overlaps*/
    if((lasIn->minB[0]<=dimage->bounds[1])&&(lasIn->maxB[0]>=dimage->bounds[0])\
      &&(lasIn->minB[1]<=dimage->bounds[3])&&(lasIn->maxB[1]>=dimage->bounds[2])){
      for(j=0;j<lasIn->nPoints;j++){  /*point loop*/
        readLasPoint(lasIn,j);
        setCoords(&xCent,&yCent,&zCent,lasIn);
        range=setRange(xCent,yCent,zCent,tlsList->scanCent[i]);  /*calculate range before realigning*/
        xCent+=tlsList->align[i][0];
        yCent+=tlsList->align[i][1];
        zCent+=tlsList->align[i][2];

        if((xCent>=dimage->bounds[0])&&(xCent<=dimage->bounds[1])&&\
           (yCent>=dimage->bounds[2])&&(yCent<=dimage->bounds[3])&&\
           (zCent>=dimage->bounds[4])&&(zCent<=dimage->bounds[5])){
          tls->x=markDo((int)tls->numb,tls->x,xCent);
          tls->y=markDo((int)tls->numb,tls->y,yCent);
          tls->z=markDo((int)tls->numb,tls->z,zCent);
          tls->rad=markFloat((int)tls->numb,tls->rad,pointSize(range,lasIn->refl,dimage->tlsTanDiv,dimage->tlsRad,14000.0,33577.0,1.0,1.0));

          /*map the point*/
          xBin=(int)((xCent-dimage->bounds[0])/tls->mapRes);
          yBin=(int)((yCent-dimage->bounds[2])/tls->mapRes);
          zBin=(int)((zCent-dimage->bounds[4])/tls->mapRes);
          mapPlace=zBin*tls->mNx*tls->mNy+yBin*tls->mNx+xBin;
          tls->map[mapPlace]=markInt(tls->nIn[mapPlace],tls->map[mapPlace],(int)tls->numb);
          tls->nIn[mapPlace]++;
          tls->numb++;
        }/*bounds check*/
      }/*point loop*/
    }/*file extent check*/
    if(lasIn->ipoo){
      fclose(lasIn->ipoo);
      lasIn->ipoo=NULL;
    }
  }/*file loop*/

  /*tidy up*/
  tidyListLas(tlsList);
  tidyLasFile(lasIn);
  return(tls);
}/*readTLS*/


/*############################################*/
/*calculate range from scan centre*/

double setRange(double x,double y,double z,double *scanCent)
{
  double range=0;
  double dx=0,dy=0,dz=0;

  dx=x-scanCent[0];
  dy=y-scanCent[1];
  dz=z-scanCent[2];

  range=sqrt(dx*dx+dy*dy+dz*dz);
  return(range);
}/*setRange*/


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

  /*options*/
  dimage->waveOut=1;
  dimage->heightOut=0;

  /*defaults*/
  dimage->nTLSwaves=4;  /*0 target profile, 1 visible projected area, 2 projected area, 3 point count*/
  dimage->maxWaves=10000000;
  strcpy(dimage->tlsNamen,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/lists/luton.all.list");
  strcpy(dimage->listNamen,"/home/sh563/data/bess/lists/las1.3.LT.dat");
  dimage->bounds[0]=507211.0;
  dimage->bounds[1]=507212.0;
  dimage->bounds[2]=225374.2;
  dimage->bounds[3]=225394.2;
  dimage->bounds[4]=120.0;
  dimage->bounds[5]=150.0;
  strcpy(dimage->outRoot,"teast");
  dimage->pointOut=0;      /*don't output individual points*/
  dimage->counter=0;
  dimage->beamRad=0.33/2.0;  /*33 cm footprints*/
  dimage->useCent=1;
  strcpy(dimage->centNamen,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/lists/rieglCentres.dat");
  strcpy(dimage->alignNamen,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/lists/luton.align.dat");
  dimage->tlsTanDiv=sin(0.35/(2.0*1000.0))/cos(0.35/(2.0*1000.0));  /*tan half of 0.35 mrad*/
  dimage->tlsRad=0.035;        /*7 mm diameter beam start*/
  dimage->appRefl=0.01;
  dimage->minGap=0.0;
  dimage->calTLS=0;
  dimage->binOff=0;
  dimage->pBuffSize=(uint64_t)200000000;

  /*denoising defaults*/
  setDenoiseDefault(&dimage->denoise);


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-lasList",8)){
        checkArguments(1,i,argc,"-lasList");
        strcpy(dimage->listNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-tlsList",8)){
        checkArguments(1,i,argc,"-tlsList");
        strcpy(dimage->tlsNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pointOut",9)){
        dimage->pointOut=1;
      }else if(!strncasecmp(argv[i],"-scanCent",9)){
        checkArguments(1,i,argc,"-scanCent");
        dimage->useCent=1;
        strcpy(dimage->centNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-align",6)){
        checkArguments(1,i,argc,"-align");
        strcpy(dimage->alignNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-appRefl",8)){
        checkArguments(1,i,argc,"-appRefl");
        dimage->appRefl=atof(argv[++i]);
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
        dimage->denoise.deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-RL",3)){
        dimage->denoise.deconMeth=1;
      }else if(!strncasecmp(argv[i],"-varNoise",9)){
        dimage->denoise.varNoise=1;
        dimage->denoise.medStats=1;
      }else if(!strncasecmp(argv[i],"-varScale",9)){
        checkArguments(1,i,argc,"-varScale");
        dimage->denoise.threshScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-binOff",7)){
        checkArguments(1,i,argc,"-binOff");
        dimage->binOff=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noWaves",8)){
        dimage->waveOut=0;
      }else if(!strncasecmp(argv[i],"-height",7)){
        dimage->heightOut=1;
      }else if(!strncasecmp(argv[i],"-fitGauss",9)){
        dimage->denoise.fitGauss=1;
      }else if(!strncasecmp(argv[i],"-gWidth",7)){
        checkArguments(1,i,argc,"-gWidth");
        dimage->denoise.gWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gaussPulse",11)){
        dimage->denoise.gaussPulse=1;
      }else if(!strncasecmp(argv[i],"-findHard",9)){
        dimage->denoise.matchHard=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\nInput output\n-lasList name;        list of input lasfiles\n-tlsList name;        list of tls filenames\n-outRoot root;        output filename root\n\n-align file;          alignment filename\n-scanCent name;       read TLS scan centres from a file\nSwitches and limits\n-bounds minX maxX minY maxY minZ maxZ;   bounds to search\n-pointOut;            output individual points\n-noWaves;      don't output waveforms\n-height;       output top ranges\n-binOff n;     number of pixels to look before\n\nALS denoising\n-thresh n;   noise threshold\n-meanN n;    mean noise level\n-minWidth n; minimum feature width to accept\n-tailThresh t; trailing edge threshold\n-smooth s;\n-preSmooth s;\n-median n;     median smoothing width\n-varNoise;     variable noise threshold\n-varScale s;   variable noise scale factor\n\nALS denconvolution\n-maxIter n;  maximum deconvolution iterations\n-tol tol;    deconvolution tolerance\n-pScale p;   pulse width scaling\n-gold;       Gold's method\n-RL;         Richarson-Lucy method\n\nGaussian fitting\n-fitGauss;   fit Gaussians\n-gWidth g;   smoothing width to count Gaussians\n-gaussPulse;     convert pulse to a Gaussian before decon\n\nTLS parameters\n-appRefl r;  apparent reflectance at 1064nm\n-minGap g;   minimum trustable gap\n-calTLS;     calibrate TLS as we go\n-findHard;         find hard targets\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*arg loop*/

  if(dimage->denoise.tailThresh<0.0)dimage->denoise.tailThresh=dimage->denoise.meanN;  /*set to meanN if left blank*/

  return(dimage);
}/*readCommands*/

/*the end*/
/*############################################*/

