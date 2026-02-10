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
/*# Decimates TLS data       #*/
/*# Reads readRXP.cpp output #*/
/*# S Hancock, 2019          #*/
/*############################*/


/*#######################################*/
/*# Copyright 2014-2019, Steven Hancock #*/
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
  char outNamen[200];
  /*options*/
  double bounds[6];  /*area of interest, minX minY minZ, maxX, maxY maxZ*/
  float vRes[3];     /*voxel resolution in each axis*/
  int subSec;        /*number of voxels on a side to have as subsection*/
  float secLen;      /*length of a section*/
  /*how to split data*/
  uint32_t nXsec;    /*number of x sections*/
  uint32_t nYsec;    /*number of y sections*/
  uint32_t nZsec;    /*number of z sections*/
}control;


/*##################################################*/
/*voxel decimation structure*/

typedef struct{
  uint32_t nX;      /*number of x voxels*/
  uint32_t nY;      /*number of y voxels*/
  uint32_t nZ;      /*number of z voxels*/
  char *inHit;      /*hit or not within each voxel*/
  float res[3];     /*voxel resolution*/
  double bounds[6]; /*voxel bounds minX minY minZ maxX maxY maxZ*/
  FILE *opoo;       /*output file pointer*/
  char hasHit;      /*does this area have any gits at all?*/
}voxDec;


/*##################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  voxDec *vox=NULL;
  voxDec *allocateSimpleVox(control *);
  void decimateByVox(control *,voxDec *);


  /*read command line*/
  dimage=readCommands(argc,argv);
                              
  /*determine how to split the data and allocate voxel*/
  vox=allocateSimpleVox(dimage);

  /*do the decimation*/
  decimateByVox(dimage,vox);

  /*tidy up*/
  if(vox){
    TIDY(vox->inHit);
    TIDY(vox);
  }
  TIDY(dimage);
  return(0);
}/*main*/


/*##################################################*/
/*decimate by voxels*/

void decimateByVox(control *dimage,voxDec *vox)
{
  uint32_t i=0,j=0,k=0;
  void simpleVoxHits(voxDec *,char **,int,char *);

  /*loop over sub sections*/
  for(i=0;i<dimage->nXsec;i++){
    vox->bounds[0]=dimage->bounds[0]+(double)i*(double)dimage->secLen;
    vox->bounds[3]=(vox->bounds[0]+(double)dimage->secLen)>dimage->bounds[3]?\
                   dimage->bounds[3]:vox->bounds[0]+(double)dimage->secLen;

    /*check bounds are usable*/
    if((vox->bounds[3]-vox->bounds[0])<vox->res[0])continue;
    /*y section loop*/
    for(j=0;j<dimage->nYsec;j++){
      vox->bounds[1]=dimage->bounds[1]+(double)j*(double)dimage->secLen;
      vox->bounds[4]=(vox->bounds[1]+(double)dimage->secLen)>dimage->bounds[4]?\
                      dimage->bounds[4]:vox->bounds[1]+(double)dimage->secLen;
      /*check bounds are usable*/
      if((vox->bounds[4]-vox->bounds[1])<vox->res[1])continue;
      /*z section loop*/
      for(k=0;k<dimage->nZsec;k++){
        vox->bounds[2]=dimage->bounds[2]+(double)k*(double)dimage->secLen;
        vox->bounds[5]=(vox->bounds[2]+(double)dimage->secLen)>dimage->bounds[5]?\
                        dimage->bounds[5]:vox->bounds[2]+(double)dimage->secLen;
        /*check bounds are usable*/
        if((vox->bounds[5]-vox->bounds[2])<vox->res[2])continue;

        /*progress message*/
        fprintf(stdout,"Section %u of %u",1+k+j*dimage->nXsec+i*dimage->nXsec*dimage->nYsec,\
                                      dimage->nXsec*dimage->nYsec*dimage->nZsec);
        fprintf(stdout," Bounds %f %f %f %f %f %f\n",vox->bounds[0],vox->bounds[3],\
                      vox->bounds[1],vox->bounds[4],vox->bounds[2],vox->bounds[5]);

        /*perform voxelisation and write hits*/
        simpleVoxHits(vox,dimage->inList,dimage->nScans,dimage->outNamen);
      }/*z section loop*/
    }/*y section loop*/
  }/*x section loop*/

  /*close file*/
  if(vox->opoo){
    fclose(vox->opoo);
    vox->opoo=NULL;
  }

  /*delete file if empty*/
  if(vox->hasHit==0)remove(dimage->outNamen);
  else fprintf(stdout,"Written to %s\n",dimage->outNamen);

  return;
}/*decimateByVox*/


/*##################################################*/
/*find hits in boxel space*/

void simpleVoxHits(voxDec *vox,char **inList,int nScans,char *outNamen)
{
  int numb=0;
  uint8_t k=0;
  uint32_t j=0,tInd=0;
  int64_t vPlace=0,xBin=0,yBin=0,zBin=0;
  int64_t nX=0,nY=0,nZ=0;
  tlsScan *tempTLS=NULL;
  double x=0,y=0,z=0;
  double sinZen=0,cosZen=0;
  double sinAz=0,cosAz=0;

  /*if not already, open output*/
  if(vox->opoo==NULL){
    if((vox->opoo=fopen(outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",outNamen);
      exit(1);
    }
  }

  /*set all to zero*/
  for(vPlace=(uint64_t)vox->nX*(uint64_t)vox->nY*(uint64_t)vox->nZ-1;vPlace>=0;vPlace--)vox->inHit[vPlace]=0;
  

  /*to save lots of recasting in loops, recast now*/
  nX=(int64_t)vox->nX;
  nY=(int64_t)vox->nY;
  nZ=(int64_t)vox->nZ;

  /*loop over scans*/
  for(numb=0;numb<nScans;numb++){
    /*initial scan read*/
    readTLSpolarBinary(inList[numb],0,&tempTLS);

    /*loop over beams*/
    for(j=0;j<tempTLS->nBeams;j++){
      readTLSpolarBinary(inList[numb],j,&tempTLS);
      tInd=j-tempTLS->pOffset;   /*update index to account for buffered memory*/


      /*calculate trig, to save time in loop*/
      if(tempTLS->beam[tInd].nHits>0){
        sinZen=sin(tempTLS->beam[tInd].zen);
        cosZen=cos(tempTLS->beam[tInd].zen);
        sinAz=sin(tempTLS->beam[tInd].az);
        cosAz=cos(tempTLS->beam[tInd].az);
      }

      /*loop over hits in beam*/
      for(k=0;k<tempTLS->beam[tInd].nHits;k++){
        x=(double)tempTLS->beam[tInd].x+tempTLS->xOff+(double)tempTLS->beam[tInd].r[k]*sinZen*sinAz;
        y=(double)tempTLS->beam[tInd].y+tempTLS->yOff+(double)tempTLS->beam[tInd].r[k]*sinZen*cosAz;
        z=(double)tempTLS->beam[tInd].z+tempTLS->zOff+(double)tempTLS->beam[tInd].r[k]*cosZen;

        /*is it within bounds?*/
        if((x>=vox->bounds[0])&&(y>=vox->bounds[1])&&(z>=vox->bounds[2])&&\
           (x<vox->bounds[3])&&(y<vox->bounds[4])&&(z<vox->bounds[5])){
          /*if so, mark voxel as filled*/
          xBin=(int64_t)((x-vox->bounds[0])/(double)vox->res[0]);
          yBin=(int64_t)((y-vox->bounds[1])/(double)vox->res[1]);
          zBin=(int64_t)((z-vox->bounds[2])/(double)vox->res[2]);

          if((xBin>=0)&&(xBin<nX)&&(yBin>=0)&&(yBin<nY)&&(zBin>=0)&&(zBin<nZ)){
            vPlace=zBin*nX*nY+yBin*nX+xBin;
            /*mark point if this is the first*/
            if(vox->inHit[vPlace]==0){
              fprintf(vox->opoo,"%.3f %.3f %.3f %f\n",x,y,z,tempTLS->beam[tInd].refl[k]);
              vox->inHit[vPlace]=1;
              vox->hasHit=1;
            }
          }/*index check*/
        }/*bounds check*/
      }/*hit loop*/
    }/*beam loop*/
    tempTLS=tidyTLScan(tempTLS);
  }/*scan loop*/

  return;
}/*simpleVoxHits*/


/*##################################################*/
/*determine how to split the data and allocate voxel*/

voxDec *allocateSimpleVox(control *dimage)
{
  int i=0;
  voxDec *vox=NULL;
  float dist[3];

  /*allocate space and start copying parameters*/
  if(!(vox=(voxDec *)calloc(1,sizeof(voxDec)))){
    fprintf(stderr,"error decimated voxel structure allocation.\n");
    exit(1);
  }
  vox->opoo=NULL;
  for(i=0;i<3;i++)vox->res[i]=dimage->vRes[i];
  vox->hasHit=0;

  /*determine size of each section*/
  for(i=0;i<3;i++)dist[i]=(float)(dimage->bounds[i+3]-dimage->bounds[i]);
  if(dimage->subSec<=0){
    for(i=0;i<3;i++)if(dist[i]>dimage->secLen)dimage->secLen=dist[i];
  }else{
    dimage->secLen=(float)dimage->subSec*dimage->vRes[0];
  }

  /*determine number of sections*/
  dimage->nXsec=(uint32_t)(dist[0]/dimage->secLen+1.0);
  dimage->nYsec=(uint32_t)(dist[1]/dimage->secLen+1.0);
  dimage->nZsec=(uint32_t)(dist[2]/dimage->secLen+1.0);

  for(i=0;i<3;i++){
    dist[i]=(float)(dimage->bounds[3+i]-dimage->bounds[i])<dimage->secLen?(float)(dimage->bounds[3+i]-dimage->bounds[i]):dimage->secLen;
  }
  vox->nX=(uint32_t)(dist[0]/vox->res[0]+0.5);
  vox->nY=(uint32_t)(dist[1]/vox->res[1]+0.5);
  vox->nZ=(uint32_t)(dist[2]/vox->res[2]+0.5);

  /*allocate voxels*/
  vox->inHit=challoc((uint64_t)vox->nX*(uint64_t)vox->nY*(uint64_t)vox->nZ,"inHit",0);

  return(vox);
}/*allocateSimpleVox*/


/*##################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;

  /*allocate space*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control structure allocation.\n");
    exit(1);
  }

  /*set defaults*/
  dimage->inList=NULL;
  dimage->nScans=0;
  strcpy(dimage->outNamen,"teast.pts");
  dimage->subSec=-1;

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
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->vRes[0]=dimage->vRes[1]=dimage->vRes[2]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-voxPer",7)){
        checkArguments(1,i,argc,"-voxPer");
        dimage->subSec=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n################################\nProgram to decimate TLS data using voxels\n################################\n\n-input name;     input TLS filename. .bin file expected\n-inList list;    input file list for multiple bin files\n-output name;    output filename\n-bounds minX minY minZ maxX maxY maxZ:    voxel bounds\n-res res;       vocel resolution\n-voxPer n;      voxels per segment, along a side (total will be cubed)\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry voxDecimate -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  /*check inputs are valid*/
  if((dimage->inList==NULL)||(dimage->nScans==0)){
    fprintf(stderr,"No input specified\n");
    exit(1);
  }

  return(dimage);
}/*readCommands*/

/*the end*/
/*##################################################*/

