#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libTLSread.h"


/*######################*/
/*# Reads libRAY hips  €*/
/*# and writes points  #*/
/*# S Hancock, 2017    #*/
/*######################*/


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


/*#################################################*/
/*control structure*/

typedef struct{
  char inNamen[200];
  char outNamen[200];
  char equalAng;       /*projection flag*/
}control;


/*#################################################*/
/*hips header*/

typedef struct{
  double origin[3];   /*scan origin*/
  int nX;
  int nY;
  float zRes;
  int lastByte;
}hipsHead;


/*#################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  void convetWoodgateData(control *);

  /*read options and set defaults*/
  dimage=readCommands(argc,argv);

  /*read data and convert to binary fornat*/
  convetWoodgateData(dimage);

  /*tidy up*/
  TIDY(dimage);
  return(0);
}/*the end*/


/*#################################################*/
/*read and write data*/

void convetWoodgateData(control *dimage)
{
  hipsHead *header=NULL;
  hipsHead *readHipsHead(control *);
  void readWriteWoodgate(hipsHead *,control *);

  /*read header*/
  header=readHipsHead(dimage);

  /*read flat binary data*/
  readWriteWoodgate(header,dimage);

  return;
}/*convetWoodgateData*/


/*#################################################*/
/*read and write data*/

void readWriteWoodgate(hipsHead *header,control *dimage)
{
  int i=0,j=0;
  uint32_t nShots=0;
  float x=0,y=0;
  float *buff=NULL;
  float r=0,refl=0;
  double zen=0,az=0,radius=0;
  char nHits=0;
  FILE *ipoo=NULL,*opoo=NULL;

  /*open files*/
  if((ipoo=fopen(dimage->inNamen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",dimage->inNamen);
    exit(1);
  }
  if((opoo=fopen(dimage->outNamen,"wb"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  }

  /*skip header*/
  fprintf(stdout,"last byte %d n %d %d\n",header->lastByte,header->nX,header->nY);
  if(fseek(ipoo,(long)header->lastByte,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  buff=falloc((uint64_t)header->nX*(uint64_t)header->nY,"data",0);
  if(fread(&(buff[0]),sizeof(float),header->nX*header->nY,ipoo)!=header->nX*header->nY){
    fprintf(stderr,"error reading data from %s\n",dimage->inNamen);
    exit(1);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }


  /*loop over points*/
  refl=1.0;
  nShots=0;
  for(i=0;i<header->nX;i++){
    x=(float)(i-header->nX/2);
    for(j=0;j<header->nY;j++){
      y=(float)(header->nY/2-j);
      az=atan2((double)y,(double)x)*180/M_PI;
      radius=(double)sqrt(x*x+y*y)/((double)(header->nX/2));
      if(dimage->equalAng)zen=radius*90.0;
      else                zen=2.0*atan2(radius,1.0+sqrt(1.0-radius*radius))*180.0/M_PI;
      r=buff[j*header->nX+i];

//if(r>0.0)fprintf(stdout,"%f %f %f r %f %f %f\n",r*sin(zen*M_PI/180.0)*cos(az*M_PI/180.0),r*sin(zen*M_PI/180.0)*sin(az*M_PI/180.0),r*cos(zen*M_PI/180.0),r,zen,az);


      if(zen<=90.0){   /*write the point*/
        if(fwrite(&zen,sizeof(double),1,opoo)!=1){
          printf("integer header not written\n");
          exit(1);
        }
        if(fwrite(&az,sizeof(double),1,opoo)!=1){
          printf("integer header not written\n");
          exit(1);
        }
        if(fwrite(header->origin,sizeof(double),3,opoo)!=3){
          printf("integer header not written\n");
          exit(1);
        }
        if(fwrite(&nShots,sizeof(uint32_t),1,opoo)!=1){
          printf("integer header not written\n");
          exit(1);
        }

        if(r>0.0)nHits=1;
        else     nHits=0;

        if(fwrite(&nHits,sizeof(char),1,opoo)!=1){
          printf("integer header not written\n");
          exit(1);
        }

        if(nHits>0){
          if(fwrite(&r,sizeof(float),1,opoo)!=1){
            printf("integer header not written\n");
            exit(1);
          }
          if(fwrite(&refl,sizeof(float),1,opoo)!=1){
            printf("integer header not written\n");
            exit(1);
          }
        }
        nShots++;
      }
    }
  }

  /*end with number of points*/
  if(fwrite(&nShots,sizeof(uint32_t),1,opoo)!=1){
    printf("integer header not written\n");
    exit(1);
  }




  /*close up*/
  TIDY(buff);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*readWriteWoodgate*/


/*#################################################*/
/*read hips header*/

hipsHead *readHipsHead(control *dimage)
{
  int i=0,j=0;
  hipsHead *header=NULL;
  char line[400];
  char *token=NULL;
  FILE *ipoo=NULL;

  /*allocate and open*/
  if(!(header=(hipsHead *)calloc(1,sizeof(hipsHead)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }
  if((ipoo=fopen(dimage->inNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",dimage->inNamen);
    exit(1);
  }

  /*read lines*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(i==0){  /*check first line is HIPS*/
      if(strncasecmp(line,"HIPS",4)){
        fprintf(stderr,"Not hips format\n");
        exit(1);
      }
    }else if(!strncasecmp(line,".",1)){   /*end of header*/
      header->lastByte=ftell(ipoo)+0;
      break;
    }else if(i==4)header->nX=atoi(line);
    else if(i==5)header->nY=atoi(line);
    else if(i==6){
      if(atoi(line)!=32){
        fprintf(stderr,"not 32 bit, %s\n",line);
        exit(1);
      }
    }else if(!strncasecmp(line,"fishStart",9)){
      j=-1;
      token=strtok(line," ");
      while(token!=NULL) {
        if(!strncasecmp(token,"-from",5))j=0;
        else{
          if((j>=0)&&(j<3))header->origin[j]=atof(token);
          if(j>=0)j++;
        }
        token=strtok(NULL," ");
      }
    }
    i++;
  }/*data reading*/

  header->zRes=((float)header->nX/2.0)/(M_PI/2.0);

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(header);
}/*readHipsHead*/


/*#################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  dimage->equalAng=1;
  strcpy(dimage->outNamen,"test.bin");

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-fishEye",8)){
        checkArguments(1,i,argc,"-fishEye");
        dimage->equalAng=0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;     input TLS csv file\n-output name;    output filename\n-fishEye;        fish eye projection\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }
  return(dimage);
}/*readCommands*/

/* the end*/
/*###########################################################*/
