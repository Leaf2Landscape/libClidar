#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"

/*######################*/
/*# Reads Riegl VZ-400 #*/
/*# S Hancock, 2014    #*/
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
  char listNamen[200];
  char outNamen[200];
  char centNamen[200];   /*scan centre file*/

  char useCent;     /*use centre switch*/
  char outPoints;   /*write point switch*/
  double bounds[6]; /*bounds to output*/
  char outPolar;    /*write polar coords switch*/
  char fitPolar;    /*fit lines to polar coordinates*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
}control;


/*#################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  listLas *lasList=NULL;
  lasFile **lasIn=NULL;          /*las files*/
  control *dimage=NULL;
  control *readCommands(int,char **);
  void writePoints(lasFile *,control *,listLas *);
  void writePolar(lasFile *,control *,listLas *);
  void fitScanLines(lasFile *,control *,listLas *);


  /*read options and set defaults*/
  dimage=readCommands(argc,argv);

  /*read input filename list*/
  lasList=readLasList(dimage->listNamen);

  /*read scan centres*/
  if(dimage->useCent)lasList->scanCent=readCoordList(lasList->nFiles,lasList->nameList,dimage->centNamen);

  /*allocate space for las files*/
  lasIn=lfalloc(lasList->nFiles);

  /*open and read file headers*/
  for(i=0;i<lasList->nFiles;i++)lasIn[i]=readLasHead(lasList->nameList[i],dimage->pBuffSize);

  /*output point cloud*/
  if(dimage->outPoints)writePoints(lasIn[0],dimage,lasList);

  /*output polar coords*/
  if(dimage->outPolar)writePolar(lasIn[0],dimage,lasList);

  /*fit to polar scan lines*/
  if(dimage->fitPolar)fitScanLines(lasIn[0],dimage,lasList);

  /*tidy up*/
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
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*#################################################*/
/*fit to polar coordinate scan lines*/

void fitScanLines(lasFile *lasIn,control *dimage,listLas *lasList)
{
  int i=0,j=0;
  int nIn=0;
  double x=0,y=0,z=0;
  double dx=0,dy=0,dz=0;
  double zen=0,az=0;
  double *tempZen=NULL,*tempAz=NULL;
  double lastZen=0;
  /*line fitting*/
  int gsl_fit_linear(const double *,const size_t,const double *,const size_t,const size_t,double *, double *,double *,double *,double *,double *);
  double c=0,m=0;   /*line parameters*/
  double cov_00=0,cov_01=0,cov_11=0,sumsq=0;   /*fit metrics*/
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.polFit",dimage->outNamen);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  /*reset counters*/
  lastZen=-400.0;
  nIn=0;
  /*loop through files*/
  for(i=0;i<lasList->nFiles;i++){
    fprintf(stdout,"File %d of %d\n",i+1,lasList->nFiles);
    for(j=0;j<lasIn[i].nPoints;j++){
      readLasPoint(&(lasIn[i]),(uint32_t)j);
      setCoords(&x,&y,&z,&(lasIn[i]));

      /*calculate angle*/
      dx=x-lasList->scanCent[i][0];
      dy=y-lasList->scanCent[i][1];
      dz=z-lasList->scanCent[i][2];

      zen=atan2(sqrt(dx*dx+dy*dy),dz);
      az=atan2(dy,dx);
//fprintf(opoo,"%d %d %f %f %f %f %f\n",i,j,zen,az,x,y,z);
      if(zen>lastZen){/*new scan line*/
        gsl_fit_linear(tempAz,1,tempZen,1,nIn,&c,&m,&cov_00,&cov_01,&cov_11,&sumsq);
        TIDY(tempAz);
        TIDY(tempZen);
        //fprintf(opoo,"%d %f %f %d\n",i,m,c,nIn);
        nIn=0;
      }/*new scan line check*/
      tempZen=markDo(nIn,tempZen,zen);
      tempAz=markDo(nIn,tempAz,az);
      nIn++;
      lastZen=zen;
    }/*point loop*/
    if(nIn>0){
      gsl_fit_linear(tempAz,1,tempZen,1,nIn,&c,&m,&cov_00,&cov_01,&cov_11,&sumsq);
      TIDY(tempAz);
      TIDY(tempZen);
      //fprintf(opoo,"%d %f %f %d\n",i,m,c,nIn);
      nIn=0;
    }
  }/*file loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Fits written to %s\n",namen);
  TIDY(tempAz);
  TIDY(tempZen);
  return;
}/*fitScanLines*/


/*#################################################*/
/*write polar coords*/

void writePolar(lasFile *lasIn,control *dimage,listLas *lasList)
{
  int i=0,j=0;
  double x=0,y=0,z=0;
  double dx=0,dy=0,dz=0;
  double zen=0,az=0,r=0;
  char namen[200];
  FILE *opoo=NULL;

  for(i=0;i<lasList->nFiles;i++){
    sprintf(namen,"%s.%d.pol",dimage->outNamen,i);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
    fprintf(opoo,"# 1 zen, 2 az, 3 range, 4 refl, 5 ret numb\n");
    for(j=0;j<lasIn[i].nPoints;j++){
      readLasPoint(&(lasIn[i]),(uint32_t)j);
      setCoords(&x,&y,&z,&(lasIn[i]));

      /*calculate angle*/
      dx=x-lasList->scanCent[i][0];
      dy=y-lasList->scanCent[i][1];
      dz=z-lasList->scanCent[i][2];

      zen=atan2(sqrt(dx*dx+dy*dy),dz);
      az=atan2(dy,dx);
      r=sqrt(dx*dx+dy*dy+dz*dz);

      fprintf(opoo,"%f %f %f %d %d\n",zen*180.0/M_PI,az*180.0/M_PI,r,(int)lasIn[i].refl,lasIn[i].field.retNumb);
    }/*point loop*/
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",namen);
  }/*file loop*/

  return;
}/*writePolar*/


/*#################################################*/
/*write point cloud*/

void writePoints(lasFile *lasIn,control *dimage,listLas *lasList)
{
  int i=0,j=0;
  double x=0,y=0,z=0;
  double range=0;
  double setRange(double,double,double,double *);
  FILE *opoo=NULL;

  if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  }

  if(dimage->bounds[0]>=9900000.0){  /*it is blank and needs resetting*/
     dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=100000000.0;
     dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=-100000000.0;
     for(i=0;i<lasList->nFiles;i++){
       for(j=0;j<3;j++){
         if(lasIn[i].minB[j]<dimage->bounds[j])dimage->bounds[j]=lasIn[i].minB[j];
         if(lasIn[i].maxB[j]>dimage->bounds[j+3])dimage->bounds[j+3]=lasIn[i].maxB[j];
       }
     }
  }

  range=-1.0;   /*blank value if it is not to be used*/

  fprintf(opoo,"# 1 x, 2 y, 3 z, 4 int, 5 range\n");
  for(i=0;i<lasList->nFiles;i++){
    fprintf(stdout,"Reading file %d of %d\n",i+1,lasList->nFiles);
    for(j=0;j<lasIn[i].nPoints;j++){
      readLasPoint(&(lasIn[i]),j);
      setCoords(&x,&y,&z,&(lasIn[i]));
      if((x>=dimage->bounds[0])&&(x<=dimage->bounds[3])&&(y>=dimage->bounds[1])&&\
         (y<=dimage->bounds[4])&&(z>=dimage->bounds[2])&&(z<=dimage->bounds[5])){

        if(dimage->useCent)range=setRange(x,y,z,lasList->scanCent[i]);
        fprintf(opoo,"%f %f %f %d %f\n",x,y,z,(int)lasIn[i].refl,range);
      }
    }/*point loop*/
  }/*file loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*writePoints*/


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


/*#################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  strcpy(dimage->listNamen,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/lists/luton.003.list");
  strcpy(dimage->outNamen,"test.pts");
  dimage->outPoints=1;
  dimage->outPolar=0;
  dimage->fitPolar=0;
  dimage->pBuffSize=(uint64_t)200000000;

  dimage->bounds[0]=10000000.0;  /*blank value*/
  dimage->useCent=0;             /*do not use scan centres*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->listNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-cent",5)){
        checkArguments(1,i,argc,"-cent");
        strcpy(dimage->centNamen,argv[++i]);
        dimage->useCent=1;
      }else if(!strncasecmp(argv[i],"-polar",6)){
        dimage->outPolar=1;
      }else if(!strncasecmp(argv[i],"-fitPolar",9)){
        dimage->fitPolar=1;
      }else if(!strncasecmp(argv[i],"-noPoints",9)){
        dimage->outPoints=0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;     list of input lasfiles\n-output name;    output filename\n-cent name;      read image centres from file\n-polar;          output polar coords\n-fitPolar;       fit to scan lines\n-noPoints;       don't write points\n-bounds minX minY minZ maxX maxY maxZ;     limit bounds\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }
  return(dimage);
}/*readCommands*/

/*the end*/
/*#################################################*/

