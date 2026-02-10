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
/*# Reads Riegl VZ-400 #*/
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
  char **inList;
  int nFiles;
  char outNamen[200];

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
  control *dimage=NULL;
  control *readCommands(int,char **);
  FILE *opoo=NULL;


  /*read options and set defaults*/
  dimage=readCommands(argc,argv);
  if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",dimage->outNamen);
    exit(1);
  }
  fprintf(opoo,"# 1 x, 2 y, 3 z, 4 refl, 5 hitN, 6 nHits, 7 zen, 8 az, 9 range\n");

  /*dump point cloud to file*/
  for(i=0;i<dimage->nFiles;i++){
    fprintf(stdout,"Doing %d of %d %s\n",i+1,dimage->nFiles,dimage->inList[i]);
    writeTLSpointFromBin(dimage->inList[i],dimage->bounds,opoo);
  }
  fprintf(stdout,"Point cloud written to %s\n",dimage->outNamen);

  /*tidy up*/
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage);
  }
  return(0);
}/*main*/


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

  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"",0);
  dimage->inList[i]=challoc(200,"",1);
  strcpy(dimage->inList[i],"/Users/stevenhancock/data/teast/mondah_tls/ScanPos009.bin");
  strcpy(dimage->outNamen,"test.pts");
  dimage->pBuffSize=(uint64_t)200000000;

  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=-1000000000.0;  /*blank value*/
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=1000000000.0;  /*blank value*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inList[0],argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;     input TLS binary file\n-inList name;    list of input TLS binary files\n-output name;    output filename\n-bounds minX minY minZ maxX maxY maxZ;     limit bounds\n\n");
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

