#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "mpfit.h"
#include "libLasProcess.h"


/*############################*/
/*# Reads ASCII waveform,    #*/
/*# denoises and deconvolves #*/
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
  char inNamen[200]; /*input filename*/
  char outNamen[200];/*output filename*/
  int nWaves;        /*number of different waves*/
  int nBins;         /*number of waveform bins*/
  char gaussFit;     /*fit Gaussian switch*/
  float gSmooth;     /*smooth when choosing turning points*/
  denPar denoise;      /*structure of denoising variables*/
}control;


/*##################################################*/
/*main*/

int main(int argc, char **argv)
{
  int i=0;
  int nGauss=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  float **data=NULL;
  float **readData(char *,int,int *);
  float *decon=NULL,*temp=NULL;
  void readPulse(denPar *);
  float *fitMultiGauss(float *,float *,int,float,int *);
  unsigned char *wave=NULL;
  void writeData(control *,float **,float *);


  /*read the command line*/
  dimage=readCommands(argc,argv);

  /*read input*/
  data=readData(dimage->inNamen,dimage->nWaves,&dimage->nBins);

  /*read pulse shape*/
  readPulse(&dimage->denoise);

  /*denoise*/
  wave=uchalloc(dimage->nBins,"unsigned char",0);
  for(i=0;i<dimage->nBins;i++)wave[i]=(unsigned char)data[1][i];

  /*process wave*/
  /*decon=processWave(wave,dimage->nBins,dimage->meanN,dimage->thresh,dimage->minWidth,pulse,\
                    dimage->pBins,dimage->res,dimage->maxIter,dimage->minChange,\
                    dimage->deconMeth,dimage->sWidth,dimage->psWidth,dimage->medLen,\
                    dimage->tailThresh,dimage->varNoise);*/
  decon=processWave(wave,dimage->nBins,&(dimage->denoise),1.0);

  TTIDY((void **)dimage->denoise.pulse,2);

  /*function fit if we are doing it*/
  if(dimage->gaussFit){
    temp=fitMultiGauss(&(data[0][0]),decon,dimage->nBins,dimage->gSmooth,&nGauss);
    TIDY(decon);
    decon=temp;
    temp=NULL;
  }

  /*output results*/
  writeData(dimage,data,decon);

  /*tidy up arrays*/
  TIDY(decon);
  TTIDY((void **)data,dimage->nWaves);
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##################################################*/
/*write out results*/

void writeData(control *dimage,float **data,float *decon)
{
  int i=0,j=0;
  FILE *opoo=NULL;

  if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  }

  fprintf(opoo,"# 1 range, 2 deconvolved, 3 smoothed, 4 denoised, 5 raw, 6 TLS area, 7 TLS vis, 8+ bits\n");
  for(i=0;i<dimage->nBins;i++){
    fprintf(opoo,"%f %f %f %f",data[0][i],decon[i],-1.0,-1.0);
    for(j=1;j<dimage->nWaves;j++)fprintf(opoo," %f",data[j][i]);
    fprintf(opoo,"\n");
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*writeData*/


/*##################################################*/
/*read data*/

float **readData(char *inNamen,int nWaves,int *nBins)
{
  int i=0,j=0;
  float **data=NULL;
  char line[200];
  char temp[6][100];
  FILE *ipoo=NULL;

  /*allocate space and open*/
  data=fFalloc(nWaves,"data",0);
  if((ipoo=fopen(inNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",inNamen);
    exit(1);
  }

  *nBins=0;
  /*count number of bins*/
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))(*nBins)++;
  }/*line counting*/

  for(j=0;j<nWaves;j++)data[j]=falloc(*nBins,"",j+1);

  /*rewind*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s %s %s",temp[0],temp[1],temp[2],temp[3],temp[4],temp[5])==6){
        if(i>*nBins){
          fprintf(stderr,"Error\n");
          exit(1);
        }
        for(j=0;j<nWaves;j++)data[j][i]=atof(&(temp[j][0]));
        i++;
      }
    }/*comment check*/
  }/*data reading*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readData*/


/*##################################################*/
/*read the command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);


  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  strcpy(dimage->inNamen,"/home/sh563/data/bess/analysis/penetration/wavecomp/waveComp.1/lutonPlot.1.1457.dat");
  strcpy(dimage->outNamen,"testDenoise.dat");
  strcpy(dimage->denoise.pNamen,"/home/sh563/data/bess/leica_shape/leicaPulse.dat");
  dimage->nWaves=6; /*0, range, 1 ALS, 2 TLS count, 3 TLS visible cover, 4 TLS bin cover*/

  /*denoising defaults*/
  setDenoiseDefault(&dimage->denoise);
  dimage->denoise.minWidth=3;  /*minimum width above noise to accept*/
  dimage->denoise.deconMeth=0; /*Gold's method*/
  dimage->denoise.deChang=0.00001;

  dimage->gaussFit=0;
  dimage->gSmooth=0.0;

  for(i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-thresh",7)){
        checkArguments(1,i,argc,"-thresh");
        dimage->denoise.thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-meanN",6)){
        checkArguments(1,i,argc,"-meanN");
        dimage->denoise.meanN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-tailThresh",11)){
        checkArguments(1,i,argc,"-tailThresh");
        dimage->denoise.tailThresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minWidth",9)){
        checkArguments(1,i,argc,"-minWidth");
        dimage->denoise.minWidth=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxIter",8)){
        checkArguments(1,i,argc,"-maxIter");
        dimage->denoise.maxIter=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-smooth",7)){
        checkArguments(1,i,argc,"-smooth");
        dimage->denoise.sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-preSmooth",10)){
        checkArguments(1,i,argc,"-preSmooth");
        dimage->denoise.psWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-median",7)){
        checkArguments(1,i,argc,"-median");
        dimage->denoise.medLen=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rlDecon",8)){
        dimage->denoise.deconMeth=1;
      }else if(!strncasecmp(argv[i],"-pScale",6)){
        checkArguments(1,i,argc,"-pScale");
        dimage->denoise.pScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minChange",10)){
        checkArguments(1,i,argc,"-minChange");
        dimage->denoise.deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gauss",6)){
        dimage->gaussFit=1;
        dimage->denoise.deconMeth=-1;
      }else if(!strncasecmp(argv[i],"-gSmooth",8)){
        checkArguments(1,i,argc,"-gSmooth");
        dimage->gSmooth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;   waveform filename\n-outut name;   output filename\n-thresh n;     noise threshold\n-meanN n;      mean noise level\n-tailThresh t; trailing edge threshold\n-minWidth n;   minimum feature width in bins\n-maxIter n;    maximum number of iterations\n-minChange t;  threshold between interations\n-smooth s;     smoothing width\n-preSmooth s;  pre-smoothing width\n-median n;     median smoothing length\n-rlDecon;      Richardson-Lucy deconvolution instead of Gold's\n-pScale f;     factor to scale pulse length by\n-gauss;        fit Gaussians\n-gSmooth sig;     smoothing width when finding turning points\n\n");
        exit(0);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  if(dimage->denoise.tailThresh<0.0)dimage->denoise.tailThresh=dimage->denoise.meanN;  /*set to meanN if left blank*/

  return(dimage);
}/*readCommands*/


/*the end*/
/*##################################################*/

