#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*##########################*/
/*# Simple program to      #*/
/*# investigate saturation #*/
/*##########################*/

int main(int argc,char **argv)
{
  int i=0,j=0,nAz=0;
  int numb=0,nBins=0,length=0;
  int start[2],end[2];
  int nSat=0;       /*number of saturated waves*/
  float maxR=0;
  char outSat=0;
  char satTest=0;
  char testSat(char *,int,int,int);
  char outRoot[200],outNamen[200];
  char inRoot[200],namen[200];
  char **saveWave(char *,int,int,int,int,int,char **);
  char **saturated=NULL;
  char *readData(char *,int *,int *,int *);
  char *data=NULL;
  void printWaves(char **,int,int,char *);

  nAz=666;
  strcpy(inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  strcpy(outRoot,"satTest");
  numb=3200;
  maxR=30.0;
  outSat=1;     /*output saturated waveforms only*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-unsat",6)){
        outSat=0;    /*output unsaturated waveforms only*/
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-outRoot root;  output filename root\n-nAz n;         number of azimuth steps (number of binary files)\n-maxR range;    maximum recorded range, metres\n-unsat;    otput unsaturated waveforms only\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/



  nBins=(int)((150.0+maxR)/0.15);

  /*band start and end bins for sampling*/
  start[0]=0;
  end[0]=(int)(maxR/0.15);
  start[1]=1000;
  end[1]=nBins;


  /*loop through azimuth files and write out points*/
  for(i=0;i<nAz;i++){ /*az loop*/
    sprintf(namen,"%s_%d.bin",inRoot,i);                  /*input filename*/
    fprintf(stdout,"%d of %d\n",i,nAz);  /*progress indicator*/
    data=readData(namen,&numb,&nBins,&length); /*read binary data into an array*/

    nSat=0;
    for(j=0;j<3200;j++){/*zen loop*/
      satTest=testSat(data,j*nBins,start[0],end[1]);  /*loop through both bands*/
      if((outSat&&satTest)||((!outSat)&&(!satTest))){  /*two different modes*/
        saturated=saveWave(data,j*nBins,start[0],end[1],nSat,nBins,saturated);
        nSat++;
      }
    }/*zen loop*/
    if(nSat){
      sprintf(outNamen,"%s.az.%d.dat",outRoot,i);
      printWaves(saturated,nSat,nBins,outNamen);
      TTIDY((void **)saturated,nSat);
    }
    TIDY(data);
  }/*az loop*/

  return(0);
}/*main*/


/*##########################################################*/
/*save saturated waves in an array*/

char **saveWave(char *data,int offset,int s,int e,int nSat,int nBins,char **old)
{
  int i=0,j=0;
  char **saturated=NULL;

  saturated=chChalloc(nSat+1,"saturated",0);
  for(i=0;i<=nSat;i++){
    saturated[i]=challoc(e-s,"saturated",i+1);
    if(i<nSat)for(j=0;j<nBins;j++)saturated[i][j]=old[i][j];
    else      for(j=0;j<nBins;j++)saturated[i][j]=data[offset+j];
  }
  if(nSat)TTIDY((void **)old,nSat);

  return(saturated);
}/*saveWave*/


/*##########################################################*/
/*print out a waveform*/

void printWaves(char **saturated,int nSat,int nBins,char *outNamen)
{
  int i=0,j=0;
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }

  for(i=0;i<nBins;i++){  /*range loop*/
    fprintf(opoo,"%g",(float)i*0.15);
    for(j=0;j<nSat;j++){ /*zen loop*/
      fprintf(opoo," %d",saturated[j][i]);
    }/*range loop*/
    fprintf(opoo,"\n");
  }/*zen loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Saturated wave %s\n",outNamen);

  return;
}/*printWaves*/


/*##########################################################*/
/*test for saturation*/

char testSat(char *data,int offset,int s,int e)
{
  int i=0,place=0;
  char satThresh=0;
  char satTest=0;

  satThresh=120;
  satTest=0;

  for(i=s;i<e;i++){
    place=offset+i;
    if(data[place]>=satThresh){
      satTest=1;
      break;
    }
  }
  return(satTest);
}/*testSat*/


/*##########################################*/
/*read data into array*/

char *readData(char *namen,int *numb,int *nBins,int *length)
{
  char *data=NULL;
  FILE *ipoo=NULL;

  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*determine the file length*/
  if(fseek(ipoo,(long)0,SEEK_END)){ /*jump to file end*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  *length=ftell(ipoo);

  /*these are now set on the command line*/
  /* *nBins=1200; */   /*Should be *length/numb;, but doesn't seem to be*/
  /* *numb=3200; */    /*(int)((*length)/(*nBins));*/
  if(((*nBins)*(*numb))>(*length)){
    fprintf(stderr,"File size mismatch\n");
    exit(1);
  }
  data=challoc(*length,"data",0);

  /*now we know hoe long, read the file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to file start*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  if(fread(&(data[0]),sizeof(char),*length,ipoo)!=*length){
    fprintf(stderr,"error reading data\n");
    exit(1);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readData*/


/*the end*/
/*##########################################################*/


