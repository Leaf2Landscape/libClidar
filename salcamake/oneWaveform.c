#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*##########################*/
/*# Reads raw SALCA binary #*/
/*# and outputs ascii      #*/
/*# for one waveform       #*/
/*##########################*/

int main(int argc,char **argv)
{
  int numb=0,nBins=0,length=0;
  int i=0,j=0;
  float maxR=0;
  char *data=NULL;
  char *readData(char *,int *,int,int *);
  char namen[200],outNamen[200];
  void writeWave(char *,int,int,int,char *,int);

  strcpy(namen,"/home/server/users/bakgrp3/nsh103/data/salca_play/hedge_plot34_f1_1936_166.bin");
  strcpy(outNamen,"oneWave.dat");
  j=-1;  /*do all by default*/
  maxR=30.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",7)){
        checkArguments(1,i,argc,"-input");
        strcpy(namen,argv[++i]); 
      }else if(!strncasecmp(argv[i],"-zenInd",7)){
        checkArguments(1,i,argc,"-zenInd");
        j=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;   input filename\n-output name;  output filename\n-zenInd n;     zenith index. Set negative for all\n-maxR r;     maximum range\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  nBins=(int)((150.0+maxR)/0.15);  /*150m of 1550 nm plus maxR of 1040 nm*/
  data=readData(namen,&numb,nBins,&length);

  writeWave(data,numb,nBins,j,outNamen,length);

  TIDY(data);

  return(0);
}/*main*/


/*##########################################*/
/*write out a single wave*/

void writeWave(char *data,int numb,int nBins,int j,char *outNamen,int length)
{
  int bin=0,place=0,start=0;
  float zen=0;
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }
  fprintf(opoo,"# 1 range, 2 intensity\n");

  if(j>=0){   /*print out a single waveform*/
    for(bin=-15;bin<nBins;bin++){
      place=j*nBins+bin;
      if(place>0)fprintf(opoo,"%g %d\n",(float)bin*0.15,data[place]);
    }
  }else{     /*print out all waveforms*/
    zen=-91.68;   /*lowest zenith*/
    for(bin=0;bin<length;bin++){
      if((bin%nBins)==0){
        zen+=0.057295779;  /*zen step in degs*/
        start=bin;
      }
      fprintf(opoo,"%g %d %g %g\n",(float)bin*0.15,data[bin],zen,(float)(bin-start)*0.15);
    }
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);
  return;
}/*writeWave*/


/*##########################################*/
/*read data into array*/

char *readData(char *namen,int *numb,int nBins,int *length)
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
  *numb=3200;    /*(int)((*length)/(*nBins));*/
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
/*########################################*/

