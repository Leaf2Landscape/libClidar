#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*###########################*/
/*# Reads raw SALCA binary  #*/
/*# and outputs a histogram #*/
/*# of return range         #*/
/*###########################*/


int main(int argc,char **argv)
{
  int i=0;          /*loop controls*/
  int numb=0;       /*number of zen zteps*/
  int nBins=0;      /*number of range bins*/
  int length=0;     /*total file length*/
  int nAz=0;        /*number of azimuth steps*/
  int *histogram=NULL;
  int *width=NULL;
  int *refl=NULL;
  float maxR=0;     /*maximum range*/
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int *,int *,int *);
  char inRoot[200],namen[200];   /*input filename root and filename*/
  char output[200];             /*output filename*/
  char writeDodgy=0;            /*output switch*/
  void makeHist(char *,int,int,int *,int *,int *,char,char *,int);
  void writeHist(int *,int *,int *,int,char *);


  /*the defaults are for the hedge test*/
  strcpy(inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  nAz=666;
  strcpy(output,"rangeHist.dat");
  maxR=30.0;
  numb=3200;   /*number of azimuth steps*/
  nBins=1200;
  writeDodgy=0;  /*don't write out dodgy by default*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(output,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-writeDodgy",5)){
        writeDodgy=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-output name;  output filename\n-nAz n;         number of azimuth steps (number of binary files)\n-maxR range;    maximum recorded range, metres\n-writeDodgy;    write out wavefrms with returns where we don't expect\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  /*set scan array dimensions*/
  nBins=(int)((150.0+maxR)/0.15);  /*150m of 1550 nm plus maxR of 1040 nm*/

  histogram=ialloc(nBins,"histogram",0);
  width=ialloc(nBins,"width",0);
  refl=ialloc(nBins,"refl",0);


  /*loop through azimuth files and write out points*/
  for(i=0;i<nAz;i++){ /*az loop*/
    sprintf(namen,"%s_%d.bin",inRoot,i);                  /*input filename*/
    fprintf(stdout,"%d of %d Reading %s\n",i,nAz,namen);  /*progress indicator*/
    data=readData(namen,&numb,&nBins,&length); /*read binary data into an array*/

    if(data)makeHist(data,numb,nBins,histogram,width,refl,writeDodgy,output,i);

    TIDY(data);   /*clean up arrays as we go along*/
  }/*az loop*/

  writeHist(histogram,width,refl,nBins,output);
  TIDY(histogram);
  TIDY(width);
  TIDY(refl);

  return(0);
}/*main*/


/*##########################################*/
/*write the histogram*/

void writeHist(int *histogram,int *width,int *refl,int nBins,char *output)
{
  int i=0,tot=0;
  float w=0,I=0;   /*width and intensity*/
  FILE *opoo=NULL;

  /*tot up*/
  for(i=0;i<nBins;i++)tot+=histogram[i];
  if((opoo=fopen(output,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",output);
    exit(1);
  }
  fprintf(opoo,"# 1 range, 2 number, 3 width, 4 intensity, 5 total\n");
  for(i=0;i<nBins;i++){
    if(histogram[i]){
      w=(float)width[i]/(float)histogram[i];
      I=(float)refl[i]/(float)histogram[i];
    }else  w=I=0.0;
    fprintf(opoo,"%f %d %f %f %d\n",(float)i*0.15,histogram[i],w,I,tot);
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",output);

  return;
}/*writeHist*/


/*##########################################*/
/*make a histogram*/

void makeHist(char *data,int numb,int nBins,int *histogram,int *width,int *refl,char writeDodgy,char *outRoot,int az)
{
  int i=0,j=0;  /*loop variables*/
  int s=0;      /*return start*/
  int place=0;  /*array place*/
  int range=0;  /*range index*/
  int minR=0,maxR=0;  /*dodgy bounds*/
  char thresh=0,max=0; /*intensity threshold and max value*/
  char dodgy=0;
  void writeWave(int,char *,int,char *,int);

  minR=(int)(140.0/0.15);
  maxR=(int)(147.5/0.15);

  thresh=-110;        /*chosen by RG, looks good*/
  for(i=0;i<numb;i++){ /*zenith loop*/
    dodgy=0;
    for(j=0;j<nBins;j++){     /*loop along waveform*/
      place=i*nBins+j;            /*array place*/
      if(data[place]>thresh){   /*signal above noise level*/
        max=data[place];          /*reset max*/
        s=j;                      /*mark start*/
        for(;j<nBins;j++){        /*step to end of feature*/
          place=i*nBins+j;        /*array place*/
          if(data[place]>max){  /*record maximum intensity and position*/
            max=data[place];      /*set max*/
            range=j;
          }                     /*max test*/
          if(data[place]<=thresh){
            histogram[range]++;
            width[range]+=j-s;
            refl[range]+=max+125;
            if(writeDodgy){   /*do we want to write this out*/
              if((range>=minR)&&(range<=maxR)){
                dodgy=1; /*we'll want to write this out*/
                fprintf(stdout,"%d %d %d\n",range,minR,maxR);
              }
            }
            break;
          }/*left feature*/
        }/*point end*/
      }/*point start*/
    }/*bin loop*/
    if(dodgy)writeWave(i,data,nBins,outRoot,az);
  }/*zenith loop*/

  return;
}/*makeHist*/


/*##########################################*/
/*output waveform with unexpected returns*/

void writeWave(int i,char *data,int nBins,char *outRoot,int az)
{
  int j=0;
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.az.%d.beam.%d.dat",outRoot,az,i);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  for(j=0;j<nBins;j++)fprintf(opoo,"%f %d\n",(float)j*0.15,data[i*nBins+j]);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }

  return;
}/*writeWave*/


/*##########################################*/
/*read data into array*/

char *readData(char *namen,int *numb,int *nBins,int *length)
{
  char *data=NULL;
  FILE *ipoo=NULL;

  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    return(NULL);
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
/*##########################################*/

