#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*##########################*/
/*# Reads raw SALCA binary #*/
/*# and outputs a histogram#*/
/*##########################*/

int sOffset;   /*wavestart offset*/


int main(int argc,char **argv)
{
  int i=0,band=0;   /*loop controls*/
  int numb=0;       /*number of zen zteps*/
  int nBins=0;      /*number of range bins*/
  int length=0;     /*total file length*/
  int nAz=0;        /*number of azimuth steps*/
  int start[2],end[2]; /*array start and end bounds*/
  int nHist=0;
  float maxR=0;     /*maximum range*/
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int *,int *,int *);
  char inRoot[200],namen[200];  /*input filename root and filename*/
  char output[200];             /*output filename*/
  char subBin=0;   /*sub bin range resolution flag*/
  void totHistogram(char *,long int *,int, int,int);
  void writeHisto(long int *,int,char *);
  long int *hist=NULL;

  /*the defaults are for the hedge test*/
  strcpy(inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  nAz=666;
  strcpy(output,"salcaHist.hist");
  maxR=30.0;
  numb=3200;   /*number of azimuth steps*/
  subBin=0;    /*use peak intensity bin only*/
  nHist=257;

  hist=lIalloc(nHist,"histogram",0);

  sOffset=7;   /*1.05 m*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",8)){
        checkArguments(1,i,argc,"-output");
        strcpy(output,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-output name;  output filename\n-nAz n;         number of azimuth steps (number of binary files)\n-maxR range;    maximum recorded range, metres\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  /*set scan array dimensions*/
  nBins=(int)((150.0+maxR)/0.15);  /*150m of 1550 nm plus maxR of 1040 nm*/

  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(maxR/0.15);
  start[1]=1000+sOffset;   /*7 to avoid outgoing pulses*/
  end[1]=nBins;

  /*loop through azimuth files and write out points*/
  for(i=0;i<nAz;i++){ /*az loop*/
    sprintf(namen,"%s_%d.bin",inRoot,i);                  /*input filename*/
    fprintf(stdout,"%d of %d Reading %s\n",i,nAz,namen);  /*progress indicator*/
    data=readData(namen,&numb,&nBins,&length); /*read binary data into an array*/
    if(data){ /*does the file exist?*/
      totHistogram(data,hist,nBins,numb,nHist);
    }/*does file exist?*/
    TIDY(data);   /*clean up arrays as we go along*/
  }/*az loop*/

  /*output histogram*/
  writeHisto(hist,nHist,output);

  TIDY(hist);
  return(0);
}/*main*/


/*##########################################*/
/*add up histogram*/

void totHistogram(char *data,long int *hist,int nBins, int numb,int nHist)
{
  int i=0,bin=0;
  int place=0;
  int ind=0;

  for(i=0;i<numb;i++){/*zen loop*/
    for(bin=0;bin<nBins;bin++){ /*bin loop*/
      place=i*nBins+bin;
      ind=(int)(data[place]+128);
      if((ind>=0)&&(ind<nHist))hist[ind]++;
      else                     fprintf(stderr,"Gone beyond\n");
    }/*bin loop*/
  }/*zen loop*/

  return;
}/*totHistogram*/


/*##########################################*/
/*find outgoing pulse and check saturation*/

int findStart(int start,int end,char *satTest,char *data,int offset)
{
  int i=0,b=0,e=0;  /*loop control and bounds*/
  int place=0;      /*array index*/
  int waveStart=0;  /*outgoing pulse bin*/
  char max=0;       /*max intensity*/
  char satThresh=0; /*saturation threshold*/
  char thresh=0;    /*noise threshold*/

  waveStart=-1; /*start-sOffset;*/

  satThresh=127;
  thresh=-110;
  max=-125;

  *satTest=0;  /*not saturated by default*/
  b=start-50;   /*4.5m, from histograms*/
  if(b<0)b=0;  /*truncate at 0*/
  e=start+sOffset;
  if(e>end)e=end;

  max=-125;
  for(i=b;i<e;i++){  /*loop around where we think the pulse might be*/
    place=offset+i;
    if((data[place]>thresh)&&(data[place]>max)){ /*outgoing peak*/
      max=data[place];
      waveStart=i;
    }                    /*outgoing peak test*/
    if((max>-125)&&(data[place]<=thresh))break;
  }/*range loop*/

  for(i=start;i<end;i++){  /*loop through full waveform to test for saturation*/
    place=offset+i;
    if(data[place]>=satThresh){   /*saturated*/
      *satTest=1;
      break;
    }  
  }/*saturation test loop*/

  return(waveStart);
}/*findStart*/


/*##########################################*/
/*read data into array*/

char *readData(char *namen,int *numb,int *nBins,int *length)
{
  char *data=NULL;
  FILE *ipoo=NULL;

  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    return(NULL);   /*no file, return NULL pointer*/
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


/*##########################################*/
/*write histogram*/

void writeHisto(long int *hist,int nHist,char *output)
{
  int i=0;
  long int total=0;
  FILE *opoo=NULL;

  if((opoo=fopen(output,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",output);
    exit(1);
  }

  for(i=0;i<nHist;i++)total+=(long int)hist[i];
  for(i=0;i<nHist;i++)fprintf(opoo,"%d %f\n",i-128,(float)hist[i]/(float)total);

  fprintf(stdout,"Written to %s\n",output);

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  return;
}/*writeHisto*/


/*the end*/
/*##########################################*/

