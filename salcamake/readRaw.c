#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*##########################*/
/*# Reads raw SALCA binary #*/
/*# and outputs ascii      #*/
/*##########################*/

int main(int argc,char **argv)
{
  int i=0;
  int numb=0;       /*number of zen zteps*/
  int nBins=0;      /*number of range bins*/
  int length=0;     /*total file length*/
  float maxR=0;
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int,int,int *);
  char namen[200],outRoot[100];
  char mode=0;      /*output mode*/
  void outputData(char *,int,int,char *,int,char);

  mode=0;
  numb=4422;   /*hard wired into SALCA*/
  strcpy(namen,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936y/hedge_plot34_f1_1936_246.bin");
  strcpy(outRoot,"salcaTest");
  maxR=30.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(namen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-single",7)){ /*output a single file*/
        mode=1;
      }else if(!strncasecmp(argv[i],"-multi",6)){ /*output only multiple returns*/
        mode=2;
      }else if(!strncasecmp(argv[i],"-early",6)){ /*output only early returns*/
        mode=3;
      }else if(!strncasecmp(argv[i],"-maxR",5)){ /*max range*/
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;   binary input filename\n-outRoot root;   output filename root\n-single;\n-multi;\n-early;\n-maxR range;   maximum recorded range\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  nBins=(int)((150.0+maxR)/0.15);

  /*read binary data into an array*/
  data=readData(namen,numb,nBins,&length);

  outputData(data,numb,nBins,outRoot,length,mode);

  TIDY(data);
  return(0);
}/*main*/


/*##########################################*/
/*output the data*/

void outputData(char *data,int numb,int nBins,char *outRoot,int length,char mode)
{
  int i=0,j=0,offset=0;
  int minR=0,maxR=0,r=0;
  char interesting=0;
  char namen[200];
  char thresh=0;
  char max=0;
  FILE *opoo=NULL;

  thresh=-110;

  if((mode==0)||(mode==2)||(mode==3)){
    for(i=0;i<3200;i++){  /*zenith loop*/
      if(mode==0)interesting=1;  /*print out everything*/
      else if(mode==2){          /*only print out multi-returns*/
        interesting=0;
        for(j=0;j<nBins;j++){
          offset=i*nBins+j;
          if(data[offset]>thresh){
            for(;j<nBins;j++){   /*step to the end of the feature*/
              offset=i*nBins+j;
              if(data[offset]<=thresh){
                interesting++;
                break;
              }
            } /*step to the end of the feature*/
          } /*found something*/
        }/*bin loop*/
        if(interesting>4)interesting=1;
        else             interesting=0;
      }else if(mode==3){
        interesting=0;
        minR=(int)(140.0/0.15);
        maxR=(int)(147.0/0.15);
        for(j=minR-5;j<maxR+5;j++){
          offset=i*nBins+j;
          if(data[offset]>thresh){
            max=-125;
            for(;j<maxR+5;j++){   /*step to the end of the feature*/
              offset=i*nBins+j;
              if(data[offset]>max){
                max=data[offset];
                r=j;
              }
              if(data[offset]<=thresh){
                if((r>=minR)&&(r<=maxR))interesting=1;
                break;
              }
            } /*step to the end of the feature*/
          } /*found something*/
        }/*bin loop*/
      } /*interesting check*/

      if(interesting){   /*only print out intereting waveforms*/
        sprintf(namen,"%s.%d.dat",outRoot,i);
        if((opoo=fopen(namen,"w"))==NULL){
          fprintf(stderr,"Error opening output file %s\n",namen);
          exit(1);
        }
        fprintf(opoo,"# 1 index, 2 intensity\n");
        for(j=0;j<nBins;j++){ /*bin loop*/
          offset=i*nBins+j;
          if(offset>=length){
            fprintf(stderr,"Offset error %d %d\n",offset,length);
            exit(1);
          }
          fprintf(opoo,"%g %d\n",(float)j*0.15,data[offset]);
        }/*bin loop*/
        if(opoo){
          fclose(opoo);
          opoo=NULL;
        }
        fprintf(stdout,"Written to %s\n",namen);
      }/*interest test*/
    }/*zenith loop*/
  }else if(mode==1){
    sprintf(namen,"%s.dat",outRoot);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
    for(i=0;i<length;i++)fprintf(opoo,"%d %d\n",i,data[i]);
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",namen);
  }else{
    fprintf(stderr,"No mode %d yet\n",mode);
    exit(1);
  }

  return;
}/*outputData*/


/*##########################################*/
/*read data into array*/

char *readData(char *namen,int numb,int nBins,int *length)
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
  fprintf(stdout,"File is %d long\n",*length);

  nBins=1200;   /**length/numb;*/
  data=challoc(*length,"data",0);

  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to file start*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  if(fread(&(data[0]),sizeof(char),*length,ipoo)!=*length){
    fprintf(stderr,"error reading data\n");
    exit(1);
  }
  fprintf(stdout,"Read\n");
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readData*/


/*the end*/
/*##########################################*/

