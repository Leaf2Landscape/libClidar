#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*##########################*/
/*# Reads raw SALCA binary #*/
/*# and outputs an ascii   #*/
/*# point cloud            #*/
/*##########################*/

int sOffset;   /*wavestart offset*/
char backSig;  /*background signal*/


typedef struct{ /*to hold calibration data*/
  float *LUT;   /*DN to reflectance array*/
  float minDN;  /*min DN in LUT*/
  float maxDN;  /*max DN in LUT*/
  float tran;   /*filter transmission*/
  int numb;     /*number of elements in LUT*/
}calibration;

typedef struct{ /*to hold control options*/
  int coarsen;  /*aggragate beams*/
  char dualOut;     /*combined/seperate waveband output switch*/
  char calibrate;   /*calibration switch*/
  char intTrack; /*track intensity to background or not*/
  char subBin;   /*sub bin range resolution flag*/
  int nAz;
  float azStep;
  float zStep;
  float maxZen;
  int zenOffset[2]; /*zenith offset*/
}control;



/*######################################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0,band=0;   /*loop controls*/
  int numb=0;       /*number of zen zteps*/
  int nBins=0;      /*number of range bins*/
  int length=0;     /*total file length*/
  int nAz=0;        /*number of azimuth steps*/
  int start[2],end[2]; /*array start and end bounds*/
  float maxR=0;     /*maximum range*/
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int,int *,int *,int *,control *);
  char inRoot[200];   /*input filename root*/
  char output[200];             /*output filename root*/
  void findPulses(char *,int,int,FILE *,float,int,int,int,control *,int,int);
  control *options=NULL;
  FILE *opoo=NULL;              /*pointer to output files*/

  /*the defaults are for the hedge test*/
  strcpy(inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  nAz=666;
  strcpy(output,"pulses.dat");
  maxR=30.0;
  nBins=1200;

  sOffset=7;   /*1.05 m*/
  backSig=-116;

  if(!(options=(control *)calloc(2,sizeof(control)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  options->coarsen=1;   /*use at native resolution*/
  options->dualOut=0;
  options->calibrate=0;/*output DN rather than reflectance*/
  options->intTrack=0;  /*use threshold rather than background*/
  options->subBin=0;    /*use peak intensity bin only*/
  options->azStep=0.001;   /*in radians*/
  options->zenOffset[0]=options->zenOffset[1]=0;
  options->maxZen=190.0;   /*full hemispheric scan*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(output,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azStep",7)){
        checkArguments(1,i,argc,"-azStep");
        options->azStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-zenOffset",10)){
        options->zenOffset[1]=1;
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        options->maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-outRoot root;  output filename root\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in radians\n-maxR range;    maximum recorded range, metres\n-subBin;    sub-bin range resolution;\n-intTrack;      track to background noise\n-calibrate;     calibrate to reflectance, need calibration file\n-calFile name;   calibrate to reflectance, need calibration file\n-filt n;     filter optical depth (0, 0.6, 1 or 1.6)\n-nFilt t1 t2;     specify transmissions rather than use defaults\n-dualOut;      output a combined point cloud\n-coarsen n;     coarsen by a factor. CAUTIION, intensities will not scale properly yet\n-zenOffset;      offset 1064nm zenith by 1\n-maxZen zen;    maximum zenith angle, degrees. 190 by default\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/


  //numb=3200;   /*number of zenith steps*/
  numb=(int)(options->maxZen/0.059375);

  /*set scan array dimensions*/
  nBins=(int)((150.0+maxR)/0.15);  /*150m of 1550 nm plus maxR of 1040 nm*/

  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(maxR/0.15);
  start[1]=1000+sOffset;   /*7 to avoid outgoing pulses*/
  end[1]=nBins;

  /*open the two output files, one for each band*/
  if((opoo=fopen(output,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",output);
    exit(1);
  }
  fprintf(opoo,"# 1 zen, 2 az, 3 zenInd, 4 azInd, 5 peak DN, 6 integral, 7 position, 8 number\n");

  /*if coarsened*/
  nAz/=options->coarsen;
  numb/=options->coarsen;
  options->azStep*=(float)options->coarsen;
  options->nAz=nAz;
  

  options->zStep=0.001047198*(float)options->coarsen;  /*fixed for SALCA*/


  /*loop through azimuth files and write out points*/
  band=1;  /*only look 1t 1064nm*/
  for(i=0;i<=nAz;i++){ /*az loop*/
    data=readData(inRoot,i,&numb,&nBins,&length,options); /*read binary data into an array*/

    if(data){ /*does the file exist?*/
      findPulses(data,numb,nBins,opoo,(float)i*options->azStep,start[band],end[band],i,options,options->zenOffset[band],band);
    }/*does file exist?*/
    TIDY(data);   /*clean up arrays as we go along*/
  }/*az loop*/

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Data written to %s\n",output);
  TIDY(options);
  return(0);
}/*main*/


/*##########################################*/
/*find outgoing 1064 pulses*/

void findPulses(char *data,int numb,int nBins,FILE *opoo,float az,int start,int end,int azInd,control *options,int zOff,int band) 
{
  int i=0,bin=0,place=0;
  int n=0;
  int waveStart=0; /*outgoing pulse bin*/
  int findStart(int,int,char *,char *,int);
  int integral=0;
  char ringTest=0;        /*test for ringing*/
  char max=0;
  char thresh=0;
  float zen=0;

  thresh=-110;

  for(i=0;i<numb;i++){ /*zenith loop*/
    zen=((float)(3200/2)-(float)(i+zOff))*options->zStep;
    waveStart=findStart(start,end,&ringTest,data,i*nBins);
      if(waveStart>=0){
      /*track back and for*/
      integral=0;
      max=-127;
      n=0;
      for(bin=waveStart;bin>=0;bin--){  /*step forwards*/
        place=i*nBins+bin;            /*array place*/
        if(data[place]>max)max=data[place];
        if(data[place]<=thresh)break;
        integral+=data[place]-backSig;
        n++;
      }
      for(bin=waveStart+1;bin<nBins;bin++){
        place=i*nBins+bin;            /*array place*/
        if(data[place]>max)max=data[place];
        if(data[place]<=thresh)break;
        integral+=data[place]-backSig;
        n++;
      }

      /*output results*/
      fprintf(opoo,"%f %f %d %d %d %d %d %d\n",zen,az,i,azInd,max,integral,waveStart,n);
    }/*found start test*/
  }/*zen loop*/

  return;
}/*findPulses*/


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

char *readData(char *inRoot,int i,int *numb,int *nBins,int *length,control *options)
{
  int j=0,k=0,m=0,bin=0;
  int place=0,cPlace=0;
  int *contN=NULL;
  int *coarse=NULL;
  char *data=NULL;
  char namen[200];
  FILE *ipoo=NULL;

  coarse=ialloc((*numb)*(*nBins),"data",0);
  contN=ialloc((*numb)*(*nBins),"counter",0);

  for(j=0;j<options->coarsen;j++){  /*coarsen*/
    sprintf(namen,"%s_%d.bin",inRoot,i*options->coarsen+j);                  /*input filename*/
    fprintf(stdout,"%d of %d Reading %s\n",i*options->coarsen+j,options->nAz*options->coarsen,namen);  /*progress indicator*/

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
    if(options->coarsen>1){ /*if coarsening, copy data about*/
      for(k=0;k<*numb;k++){
        for(bin=0;bin<*nBins;bin++){
          cPlace=k*(*nBins)+bin;
          for(m=0;m<options->coarsen;m++){
            place=(k*options->coarsen+m)*(*nBins)+bin;
            coarse[cPlace]+=data[place]+127;
            contN[cPlace]++;
          }
        }
      }
      TIDY(data);
    }/*if coarsening, copy data about*/
  }/*coarsening loop*/

  if(options->coarsen>1){ /*if coarsening, copy data about*/
    data=challoc(*numb*(*nBins),"data",0);
    for(k=0;k<*numb;k++){
      for(bin=0;bin<*nBins;bin++){
        place=k*(*nBins)+bin;
        if(contN[place]>0)data[place]=(int)((float)coarse[place]/(float)contN[place]-127.0);
      }
    }
  }/*if coarsening, copy data about*/


  TIDY(coarse);
  TIDY(contN);

  return(data);
}/*readData*/


/*the end*/
/*##########################################*/

