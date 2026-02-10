#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*#######################*/
/*# Reads lots of SALCA #*/
/*# and tries to find   #*/
/*# the pulse shape     #*/
/*#######################*/


int nZens;  /*number of zeniths in a slice*/
float sRes;  /*SALCA res*/

typedef struct{
  char **fileRoots;
  float *maxRlist;
  float **rBounds;
  int **azBounds;
  int **zenBounds;
  int nFiles;
}inList;

typedef struct{
  float res;   /*sampling resolution*/
  int length;  /*pulse array length*/
  float maxL;  /*maximum likely pulse length*/
  float ***hist; /*combined pulse shape*/
  int maxW;
  char nThresh;
  char background;
  float sWidth;
  float *smoother;
  int sLength;
  char pSmooth;
}pulseShape;

#define TOLERANCE 0.000001


int main(int argc, char **argv)
{
  int i=0;
  inList *control=NULL;
  inList *readInput(char *n);
  pulseShape *pulse=NULL;
  void collatePulses(inList *,pulseShape *);
  void tidyUp(inList *,pulseShape *);
  void outputPulse(pulseShape *,char *);
  void setSmoother(pulseShape *);
  void postSmooth(inList *,pulseShape *);
  char namen[200];  /*input list name*/
  char outNamen[200];


  strcpy(outNamen,"pulseShape2");
  nZens=3200;
  sRes=0.15;
  if(!(pulse=(pulseShape *)calloc(1,sizeof(pulseShape)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  strcpy(namen,"/home/server/users/bakgrp3/nsh103/word/word/cockle/scripts/coords/hardTarg.polar");

  pulse->res=0.01;  /*1 cm*/
  pulse->maxL=10.0; /*10m, why ever not?*/
  pulse->maxW=20;    /*maximum pulse width in bins*/
  pulse->nThresh=-110;
  pulse->background=-115.7;
  pulse->sWidth=1.0;
  pulse->smoother=NULL;
  pulse->pSmooth=0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if (!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(namen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        pulse->res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sWidth",7)){
        checkArguments(1,i,argc,"-sWidth");
        pulse->sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pSmooth",7)){
        pulse->pSmooth=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;     input filename\n-outRoot root;   output filename root\n-res x;          pulse sample resolution\n-sWidth sigma;   smoothinmg function width\n-pSmooth;        smooth histogram after collation\n\n");
        exit(0);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry echiMake -help\n",argv[0],argv[i]);   /* Program_name = argv[0] */
        exit(1);
      }
    }
  }/*read command line*/


  if(fabs(pulse->sWidth)>TOLERANCE)setSmoother(pulse);


  /*read input file list*/
  control=readInput(namen);

  /*add up all the pulses*/
  collatePulses(control,pulse);

  /*post smooth the histogram*/
  if((fabs(pulse->sWidth)>TOLERANCE)&&(pulse->pSmooth))postSmooth(control,pulse);

  /*output*/
  outputPulse(pulse,outNamen);

  tidyUp(control,pulse);

  return(0);
}/*main*/


/*################################################*/
/*smooth histogram*/

void postSmooth(inList *control,pulseShape *pulse)
{
  int band=0,w=0,i=0;
  float *temp=NULL;
  float *smoothFloat(pulseShape *,float *,int);

  for(band=0;band<2;band++){
    //for(w=0;w<pulse->maxW;w++){
    for(w=6;w<=6;w++){
      temp=smoothFloat(pulse,pulse->hist[band][w],pulse->length);
      for(i=0;i<pulse->length;i++){
        pulse->hist[band][w][i]=temp[i];
//fprintf(stdout,"%d %d %d %f %f\n",band,w,i,temp[i],pulse->hist[band][w][i]);
      }
      TIDY(temp);
    }
  }

  return;
}/*postSmooth*/


/*################################################*/
/*output data*/

void outputPulse(pulseShape *pulse,char *outNamen)
{
  int i=0,j=0;
  char namen[200];
  FILE *opoo=NULL;

  for(j=0;j<pulse->maxW;j++){
    sprintf(namen,"%s.width.%d.dat",outNamen,j);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }

    fprintf(opoo,"#\n");
    for(i=0;i<pulse->length;i++){
      fprintf(opoo,"%f %f %f\n",(float)i*pulse->res,pulse->hist[0][j][i],pulse->hist[1][j][i]);
    }

    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",namen);
  }/*width loop*/

  return;
}/*outputPulse*/


/*################################################*/
/*collate pulse shapes*/

void collatePulses(inList *control,pulseShape *pulse)
{
  int i=0,j=0,zen=0;
  int length=0,nBins=0;
  int band=0,offset=0;
  int start=0,end=0;
  int findStart(int,char *);
  char *fullWave=NULL;
  char *readWave(char *,float,int *,int *);
  char namen[200];
  float *smoothed=NULL;
  float *smooth(pulseShape *,char *,int);
  void totPulses(pulseShape *,float *,int,int,int);

  pulse->length=(int)(pulse->maxL/pulse->res)+1;
  if(!(pulse->hist=(float ***)calloc(2,sizeof(float **)))){
    fprintf(stderr,"error in pulse histogram.\n");
    exit(1);
  }
  for(j=0;j<2;j++){
    pulse->hist[j]=fFalloc(pulse->maxW,"pulse histogram",0);
    for(i=0;i<pulse->maxW;i++)pulse->hist[j][i]=falloc(pulse->length,"pulse histogram",i+1);
  }

  for(i=0;i<control->nFiles;i++){
    for(j=control->azBounds[0][i];j<=control->azBounds[1][i];j++){ /*az loop*/
      sprintf(namen,"%s_%d.bin",control->fileRoots[i],j+control->azBounds[0][i]);
      fullWave=readWave(namen,control->maxRlist[i],&length,&nBins);  /*read all zeniths*/
      if(fullWave){
        /*smooth the wave*/
        for(band=0;band<2;band++){
          for(zen=control->zenBounds[0][i];zen<=control->zenBounds[1][i];zen++){ /*zen loop*/
            start=zen*nBins+(int)(control->rBounds[0][i]/sRes);
            end=zen*nBins+(int)(control->rBounds[1][i]/sRes);
            if(band==1){
               offset=findStart(zen*nBins+1000,fullWave);
//fprintf(stdout,"az %d zen %d offset %d file %s\n",j,zen,offset,namen);
               //if(offset<0)break;
               start+=offset;
               end+=offset;
             }

            if((start<0)||(start>=length)||(end<0)||(end>=length)){
              fprintf(stderr,"Array error %d %d for zen %d band %d offset %d\n",start,end,zen,band,offset);
              exit(1);
            }
            smoothed=smooth(pulse,&(fullWave[start]),end-start+1);
            totPulses(pulse,smoothed,0,end-start+1,band);
       
            TIDY(smoothed);
          }/*zen loop*/
        }/*band loop*/
      }
      TIDY(fullWave);
    }/*az loop*/
  }/*file loop*/

  return;
}/*collatePulses*/


/*################################################*/
/*smooth the waveform*/

float *smooth(pulseShape *pulse,char *fullWave,int length)
{
  int i=0,j=0,place=0;
  float *smoothed=NULL;

  smoothed=falloc(length,"smoothed",0);
  if((pulse->sWidth<TOLERANCE)||(pulse->pSmooth)){
    for(i=0;i<length;i++)smoothed[i]=(float)fullWave[i];
  }else{
    for(i=0;i<length;i++){
      smoothed[i]=0.0;
      for(j=0;j<pulse->sLength;j++){
        place=i+(j-pulse->sLength/2);
        if((place>=0)&&(place<length))smoothed[i]+=(float)fullWave[place]*pulse->smoother[j];
      }
    }
  }

  return(smoothed);
}/*smooth*/


/*################################################*/
/*smooth the waveform*/

float *smoothFloat(pulseShape *pulse,float *fullWave,int length)
{
  int i=0,j=0,place=0;
  float *smoothed=NULL;
  float contN=0;

  smoothed=falloc(length,"smoothed",0);
  if(pulse->sWidth<TOLERANCE){
    for(i=0;i<length;i++)smoothed[i]=fullWave[i];
  }else{
    for(i=0;i<length;i++){
      smoothed[i]=0.0;
      contN=0.0;
      for(j=0;j<pulse->sLength;j++){
        place=i+(j-pulse->sLength/2);
        if((place>=0)&&(place<pulse->sLength)){
          smoothed[i]+=fullWave[place]*pulse->smoother[j];
          contN+=pulse->smoother[j];
//fprintf(stdout,"%d %d %d\n",i,j,place);
        }
      }
      if(contN>0.0)smoothed[i]/=contN;
    }
  }

  return(smoothed);
}/*smoothFloat*/


/*################################################*/
/*add up pulses*/

void totPulses(pulseShape *pulse,float *wave,int start,int end,int band)
{
  int i=0,nFeat=0;
  int numb=0,bin=0;
  float contN=0;
  char within=0;
  float mean=0;  /*mean range*/
  float pos=0;   /*bin position*/

  mean=0.0;
  /*there should be only one feature per waveform*/
  within=0;
  numb=0;
  contN=0.0;
  for(i=start;i<end;i++){
    if(wave[i]>(float)pulse->nThresh){
      within=1;
      mean+=sRes*(float)(i-start)*(wave[i]-(float)pulse->background);
      contN+=(wave[i]-(float)pulse->background);
      numb++;
    }
    else{
      if(within)nFeat++;
      within=0;
    }
  }/*number of features check*/

  if((pulse->sWidth>TOLERANCE)&&(!pulse->pSmooth))numb=0;

  if(nFeat>1){
    //fprintf(stdout,"Too many features %d\n",nFeat);
  }else{
    if(numb>=pulse->maxW){
      fprintf(stderr,"Pulse longer than expected %d\n",numb);
      exit(1);
    }
    mean/=contN;  /*this is the pulse centre*/

    /*now add the pulse to the total histogram*/
    /*we should normalise for relfectance, but I'm not sure how to do that without function fitting*/
    for(i=start;i<end;i++){
//if(band==1)fprintf(stdout,"%f\n",wave[i]);
      if(wave[i]>(float)pulse->nThresh){
        pos=(float)(i-start)*sRes-mean;
        bin=(int)(pos/pulse->res+(float)pulse->length/2.0);
        if((bin<0)||(bin>=pulse->length)){
          fprintf(stderr,"Bounding mistake %d of %d\n",bin,pulse->length);
          exit(1);
        }
        pulse->hist[band][numb][bin]+=wave[i]-pulse->background;
      }
    }
  }/*number of features check*/

  return;
}/*totPulses*/


/*################################################*/
/*read a single waveform*/

char *readWave(char *namen,float maxR,int *length,int *nBins)
{
  char *waveform=NULL;
  FILE *ipoo=NULL;

  *nBins=(int)((150.0+maxR)/sRes+0.5);
  *length=*nBins*nZens;
  waveform=challoc(*length,"waveform",0);

  if((ipoo=fopen(namen,"rb"))==NULL){
    //fprintf(stderr,"Error opening input file %s\n",namen);
    return(NULL);
  }
  if(fread(&(waveform[0]),sizeof(char),*length,ipoo)!=*length){
    fprintf(stderr,"error reading data %d from %s\n",*length,namen);
    exit(1);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(waveform);
}/*readWave*/


/*################################################*/
/*read input list*/

inList *readInput(char *namen)
{
  int i=0;
  inList *control=NULL;
  int countLines(FILE *);
  char line[1600];
  char temp1[100],temp2[100],temp3[100],temp4[100],temp5[100],temp6[100],temp7[100],temp8[100];
  FILE *ipoo=NULL;

  if(!(control=(inList *)calloc(1,sizeof(inList)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  control->nFiles=0;
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*read file size*/
  control->nFiles=countLines(ipoo);

  control->fileRoots=chChalloc(control->nFiles,"file root",0);
  control->maxRlist=falloc(control->nFiles,"file root",0);
  control->rBounds=fFalloc(2,"range bounds",0);
  control->azBounds=iIalloc(2,"azimuth bounds",0);
  control->zenBounds=iIalloc(2,"zenith bounds",0);
  for(i=0;i<2;i++){
    control->rBounds[i]=falloc(control->nFiles,"range bounds",i+1);
    control->azBounds[i]=ialloc(control->nFiles,"azimuth bounds",i+1);
    control->zenBounds[i]=ialloc(control->nFiles,"zenith bounds",i+1);
  }


  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(sscanf(line,"%s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8)==8){
      if(strncasecmp(temp1,"#",1)){
        control->fileRoots[i]=challoc(strlen(temp1)+1,"file roots",i+1);
        strcpy(control->fileRoots[i],temp1);
        control->maxRlist[i]=atof(&(temp2[0]));
        control->rBounds[0][i]=atof(&(temp3[0]));
        control->rBounds[1][i]=atof(&(temp4[0]));
        control->azBounds[0][i]=atoi(&(temp5[0]));
        control->azBounds[1][i]=atoi(&(temp6[0]));
        control->zenBounds[0][i]=atoi(&(temp7[0]));
        control->zenBounds[1][i]=atoi(&(temp8[0]));
        i++;
      }
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  return(control);
}/*readInput*/


/*################################################*/
/*count number of lines in file*/

int countLines(FILE *ipoo)
{
  int nLines=0;
  char line[200],temp1[200];

  while(fgets(line,200,ipoo)!=NULL){
    if(sscanf(line,"%s",temp1)==1){
      if(strncasecmp(temp1,"#",1))nLines++;
    }
  }
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  return(nLines);
}/*countLines*/


/*################################################*/
/*tidy control array*/

void tidyUp(inList *control,pulseShape *pulse)
{
  int i=0;

  if(control){
    TTIDY((void **)control->fileRoots,control->nFiles);
    TIDY(control->maxRlist);
    TTIDY((void **)control->zenBounds,2);
    TTIDY((void **)control->azBounds,2);
    TTIDY((void **)control->rBounds,2);
    TIDY(control);
  }

  if(pulse){
    for(i=0;i<2;i++)TTIDY((void **)pulse->hist[i],pulse->maxW);
    TIDY(pulse->smoother);
    TIDY(pulse->hist);
    TIDY(pulse);
  }
  return;
}/*tidyUp*/


/*################################################*/
/*make smoothing function*/

void setSmoother(pulseShape *pulse)
{
  int i=0;
  float x=0,y=0;
  float mu=0;
  float minY=0;

  minY=0.0001;

  /*determine length*/
  for(i=0;;i++){
    x=(float)i*pulse->res;
    y=(float)gaussian((double)x,(double)pulse->sWidth,0.0);
    if(y<minY)break;
  }
  pulse->sLength=2*i;
  pulse->smoother=falloc(pulse->sLength,"smoother",0);

  mu=pulse->res*(float)i;
  for(i=0;i<pulse->sLength;i++){
    x=(float)i*pulse->res;
    pulse->smoother[i]=(float)gaussian((double)x,(double)pulse->sWidth,(double)mu);
  }

  return;
}/*setSmoother*/


/*##########################################*/
/*find outgoing pulse and check saturation*/

int findStart(int start,char *data)
{
  int i=0,b=0,e=0;  /*loop control and bounds*/
  int waveStart=0;  /*outgoing pulse bin*/
  char max=0;       /*max intensity*/
  char thresh=0;    /*noise threshold*/

  waveStart=-1;

  thresh=-110;
  max=-125;

  b=start-50;   /*4.5m, from histograms*/
  if(b<0)b=0;  /*truncate at 0*/
  e=start+7;

  max=-125;
  for(i=b;i<e;i++){  /*loop around where we think the pulse might be*/
//fprintf(stdout,"%d %d\n",i,data[i]);
    if((data[i]>thresh)&&(data[i]>max)){ /*outgoing peak*/
      max=data[i];
      waveStart=i;
    }                    /*outgoing peak test*/
    if((max>-125)&&(data[i]<=thresh))break;
  }/*range loop*/

  return(waveStart-start);
}/*findStart*/

/*the end*/
/*################################################*/

