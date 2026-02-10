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

int main(int argc,char **argv)
{
  int i=0,band=0;   /*loop controls*/
  int numb=0;       /*number of zen zteps*/
  int nBins=0;      /*number of range bins*/
  int length=0;     /*total file length*/
  int nAz=0;        /*number of azimuth steps*/
  int start[2],end[2]; /*array start and end bounds*/
  int **histogram=NULL;
  float azStep=0;   /*azimuth step*/
  float maxR=0;     /*maximum range*/
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int *,int *,int *);
  char inRoot[200],namen[200];   /*input filename root and filename*/
  char outNamen[200];             /*output filename root*/
  void markHisto(char *,int,int,int,int,int *);
  void writeHisto(int **,char *);

  /*the defaults are for the hedge test*/
  strcpy(inRoot,"/Users/growler/data/SALCA/raw/hedge/hedge_plot34_f1_1936");
  nAz=666;
  strcpy(outNamen,"salcaTest.hist");
  azStep=0.001;   /*in radians*/
  maxR=30.0;
  numb=3200;   /*number of azimuth steps*/
  nBins=1200;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",8)){
        checkArguments(1,i,argc,"-output");
        strcpy(outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azStep",7)){
        checkArguments(1,i,argc,"-azStep");
        azStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-output name;  output filename\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in radians\n-maxR range;    maximum recorded range, metres\n\n");
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
  start[0]=7;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(maxR/0.15);
  start[1]=1000+7;
  end[1]=nBins;

  /*allocate histogram*/
  histogram=iIalloc(2,"histogram",0);
  for(band=0;band<2;band++){
    histogram[band]=ialloc(256,"histogram",band+1);
    for(i=0;i<255;i++)histogram[band][i]=0;
  }

  /*loop through azimuth files and write out points*/
  for(i=0;i<nAz;i++){ /*az loop*/
    sprintf(namen,"%s_%d.bin",inRoot,i);                  /*input filename*/
    fprintf(stdout,"%d of %d Reading %s\n",i,nAz,namen);  /*progress indicator*/
    data=readData(namen,&numb,&nBins,&length); /*read binary data into an array*/

    if(data){  /*check for usable data*/
      for(band=0;band<2;band++){   /*loop through wavebands*/
        markHisto(data,numb,nBins,start[band],end[band],histogram[band]);  /*write out point cloud*/
      }/*band loop*/
    }
    TIDY(data);   /*clean up arrays as we go along*/
  }/*az loop*/


  /*write out histogram*/
  writeHisto(histogram,outNamen);

  TTIDY((void **)histogram,2);
  return(0);
}/*main*/


/*##########################################*/
/*write out histogram*/

void writeHisto(int **histogram,char *outNamen)
{
  int i=0;
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }

  for(i=0;i<=255;i++)fprintf(opoo,"%d %d %d\n",i,histogram[0][i],histogram[1][i]);

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);
  return;
}/*writeHisto*/


/*##########################################*/
/*output a point cloud*/

void markHisto(char *data,int numb,int nBins,int start,int end,int *histogram)
{
  int i=0,j=0;  /*loop variables*/
  int place=0;  /*array place*/
  int index=0;  /*range index*/
  int nIn=0;    /*number of features*/
  int waveStart=0; /*outgoing pulse bin*/
  int findStart(int,int,char *,char *,int);
  float zStep=0,zen=0; /*zenith and step in radians*/
  char thresh=0,max=0; /*intensity threshold and max value*/
  char satTest=0;      /*saturation test indicator*/

  thresh=-110;        /*chosen by RG, looks good*/
  zStep=0.001047198;  /*fixed for SALCA*/

  for(i=0;i<numb;i++){ /*zenith loop*/
    //zen=((float)(numb/2+200)-(float)i)*zStep;   /*I have no idea why it's -200*/
    zen=((float)(numb/2)-(float)i)*zStep;
    nIn=1;                      /*number of points per waveform*/

    /*find the waveform start and check for saturation*/
    waveStart=findStart(start,end,&satTest,data,i*nBins);

    for(j=start;j<end;j++){     /*loop along waveform*/
      place=i*nBins+j;            /*array place*/
      if(data[place]>thresh){   /*signal above noise level*/
        max=data[place];          /*reser max*/
        for(;j<end;j++){        /*step to end of feature*/
          place=i*nBins+j;        /*array place*/
          if(data[place]>max){  /*record maximum intensity and position*/
            max=data[place];      /*set max*/
            index=place;
          }                     /*max test*/
          if(data[place]<=thresh){
            histogram[(int)(data[index]+125)]++;
            nIn++;    /*record number of points per waveform*/
            break;
          }/*left feature*/
        }/*point end*/
      }/*point start*/
    }/*bin loop*/
  }/*zenith loop*/
  return;
}/*pointOut*/


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

  satThresh=127;
  thresh=-110;
  max=-125;

  *satTest=0;  /*not saturated by default*/
  b=start-15;   /*1.05m, chosen by RG*/
  if(b<0)b=0;  /*truncate at 0*/
  e=start+7;
  if(e>end)e=end;

  for(i=b;i<e;i++){  /*loop around where we think the pulse might be*/
    place=offset+i;
    if((data[place]>thresh)&&(data[place]>max)){ /*outgoing peak*/
      max=data[place];
      waveStart=i;
    }                    /*outgoing peak test*/
  }/*range loop*/

  //for(i=start;i<end;i++){  /*loop through full waveform to test for saturation*/
  //  place=offset+i;
  //  if(data[place]>=satThresh){   /*saturated*/
  //    *satTest=1;
  //    break;
  //  }  
  //}/*saturation test loop*/

  return(waveStart);
}/*findStart*/


/*##########################################*/
/*mark a point down*/

void markPoint(float zen,float az,float r,FILE *opoo,char max,int nIn,int azInd,int zenInd)
{
  float x=0,y=0,z=0;  /*cartesian coordinates*/
  float offset=0;     /*for the wonky beam correction*/
  float cAz=0,cZen=0; /*for the wonky beam correction*/

  /*the 1.6 deg wonky beam correction*/
  /*we could pre-calculate the offset and trig values*/
  /*to save time. It would need storing in a structure*/
  offset=1.6*M_PI/180.0;
  cAz=az-atan2(tan(offset),sin(zen));
  cZen=acos(cos(offset)*cos(zen));

  x=r*sin(cZen)*cos(cAz);
  y=r*sin(cZen)*sin(cAz);
  z=r*cos(cZen);   /*I don't know why this should be negative, unless my offset is wrong. DO a trip*/

  fprintf(opoo,"%g %g %g %d %g %g %g %d %d %d\n",x,y,z,max,zen*180.0/M_PI,az*180.0/M_PI,r,nIn,azInd,zenInd);

  return;
}/*markPoint*/


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
    return(NULL);
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
/*open two output pointers*/

FILE **openOutput(char *outRoot)
{
  int band=0;
  FILE **opoo=NULL;
  char outNamen[200];

  if(!(opoo=(FILE **)calloc(2,sizeof(FILE *)))){
    fprintf(stderr,"error in file pointer.\n");
    exit(1);
  }

  /*band loop*/
  for(band=0;band<2;band++){
    sprintf(outNamen,"%s.band.%d.coords",outRoot,band);
    if((opoo[band]=fopen(outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",outNamen);
      exit(1);
    }
    fprintf(opoo[band],"# 1 x, 2 y, 3 z, 4 int, 5 zen, 6 az, 7 range, 8 number, 9 azInd, 10 zenInd\n");
  }
  return(opoo);
}/*openOutput*/


/*##########################################*/
/*close output files*/

void closeFiles(FILE **opoo,char *outRoot)
{
  int band=0;

  for(band=0;band<2;band++){
    if(opoo[band]){
      fclose(opoo[band]);
      opoo[band]=NULL;
    } 
    fprintf(stdout,"Written to %s.band.%d.coords\n",outRoot,band);
  }
  TIDY(opoo);

  return;
}/*closeFiles*/


/*the end*/
/*##########################################*/

