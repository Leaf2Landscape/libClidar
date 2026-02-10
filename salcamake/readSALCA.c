#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

/*##########################*/
/*# Reads raw SALCA binary #*/
/*# for P Bunting to get   #*/
/*# into SPD lib.          #*/
/*# Newcastle Nov 2013     #*/
/*##########################*/


/*library functions*/
#define TIDY(arr) if((arr)){free((arr));(arr)=NULL;}  /*free an array*/
float *falloc(int,char *,int);
char *challoc(int,char *,int);
int *ialloc(int,char *,int);
void checkArguments(int,int,int,char *);



/*######################################################################*/

typedef struct{ /*to hold control options*/
  float maxR;
  int nAz;
  int nZen;
  float azStep;
  float zStep;
  float maxZen;
  float azStart;
  float azSquint;
  float zenSquint;
  float *zen;       /*true zenith*/
  float *azOff;     /*azimuth offset*/
  float omega;      /*mirror slope angle*/
  char func;        /*perform function fitting*/
}control;


/*######################################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;          /*loop controls*/
  char inRoot[200];   /*input filename root*/
  void splitShots(control *,char *);
  //void pointOut(char *,int,int,FILE *,float,int,int,int,char,char,calibration *,control *,int,int,float *);
  control *options=NULL;


  /*the defaults are for the hedge test*/
  strcpy(inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  if(!(options=(control *)calloc(2,sizeof(control)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  options->maxR=30.0;
  options->azStep=0.06;   /*in degrees*/
  options->maxZen=190.0;   /*full hemispheric scan*/
  options->azSquint=1.6*M_PI/180.0;
  options->zenSquint=0.0;
  options->omega=M_PI/4.0;  /*45 degrees*/
  options->azStart=0.0;
  options->nZen=3200;
  options->nAz=666;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        options->nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azStep",7)){
        checkArguments(1,i,argc,"-azStep");
        options->azStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        options->maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        options->maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azSquint",9)){
        checkArguments(1,i,argc,"-azSquint");
        options->azSquint=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-zenSquint",10)){
        checkArguments(1,i,argc,"-zenSquint");
        options->zenSquint=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-omega",6)){
        checkArguments(1,i,argc,"-omega");
        options->omega=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-azStart",8)){
        checkArguments(1,i,argc,"-azStart");
        options->azStart=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in degrees\n-maxR range;    maximum recorded range, metres\n-maxZen zen;    maximum zenith angle, degrees. 190 by default\n-azSquint angle;    azimuth squint angle in degrees\n-zenSquint angle;   zenith squint angle, degrees\n-omega angle;       mirror angle in degrees\n-azStart angle;     azimuth start, degrees\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry readSALCA -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/


  /*split shots up*/
  splitShots(options,inRoot);


  if(options){
    TIDY(options->azOff);
    TIDY(options->zen);
    TIDY(options);
  }
  return(0);
}/*main*/


/*#########################################################*/
/*read through data and split into shots*/

void splitShots(control *options,char *inRoot)
{
  int a=0,z=0,band=0;  /*loop control*/
  int numb=0,nBins=0;  /*number of zenith steps and range bins per file (one file per azimuth)*/
  int sOffset=0;       /*bounds to avoid outgoing pulse*/
  int start[2],end[2]; /*array start and end bounds*/
  int length=0;        /*total file length*/
  int waveStart=0;     /*waveform start in data array*/
  int findStart(int,int,char *,char *,int);
  float zen=0,az=0;    /*true zenith and azimuth*/
  char *data=NULL;     /*waveform data. Holds whole file.*/
  char *readData(char *,int,int *,int *,int *,control *);
  char ringTest=0;     /*ringing and saturation indicator*/

  /*for bad SALCA geometry*/
  void translateSquint(control *);
  void setSquint(control *,int);



  sOffset=7;   /*1.05 m*/


  numb=(int)(options->maxZen/0.059375); /*number of zenith steps. Zenith res is hard wired for SALCA*/
  nBins=(int)((150.0+options->maxR)/0.15);       /*150m of 1550 nm plus maxR of 1040 nm*/

  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(options->maxR/0.15);
  start[1]=1000+sOffset;   /*7 to avoid outgoing pulses*/
  end[1]=nBins;

  /*set up squint angles*/
  translateSquint(options);
  setSquint(options,numb);

  /*loop through azimuth files and write out points*/
  for(a=0;a<=options->nAz;a++){ /*azimuth loop*/
    data=readData(inRoot,a,&numb,&nBins,&length,options); /*read binary data into an array*/

    if(data){ /*does the file exist?*/
      for(z=0;z<numb;z++){ /*zenith loop*/
        zen=options->zen[z];
        az=options->azOff[z]+(float)a*options->azStep+options->azStart;

        for(band=0;band<2;band++){   /*loop through wavebands*/
          /*determine shot start position*/
          waveStart=findStart(start[band],end[band],&ringTest,data,z*nBins);
          if((band==0)&&(waveStart<0))waveStart=0; /*as the 1545 pulse is often lost off the end*/


          /*the wave runs from data[z*nBins+waveStart] to data[z*nBins+end[band]] */
          fprintf(stdout,"Wave az %d zen %d band %d runs from %d to %d\n",a,z,band,z*nBins+waveStart,z*nBins+end[band]);


        }/*band loop*/
      }/*zenith loop*/
      TIDY(data);   /*clean up arrays as we go along*/
    }/*does file exist?*/
  }/*azimuth loop*/

  return;
}/*splitShots*/


/*##########################################*/
/*find outgoing pulse and check saturation*/

int findStart(int start,int end,char *satTest,char *data,int offset)
{
  int i=0,b=0,e=0;  /*loop control and bounds*/
  int place=0;      /*array index*/
  int waveStart=0;  /*outgoing pulse bin*/
  char max=0;       /*max intensity*/
  char satThresh=0; /*saturation threshold*/
  char thresh=0;
  int sOffset=0;


  sOffset=7;   /*1.05 m*/
  waveStart=-1; /*start-sOffset;*/

  thresh=-110;  /*noise threshold, chosen from histograms*/
  satThresh=127;
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
  char *data=NULL;
  char namen[200];
  FILE *ipoo=NULL;

  sprintf(namen,"%s_%d.bin",inRoot,i);                  /*input filename*/
  fprintf(stdout,"%d of %d Reading %s\n",i,options->nAz,namen);  /*progress indicator*/

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

  if(((*nBins)*(*numb))>(*length)){
  fprintf(stderr,"File size mismatch\n");
    exit(1);
  }
  data=challoc(*length,"data",0);

  /*now we know how long, read the file*/
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



/*#####################################################################*/
/*precalculate squint angles*/

void setSquint(control *options,int numb)
{
  int i=0;
  float zen=0,az=0;
  float cZen=0,cAz=0;
  void squint(float *,float *,float,float,float,float,float);

  options->zen=falloc(numb,"zenith squint",0);
  options->azOff=falloc(numb,"azimuth squint",0);

  az=0.0;
  for(i=0;i<numb;i++){
    zen=((float)(options->nZen/2)-(float)i)*options->zStep;
    squint(&(cZen),&(cAz),zen,az,options->zenSquint,options->azSquint,options->omega);
    options->zen[i]=cZen;
    options->azOff[i]=cAz;
  }/*zenith loop*/

  return;
}/*setSquint*/


/*#####################################################################*/
/*caluclate squint angle*/

void squint(float *cZen,float *cAz,float zM,float aM,float zE,float aE,float omega)
{
  float inc=0;  /*angle of incidence*/
  void rotateX(float *,float);
  void rotateY(float *,float);
  void rotateZ(float *,float);
  float *vect=NULL;
  /*working variables*/
  float mX=0,mY=0,mZ=0; /*mirror vector*/
  float lX=0,lY=0,lZ=0; /*incoming laser vector*/
  float rX=0,rY=0,rZ=0; /*vector orthogonal to m and l*/
  float thetaZ=0;       /*angle to rotate to mirror surface about z axis*/
  float thetaX=0;       /*angle to rotate about x axis*/
  float slope=0;        /*rotstion vector slope angle, for rounding issues*/
  /*trig*/
  float coszE=0,sinzE=0;
  float cosaE=0,sinaE=0;
  float coszM=0,sinzM=0;
  float cosW=0,sinW=0;

  coszE=cos(zE);
  sinzE=sin(zE);
  cosaE=cos(aE);
  sinaE=sin(aE);
  cosW=cos(omega);
  sinW=sin(omega);
  coszM=cos(zM);
  sinzM=sin(zM);

  mX=cosW;        /*mirror normal vector*/
  mY=sinW*sinzM;
  mZ=sinW*coszM;
  lX=-1.0*coszE;  /*laser Poynting vector*/
  lY=sinaE*sinzE;
  lZ=cosaE*sinzE;
  rX=lY*mZ-lZ*mY; /*cross product of mirror and laser*/
  rY=lZ*mX-lX*mZ; /*ie the vector to rotate about*/
  rZ=lX*mY-lY*mX;

  inc=acos(-1.0*mX*lX+mY*lY+mZ*lZ);   /*angle of incidence. Reverse x to get acute angle*/
  thetaZ=-1.0*atan2(rX,rY);
  thetaX=atan2(rZ,sqrt(rX*rX+rY*rY));

  vect=falloc(3,"vector",0);
  vect[0]=lX;
  vect[1]=lY;
  vect[2]=lZ;

  /*to avoid rounding rotate to z or y axis as appropriate*/
  slope=atan2(sqrt(rX*rX+rY*rY),fabs(rZ));
  if(fabs(slope)<(M_PI/4.0)){  /*rotate about z axis*/
    thetaX=-1.0*atan2(sqrt(rX*rX+rY*rY),rZ);
    thetaZ=-1.0*atan2(rX,rY);
    rotateZ(vect,thetaZ);
    rotateX(vect,thetaX);
    rotateZ(vect,-2.0*inc);
    rotateX(vect,-1.0*thetaX);
    rotateZ(vect,-1.0*thetaZ);
  }else{                        /*rotate about y axis*/
    thetaZ=-1.0*atan2(rX,rY);
    thetaX=atan2(rZ,sqrt(rX*rX+rY*rY));
    rotateZ(vect,thetaZ);
    rotateX(vect,thetaX);
    rotateY(vect,-2.0*inc);
    rotateX(vect,-1.0*thetaX);
    rotateZ(vect,-1.0*thetaZ);
  }

  *cZen=atan2(sqrt(vect[0]*vect[0]+vect[1]*vect[1]),vect[2]);
  if(vect[1]!=0.0)*cAz=atan2(vect[0],vect[1])+aM;
  else            *cAz=aM;

  TIDY(vect);
  return;
}/*squint*/


/*########################################################################*/
/*rotate about x axis*/

void rotateX(float *vect,float theta)
{
  int i=0;
  float temp[3];

  temp[0]=vect[0];
  temp[1]=vect[1]*cos(theta)+vect[2]*sin(theta);
  temp[2]=vect[2]*cos(theta)-vect[1]*sin(theta);

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateX*/


/*########################################################################*/
/*rotate about y axis*/

void rotateY(float *vect,float theta)
{
  int i=0;
  float temp[3];

  temp[0]=vect[0]*cos(theta)-vect[1]*sin(theta);
  temp[1]=vect[1];
  temp[2]=vect[0]*sin(theta)+vect[2]*cos(theta);

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateX*/


/*########################################################################*/
/*rotate about z axis*/

void rotateZ(float *vect,float theta)
{
  int i=0;
  float temp[3];

  temp[0]=vect[0]*cos(theta)+vect[1]*sin(theta);
  temp[1]=vect[1]*cos(theta)-vect[0]*sin(theta);
  temp[2]=vect[2];

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateX*/


/*########################################################################*/
/*translate from nice squint angles to those used in equations*/

void translateSquint(control *options)
{
  float sinAz=0,sinZen=0;

  sinZen=sin(options->zenSquint);
  sinAz=sin(options->azSquint);

  options->azSquint=atan2(sinAz,sinZen);
  options->zenSquint=atan2(sqrt(sinAz*sinAz+sinZen*sinZen),1.0);

  return;
}/*translateSquint*/


/*########################################################################*/
/*########################################################################*/
/*library functions*/

/*#########################################################################*/
/*allocate an int arrays, checking for errors*/

int *ialloc(int length,char *namen,int n)
{
  int *jimlad=NULL;
  if(!(jimlad=(int *)calloc(length,sizeof(int)))){
    fprintf(stderr,"error in %s array allocation %d\n",namen,n);
    fprintf(stderr,"allocating %d\n",length);
    exit(1);
  }
  return(jimlad);
}/*ialloc*/

/*#########################################################################*/
/*allocate a char arrays, checking for errors*/

char *challoc(int length,char *namen,int n)
{
  char *jimlad=NULL;
  if(!(jimlad=(char *)calloc(length,sizeof(char)))){
    fprintf(stderr,"error in %s array allocation %d\n",namen,n);
    fprintf(stderr,"allocating %d\n",length);
    exit(1);
  }
  return(jimlad);
}/*challoc*/


/*########################################################################*/
/*allocate a float array, checking for errors*/

float *falloc(int length,char *namen,int n)
{
  float *jimlad=NULL;
  if(!(jimlad=(float *)calloc(length,sizeof(float)))){
    fprintf(stderr,"error in %s array allocation %d\n",namen,n);
    fprintf(stderr,"allocating %d\n",length);
    exit(1);
  }
  return(jimlad);
}/*falloc*/


/*#########################################################################*/
/*check the number of command line arguments*/

void checkArguments(int numberOfArguments,int thisarg,int argc,char *option)
{
  int i=0;

  for(i=0;i<numberOfArguments;i++){
    if(thisarg+1+i>=argc){
      fprintf(stderr,"error in number of arguments for %s option: %d required\n",option,numberOfArguments);
      exit(1);
    }
  }

  return;
}/*checkArguments*/

/*the end*/
/*##########################################*/

