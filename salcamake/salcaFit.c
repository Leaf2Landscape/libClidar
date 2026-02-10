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

#define TOL 0.00001
#define MAX_ITER 30
int sOffset;   /*wavestart offset*/
char backSig;  /*background signal*/
char thresh;   /*noise threshold*/

int failNumb;


typedef struct{ /*to hold calibration data*/
  float *LUT;   /*DN to reflectance array*/
  float minDN;  /*min DN in LUT*/
  float maxDN;  /*max DN in LUT*/
  float tran;   /*filter transmission*/
  int numb;     /*number of elements in LUT*/
}calibration;

/*######################################################################*/

typedef struct{ /*to hold control options*/
  int coarsen;  /*aggragate beams*/
  char dualOut;     /*combined/seperate waveband output switch*/
  char calibrate;   /*calibration switch*/
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
  char oldSquint;
  char joy;
}control;

/*function to optimise to*/
void gaussErr(float,float *,float *,float *,int);


/*######################################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0,band=0;   /*loop controls*/
  int numb=0;       /*number of zen zteps*/
  int nBins=0;      /*number of range bins*/
  int length=0;     /*total file length*/
  int start[2],end[2]; /*array start and end bounds*/
  float maxR=0;     /*maximum range*/
  float filt=0,filt1=0,filt2=0;     /*filter strength (0, 0.6, 1.0 or 1.6)*/
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int,int *,int *,int *,control *);
  char inRoot[200];   /*input filename root*/
  char outRoot[200];             /*output filename root*/
  char calFile[100]; /*DN to refl calibration file*/
  void pointOut(char *,int,int,FILE *,float,int,int,int,calibration *,control *,int);
  void bothOut(char *,int,FILE *,int,float,float,calibration *,int,control *);
  void closeFiles(FILE **,char *);
  void tidyCal(calibration *);
  calibration *readCalibration(char *,float,float,float);
  calibration *cal=NULL;
  control *options=NULL;
  FILE **openOutput(char *);
  FILE **openDualOutput(char *);
  FILE **opoo=NULL;              /*pointer to output files*/
  void translateSquint(control *);

  /*squint*/
  void setSquint(control *,int);

  /*the defaults are for the hedge test*/
  strcpy(inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  strcpy(outRoot,"salcaTest");
  maxR=30.0;
  nBins=1200;
  strcpy(calFile,"/mnt/geodesy38/nsh103/SALCA/calibration/calRGlambLUT.dat");
  filt=0.0;
  filt1=filt2=-1.0;

  sOffset=7;   /*1.05 m*/
  backSig=-116;
  thresh=-110;

  if(!(options=(control *)calloc(2,sizeof(control)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  options->coarsen=1;   /*use at native resolution*/
  options->dualOut=0;
  options->calibrate=0;/*output DN rather than reflectance*/
  options->azStep=0.06;   /*in degrees*/
  options->maxZen=190.0;   /*full hemispheric scan*/
  options->azSquint=1.6*M_PI/180.0;
  options->zenSquint=0.0;
  options->omega=M_PI/4.0;  /*45 degrees*/
  options->azStart=0.0;
  options->func=0;        /*don't fit a function*/
  options->nZen=3200;
  options->nAz=666;
  options->oldSquint=0;  /*use new suint angles*/
  options->joy=1;        /*print out status*/

  failNumb=0;

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
        options->nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azStep",7)){
        checkArguments(1,i,argc,"-azStep");
        options->azStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-calibrate",10)){
        options->calibrate=1;
      }else if(!strncasecmp(argv[i],"-calFile",8)){
        checkArguments(1,i,argc,"-calFile");
        strcpy(calFile,argv[++i]);
        options->calibrate=1;
      }else if(!strncasecmp(argv[i],"-filt",5)){
        checkArguments(1,i,argc,"-filt");
        filt=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nFilt",6)){
        checkArguments(2,i,argc,"-nFilt");
        filt1=atof(argv[++i]);
        filt2=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-dualOut",8)){
        options->dualOut=1;
        fprintf(stderr,"-dualOut not working yet. Geometry needs sorting\n");
        exit(1);
      }else if(!strncasecmp(argv[i],"-coarsen",8)){
        checkArguments(1,i,argc,"-coarsen");
        options->coarsen=atoi(argv[++i]);
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
      }else if(!strncasecmp(argv[i],"-background",11)){
        checkArguments(1,i,argc,"-background");
        backSig=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-oldSquint",10)){
        options->oldSquint=1;
      }else if(!strncasecmp(argv[i],"-joyless",8)){
        options->joy=0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-outRoot root;  output filename root\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in degrees\n-maxR range;    maximum recorded range, metres\n-calibrate;     calibrate to reflectance, need calibration file\n-calFile name;  calibrate to reflectance, need calibration file\n-filt n;        filter optical depth (0, 0.6, 1 or 1.6)\n-nFilt t1 t2;   specify transmissions rather than use defaults\n-dualOut;       output a combined point cloud\n-coarsen n;     coarsen by a factor. CAUTIION, intensities will not scale properly yet\n-maxZen zen;    maximum zenith angle, degrees. 190 by default\n-azSquint angle;    azimuth squint angle in degrees\n-zenSquint angle;   zenith squint angle, degrees\n-omega angle;       mirror angle in degrees\n-azStart angle;     azimuth start, degrees\n-background DN;     background DN value\n-oldSquint;         use old squint function\n-joyless;           don't print out status\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/



  //numb=3200;   /*number of zenith steps*/
  numb=(int)(options->maxZen/0.059375);


  if(options->calibrate)cal=readCalibration(calFile,filt,filt1,filt2); /*set up calibration bits*/

  /*set scan array dimensions*/
  nBins=(int)((150.0+maxR)/0.15);  /*150m of 1550 nm plus maxR of 1040 nm*/

  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(maxR/0.15);
  start[1]=1000+sOffset;   /*7 to avoid outgoing pulses*/
  end[1]=nBins;

  /*open the two output files, one for each band*/
  if(!options->dualOut)opoo=openOutput(outRoot);
  else        opoo=openDualOutput(outRoot);


  /*if coarsened*/
  if(options->coarsen<1)options->coarsen=1;  /*so zero means uncoarsened*/
  options->nAz/=options->coarsen;
  numb/=options->coarsen;
  options->azStep*=(float)options->coarsen*M_PI/180.0; /*convert to radians as well*/
  options->nZen/=options->coarsen;
  options->zStep=0.001047198*(float)options->coarsen;  /*fixed for SALCA*/

  /*set up squint angles*/
  if(options->oldSquint==0){
    translateSquint(options);
    setSquint(options,numb);
  }


  /*loop through azimuth files and write out points*/
  for(i=0;i<=options->nAz;i++){ /*az loop*/
    data=readData(inRoot,i,&numb,&nBins,&length,options); /*read binary data into an array*/

    if(data){ /*does the file exist?*/
      if(!options->dualOut){   /*output seperate files*/
        for(band=0;band<2;band++){   /*loop through wavebands*/
          pointOut(data,numb,nBins,opoo[band],(float)i*options->azStep+options->azStart,start[band],end[band],i,&(cal[band]),options,band);  /*write out point cloud*/
        }/*band loop*/
      }else{  /*output a single file with both bands*/
        bothOut(data,numb,opoo[0],nBins,(float)i*options->azStep+options->azStart,maxR,cal,i,options);
      }
      TIDY(data);   /*clean up arrays as we go along*/
    }/*does file exist?*/
  }/*az loop*/

  if(!options->dualOut)closeFiles(opoo,outRoot);   /*close output files*/
  else{
    fprintf(stdout,"Written to %s\n",outRoot);
    if(opoo[0]){
      fclose(opoo[0]);
      opoo[0]=NULL;
    }
    TIDY(opoo);
  }

  if(options){
    if(options->calibrate)tidyCal(cal);
    TIDY(options->azOff);
    TIDY(options->zen);
    TIDY(options);
  }
  return(0);
}/*main*/


/*##########################################*/
/*output a combined point cloud*/

void bothOut(char *data,int numb,FILE *opoo,int nBins,float az,float maxR,calibration *cal,int azInd,control *options)
{
  int i=0,j=0,k=0,band=0;
  int bin=0,place=0;
  int maxN=0,nIn[2];
  int **sBin=NULL,**eBin=NULL;
  int **pBin=NULL;
  int start[2],end[2];
  int waveStart=0; /*outgoing pulse bin*/
  float **integral=NULL;
  float zen=0;
  char **satFlag=NULL;
  char max=0;    /*noise threshold and max tracker*/
  char ringTest=0;        /*test for ringing*/
  float findStart(int,int,char *,char *,int);
  void markDual(int,int,float,float,char,char,float,float,FILE *,calibration *,int,int,int,int,control *);

  fprintf(stderr,"This will not work just now\n");
  exit(1);

  maxN=20;   /*maximum likely number, can be increased later*/


  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(maxR/0.15);
  start[1]=1000+sOffset;   /*7 to avoid outgoing pulses*/
  end[1]=nBins;

  sBin=iIalloc(2,"start bin",0);
  eBin=iIalloc(2,"end bin",0);
  pBin=iIalloc(2,"peak bin",0);
  integral=fFalloc(2,"peak bin",0);
  satFlag=chChalloc(2,"saturation flag",0);
  for(band=0;band<2;band++){
    sBin[band]=ialloc(maxN,"start bin",band+1);
    eBin[band]=ialloc(maxN,"end bin",band+1);
    pBin[band]=ialloc(maxN,"peak bin",band+1);
    integral[band]=falloc(maxN,"integral",band+1);
    satFlag[band]=challoc(maxN,"saturation flag",band+1);
  }/*array allocation*/

  for(i=0;i<numb;i++){ /*zenith loop*/
    nIn[0]=nIn[1]=0;   /*reset counter*/
    zen=((float)(options->nZen/2)-(float)i)*options->zStep; //+options->zenSquint; //+(190.0-options->maxZen)*M_PI/180.0;
    for(band=0;band<2;band++){
      /*find the waveform start and check for saturation*/
      waveStart=(int)findStart(start[band],end[band],&ringTest,data,i*nBins);
      if((band==0)&&(waveStart<0))waveStart=0;
      if(waveStart>=0){  /*have we found it?*/
        for(bin=start[band];bin<end[band];bin++){  /*bin loop*/
          place=i*nBins+bin;            /*array place*/
          if(place>=0){
            if(data[place]>thresh){   /*signal above noise level*/
              max=-127;
              sBin[band][nIn[band]]=bin-waveStart;
              satFlag[band][nIn[band]]=0;
              integral[band][nIn[band]]=0.0;

              for(;bin<end[band];bin++){        /*step to end of feature*/
                place=i*nBins+bin;        /*array place*/
                if(place>=0){
                  integral[band][nIn[band]]+=(float)(data[place]-backSig);
                  if(data[place]>max){  /*record maximum intensity and position*/
                    max=data[place];      /*set max*/
                    pBin[band][nIn[band]]=bin-waveStart;
                  }                     /*max test*/
                  if(data[place]==127)satFlag[band][nIn[band]]=1;
                  if(data[place]<=thresh){
                    eBin[band][nIn[band]]=bin-waveStart;
                    break;
                  }
                }/*within data check*/
              }/*signal drops beneath threshold*/

              nIn[band]++;

              if(nIn[band]>=maxN){
                maxN++;
                for(k=0;k<2;k++){
                  if(!(sBin[k]=(int *)realloc(sBin[k],maxN*sizeof(int)))){
                    fprintf(stderr,"Balls\n");
                    exit(1);
                  }
                  if(!(eBin[k]=(int *)realloc(eBin[k],maxN*sizeof(int)))){
                    fprintf(stderr,"Balls\n");
                    exit(1);
                  }
                  if(!(pBin[k]=(int *)realloc(pBin[k],maxN*sizeof(int)))){
                    fprintf(stderr,"Balls\n");
                    exit(1);
                  }
                  if(!(integral[k]=(float *)realloc(integral[k],maxN*sizeof(float)))){
                    fprintf(stderr,"Balls\n");
                    exit(1);
                  }
                  if(!(satFlag[k]=(char *)realloc(satFlag[k],maxN*sizeof(char)))){
                    fprintf(stderr,"Balls\n");
                    exit(1);
                  }
                }/*resize arrays if necessary*/

              }
            }/*signal threshold*/
          }/*within data check*/
        }/*bin loop*/


      }/*found start test*/
    }/*band loop*/

    /*having determined returns, compare the two wavelengths*/
    for(j=0;j<nIn[0];j++){
      for(k=0;k<nIn[1];k++){
        if(((sBin[1][k]>=sBin[0][j])&&(sBin[1][k]<=eBin[0][j]))||((eBin[1][k]>=sBin[0][j])&&(eBin[1][k]<=eBin[0][j]))){
          markDual(pBin[0][j],pBin[1][k]-start[0],integral[0][j],integral[1][k],satFlag[0][j],satFlag[1][k],zen,az,opoo,cal,eBin[0][j]-sBin[0][j],eBin[1][k]-sBin[1][k],i,azInd,options);
          break;
        }/*feature overlap test*/
      }
    }/*having determined returns, compare the two wavelengths*/
  }/*zenith loop*/


  /*tidy up arrays*/
  TTIDY((void **)sBin,2);
  TTIDY((void **)eBin,2);
  TTIDY((void **)pBin,2);
  TTIDY((void **)integral,2);
  TTIDY((void **)satFlag,2);
  return;
}/*bothOut*/


/*############################################################*/
/*mark a dual wavelength point*/

void markDual(int pBin0,int pBin1,float int0,float int1,char sat0,char sat1,float zen,float az,FILE *opoo,calibration *cal,int wid0,int wid1,int zenInd,int azInd,control *options)
{
  float r=0,r0=0,r1=0; /*range*/
  float x=0,y=0,z=0;
  float offset=0;     /*for the wonky beam correction*/
  float cAz=0,cZen=0; /*for the wonky beam correction*/
  float rho0=0,rho1=0;

  if(options->oldSquint==0){
    /*new wonky beam*/
    /*squint(&(cZen),&(cAz),zen,az,options->zenSquint,options->azSquint,options->omega);*/
    cZen=options->zen[zenInd];
    cAz=options->azOff[zenInd]+az;
  }else{
    /*the 1.6 deg wonky beam correction*/
    /*we could pre-calculate the offset and trig values*/
    /*to save time. It would need storing in a structure*/
    offset=options->azSquint;
    cAz=az+atan2(tan(offset),sin(zen));
    cZen=acos(cos(offset)*cos(zen));
  }


  r0=(float)pBin0*0.15;
  r1=(float)pBin1*0.15;

  r=(r0+r1)/2.0;

  x=r*sin(cZen)*cos(cAz);
  y=r*sin(cZen)*sin(cAz);
  z=r*cos(cZen);   /*I don't know why this should be negative, unless my offset is wrong. DO a trip*/

  if(cal){  /*calibrate*/
    if((int0>=0)&&(int0<cal[0].maxDN))rho0=cal[0].LUT[(int)int0]*r0*r0;
    else if(int0<0)                   rho0=0.0;
    else if(int0>=cal[0].maxDN)       rho0=1.0;

    if((int1>=0)&&(int1<cal[1].maxDN))rho1=cal[1].LUT[(int)int1]*r0*r0;
    else if(int1<0)                   rho1=0.0;
    else if(int1>=cal[1].maxDN)       rho1=1.0;
  }else rho0=rho1=-1.0;

  fprintf(opoo,"%f %f %f %d %d %d %d %f %f %f %f %f %d %d\n",x,y,z,(int)int0,(int)int1,sat0,sat1,rho0,rho1,r1-r0,(float)wid0*0.15,(float)wid1*0.15,zenInd,azInd);

  return;
}/*markDual*/


/*##########################################*/
/*output a point cloud*/

void pointOut(char *data,int numb,int nBins,FILE *opoo,float az,int start,int end,int azInd,calibration *cal,control *options,int band)
{
  int i=0,j=0;  /*loop variables*/
  int jS=0;
  int place=0;  /*array place*/
  int nIn=0;    /*number of features*/
  float waveStart=0; /*outgoing pulse bin*/
  float findStart(int,int,char *,char *,int);
  float width=0;       /*feature sigma*/
  float zen=0; /*zenith and step in radians*/
  float range=0;         /*range index*/
  char satTest=0;         /*saturation test indicator*/
  char ringTest=0;        /*test for ringing*/
  char brEak=0;           /*break indicator*/
  char max=0;             /*maximum intensity*/
  char satThresh=0;       /*saturation threshold*/
  float integral=0.0;     /*return energy integral*/
  float rho=0;            /*target reflectance*/
  int indicator=0;        /*function fitting success indicator*/
  int fitFunction(char *,int,float *,float *,float *,float);
  void markPoint(float,float,float,FILE *,char,int,int,int,int,float,char,float,float,control *);
  void writeOne(char *,int,int,float,int);


  satThresh=120;      /*ringing thresholdf*/
  rho=-1.0;           /*nonesense value*/

  for(i=0;i<numb;i++){ /*zenith loop*/
    zen=((float)(options->nZen/2)-(float)i)*options->zStep; //+options->zenSquint; //+(190.0-options->maxZen)*M_PI/180.0;
    nIn=1;                      /*number of points per waveform*/

    /*find the waveform start and check for saturation*/
    waveStart=findStart(start,end,&ringTest,data,i*nBins);
    if((band==0)&&(waveStart<0.0))waveStart=0.0; /*as the 1545 pulse is often lost off the end*/

    if(waveStart>=-100.0){    /*found the outgoing pulse, otherwise we can't trust range*/
      for(j=start;j<end;j++){     /*loop along waveform*/
        place=i*nBins+j;            /*array place*/
        if(data[place]>thresh){   /*signal above noise level*/
          satTest=0;
          max=thresh;          /*reset max*/

          jS=j-1;                   /*mark data for function fitting*/
          for(;j<end;j++){        /*step to end of feature*/
            place=i*nBins+j;        /*array place*/
            if(data[place]>max)max=data[place];

            if(data[place]==127)satTest=1;   /*mark saturation*/

            if(data[place]<=thresh){ /*We have left feature*/
              indicator=fitFunction(&(data[i*nBins+jS]),j-jS+1,&range,&integral,&width,(float)(jS-start)*0.15-waveStart);

              if(options->calibrate){
                if(((int)integral>=0)&&((int)integral<cal->numb))rho=cal->LUT[(int)integral]*range*range;
                else{
                  if((int)integral<=0)rho=0.0;
                  else if((int)integral>=cal->numb)rho=1.0;
                }
              }

              markPoint(zen,az,range,opoo,max,nIn,azInd,i,waveStart,integral,satTest,rho,width,options);

              if(ringTest&&(max>=satThresh))brEak=1;   /*if saturated only take the first return*/

              nIn++;    /*record number of points per waveform*/
              break;
            }/*left feature*/
          }/*point end*/
        }/*point start*/
        if(brEak){   /*for saturated only take first return*/
          brEak=0;   /*this avoids ringing*/
          break;
        }
      }/*bin loop*/
    }/*found wavestart check*/
  }/*zenith loop*/
  return;
}/*pointOut*/


/*##########################################*/
/*fit a a Gaussian by Levenberg-Marquardt*/

int fitFunction(char *data,int numb,float *range,float *integral,float *width,float rOffset)
{
  int i=0;
  int nGauss=0,nParams=0;
  int *pSwitch=NULL;   /*switch to optimise each parameter*/
  int indicator=0;
  float *x=NULL;
  float *temp=NULL;
  float *filterData(char *,int);
  float chisq=0;
  float alambda=0;
  float **covar=NULL;
  float **alpha=NULL;
  float *sig=NULL;
  float *params=NULL;
  float minErr=0,maxErr=0;
  float *initialGuess(float *,float *,int,int);
  int mrqminMine(float *,float *,float *,int,float *,int *,int,float **,float **,float *,void (*)(float,float *,float *,float *,int),float *);
  void writeData(float *,float *,int,float,char *);
  char found=1;

  minErr=0.001;
  maxErr=10.0;

  temp=filterData(data,numb);
  x=falloc(numb+1,"range",0);
  for(i=1;i<=numb;i++)x[i]=(float)(i-1)*0.15+rOffset;

  nGauss=1;
  nParams=3*nGauss;
  params=initialGuess(x,temp,numb,nGauss);

  //if(numb==3){  /*one point above noise, fix max and width*/
  //  params[3]=0.11;
  //  found=1;
  //}else{        /*there are enough points to fit to*/
    pSwitch=ialloc(nParams+1,"parameters",0);

    /*set up L-M arrays*/
    sig=falloc(numb+1,"error",0);
    for(i=1;i<=numb;i++)sig[i]=1.0;
    alpha=fFalloc(nParams+1,"alpha",0);
    covar=fFalloc(nParams+1,"covar",0);
    for(i=nParams;i>=1;i--){
      alpha[i]=falloc(nParams+1,"alpha",i);
      covar[i]=falloc(nParams+1,"covar",i);
      pSwitch[i]=1;
    }

    /*initialise*/
    alambda=-1.0;
    indicator=mrqminMine(x,temp,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);

    if(indicator==0){
      /*do the iterations*/
      alambda=0.001;
      i=0;
      do{
        indicator=mrqminMine(x,temp,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
        i++;
        if((i>MAX_ITER)||indicator){
          minErr*=10.0;
          if(minErr>=maxErr){
            //writeData(x,temp,numb,chisq,"failedReturn");
            found=0;
            break;
          }
        }
      }while(chisq>minErr);

      /*final call to get bits and pieces*/
      alambda=0.0;
      indicator=mrqminMine(x,temp,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
    }else found=0;

    //if(found)writeData(x,temp,numb,-1.0,"foundReturn");

    /*force widths to be positive*/
    for(i=0;i<nGauss;i++)params[3*i+3]=fabs(params[3*i+3]);
    TTIDY((void **)covar,nParams+1);
    TTIDY((void **)alpha,nParams+1);
    TIDY(pSwitch);
    TIDY(sig);
  //}/*number of points check*/

  if(found){
    *range=params[1];
    *integral=params[2]*sqrt(2.0*M_PI)*params[3];
    *width=params[3];
  }else{
    *range=*integral=*width=-1.0;
    indicator=1;
  }

  TIDY(x);
  TIDY(temp);
  return(indicator);
}/*fitFunction*/


/*###############################################*/
/*write out failed returns*/

void writeData(float *x,float *temp,int numb,float chisq,char *root)
{
  int i=0;
  char errClass=0;
  char namen[200];
  FILE *opoo=NULL;

  if(chisq<0.0)         errClass=numb;
  if(chisq<1.0)         errClass=0;
  else if(chisq<10.0)   errClass=1;
  else if(chisq<1000.0) errClass=2;
  else if(chisq<10000.0)errClass=3;
  else                  errClass=4;

  sprintf(namen,"%s.%d.n.%d.dat",root,errClass,failNumb);

  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  for(i=1;i<=numb;i++){
    fprintf(opoo,"%f %f\n",x[i],temp[i]);
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  failNumb++;
  fprintf(stderr,"Failed to %s\n",namen);

  return;
}/*writeData*/


/*###############################################*/
/*error function*/

void gaussErr(float x,float *params,float *yfit,float *dyda,int nParams)
{
  int i=0;
  int nGauss=0;
  float A=0,sig=0,mu=0;
  float arg=0,sqrt2pi=0;
  float gauss(float,float,float);

  nGauss=nParams/3;
  sqrt2pi=sqrt(2.0*M_PI);

  *yfit=0.0;
  for(i=1;i<=nParams;i++)dyda[i]=0.0;

  for(i=0;i<nGauss;i++){  /*loop over all Gaussians and add up*/
    mu=params[i*3+1];
    A=params[i*3+2];
    sig=params[i*3+3];
    arg=-1.0*(x-mu)*(x-mu)/(2.0*sig*sig);

    *yfit+=A*gauss(x,sig,mu);
    dyda[1]+=A*(x-mu)/(sig*sig)*exp(arg);
    dyda[2]+=exp(arg);
    dyda[3]+=A*(x-mu)*(x-mu)/(sig*sig*sig)*exp(arg);
  }/*Gaussian loop*/

  return;
}/*gaussian fit error*/


/*###############################################*/
/*non normalised Gaussian*/

float gauss(float x,float sigma,float offset)
{
  float y=0;
  y=exp(-1.0*(x-offset)*(x-offset)/(2.0*sigma*sigma));
  return(y);
}/*gauss*/


/*###############################################*/
/*initial parameter guess*/

float *initialGuess(float *x,float *y,int numb,int nGauss)
{
  int i=0,maxBin=0;
  float *params=NULL;
  float contN=0,thresh=0;

  params=falloc(3*nGauss+1,"parameters",0);

  contN=0.0;
  params[1]=0.0;
  params[2]=-10.0;
  for(i=1;i<=numb;i++){
    params[1]+=x[i]*y[i];  /*position is weighted by intensity*/
    contN+=y[i];
    if(y[i]>params[2]){
      params[2]=y[i];
      maxBin=i;           /*to calculate the width from */
    }
  }
  params[1]/=contN;       /*normalise mean range*/

  thresh=params[2]*exp(-0.5);  /*1/e^2 of maximum*/
  params[3]=-1.0;              /*nonsense value*/
  for(i=maxBin;i<=numb;i++){
    if(y[i]<=thresh){
      params[3]=(x[i]-params[1])/2.0;
      break;
    }
  }

  /*to prevent daft values*/
  if(params[3]<=0.00001)params[3]=1.0;
  return(params);
}/*initialGuess*/


/*##########################################*/
/*arrange data for fitting*/

float *filterData(char *data,int numb)
{
  int i=0;
  float *temp=NULL;

  temp=falloc(numb+1,"temporary data",0);
  for(i=0;i<numb;i++){
    if(data[i]>thresh)temp[i+1]=(float)data[i]-(float)backSig;
    else              temp[i+1]=0.0;
  }

  return(temp);
}/*filterData*/


/*##########################################*/
/*find outgoing pulse and check saturation*/

float findStart(int start,int end,char *satTest,char *data,int offset)
{
  int i=0,b=0,e=0;  /*loop control and bounds*/
  int place=0;      /*array index*/
  int indicator=0;  /*function fitting success indicator*/
  int fitFunction(char *,int,float *,float *,float *,float);
  float waveStart=0;  /*outgoing pulse bin*/
  float contN=0;      /*contribution weighting counter*/
  float dummy1=0,dummy2=0;  /*to save re-writing fitFunction()*/
  char satThresh=0; /*saturation threshold*/
  char signal=0;
  char started=0;

  waveStart=-1000.0; /*daft value*/
  satThresh=127;

  *satTest=0;  /*not saturated by default*/
  b=start-50;   /*4.5m, from histograms*/
  if(b<0)b=0;  /*truncate at 0*/
  e=start+sOffset;
  if(e>end)e=end;

  /*is there a signal?*/
  signal=0;
  for(i=b;i<=e;i++){
    if(data[i]>thresh){
      signal=1;
      break;
    }
  }

  if(signal){ /*fit a function to get the start range, if there is any signal*/
    indicator=fitFunction(&(data[b]),e-b+1,&waveStart,&dummy1,&dummy2,0.0);
    if(indicator){  /*function fitting failed, backup*/
      waveStart=contN=0.0;
      started=0;
      for(i=b;i<=e;i++){
        if(data[i]>thresh){
          started=1;
//fprintf(stdout,"%d %d %f\n",i,b,(float)(i-b)
          waveStart+=(float)(i-start)*((float)data[i]-(float)backSig);  /*scale by resolution later*/
          contN+=(float)data[i]-(float)backSig;
        }else if(started)break;
      }
      waveStart*=0.15/contN;
    }/*backup algorithm*/
  }

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
/*mark a point down*/

void markPoint(float zen,float az,float r,FILE *opoo,char max,int nIn,int azInd,int zenInd,int waveStart,float integral,char satTest,float rho,float width,control *options)
{
  float x=0,y=0,z=0;  /*cartesian coordinates*/
  float offset=0;     /*for the wonky beam correction*/
  float cAz=0,cZen=0; /*for the wonky beam correction*/

  if(options->oldSquint==0){
    /*squint(&(cZen),&(cAz),zen,az,options->zenSquint,options->azSquint,options->omega);*/
    cZen=options->zen[zenInd];
    cAz=options->azOff[zenInd]+az;
  }else{
    /*the 1.6 deg wonky beam correction*/
    /*we could pre-calculate the offset and trig values*/
    /*to save time. It would need storing in a structure*/
    offset=options->azSquint;
    cAz=az+atan2(tan(offset),sin(zen));
    cZen=acos(cos(offset)*cos(zen));
  }

  x=r*sin(cZen)*cos(cAz);
  y=r*sin(cZen)*sin(cAz);
  z=r*cos(cZen);   /*I don't know why this should be negative, unless my offset is wrong. DO a trip*/

  fprintf(opoo,"%g %g %g %d %g %g %g %d %d %d %d %f %d %f %f\n",x,y,z,max,cZen*180.0/M_PI,cAz*180.0/M_PI,r,nIn,azInd,zenInd,waveStart,integral,satTest,rho,width);

  return;
}/*markPoint*/


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
    if(options->joy)fprintf(stdout,"%d of %d Reading %s\n",i*options->coarsen+j,options->nAz*options->coarsen,namen);  /*progress indicator*/

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
    fprintf(opoo[band],"# 1 x, 2 y, 3 z, 4 int, 5 zen, 6 az, 7 range, 8 number, 9 azInd, 10 zenInd, 11 wave bin start, 12 energy, 13 sat flag, 14 rho, 15 width\n");
  }
  return(opoo);
}/*openOutput*/


/*##########################################*/
/*open a single output file*/

FILE **openDualOutput(char *outRoot)
{
  FILE **opoo=NULL;

  if(!(opoo=(FILE **)calloc(1,sizeof(FILE *)))){
    fprintf(stderr,"error in file pointer.\n");
    exit(1);
  }
  if((opoo[0]=fopen(outRoot,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outRoot);
    exit(1);
  }
  fprintf(opoo[0],"# 1 x, 2 y, 3 z, 4 DN 1545, 5 DN 1064, 6 satFlag 1545, 7 satFlag 1064, 8 rho 1545, 9 rho 1064, 10 range difference, 11 width 1545, 12 width 1064, 13 zenInd, 14 azInd\n");

  return(opoo);
}/*openDualOutput*/


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


/*##########################################*/
/*write a single waveform*/

void writeOne(char *data,int i,int nBins,float az,int width)
{
  int j=0;
  char namen[200];
  FILE *opoo=NULL;

  sprintf(&(namen[0]),"wideBeam.az.%d.zen.%d.width.%d.dat",(int)(az/0.001047198),i,width);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  for(j=0;j<nBins;j++)fprintf(opoo,"%f %d\n",((float)j+0.5)*0.15,data[i*nBins+j]);
  fprintf(stdout,"Wide beam to %s\n",namen);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  return;
}/*writeOne*/


/*##########################################*/
/*read calibration LUT*/

calibration *readCalibration(char *calFile,float filt,float filt1,float filt2)
{
  int band=0,i=0;
  calibration *cal=NULL;
  char line[200];
  char temp1[100],temp2[100],temp3[100];
  void setTransmission(calibration *,float);
  FILE *ipoo=NULL;

  if(!(cal=(calibration *)calloc(2,sizeof(calibration)))){
    fprintf(stderr,"error in calibration allocation.\n");
    exit(1);
  }
  cal[0].minDN=cal[1].minDN=0;
  cal[0].maxDN=cal[1].maxDN=4000;
  cal[0].numb=cal[1].numb=cal[0].maxDN-cal[0].minDN+1;  /*hard wired for now*/

  if(filt1<0.0)setTransmission(cal,filt);
  else{
    cal[0].tran=filt1;
    cal[1].tran=filt2;
  }

  for(band=0;band<2;band++)cal[band].LUT=falloc(cal[band].numb,"calibration LUT",band);

  if((ipoo=fopen(calFile,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",calFile);
    exit(1);
  }
  while(fgets(line,200,ipoo)!=NULL){
    if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
      if(strncasecmp(temp1,"#",1)){
        i=atoi(temp1);
        if((i<0)||(i>=cal[0].numb)){
          fprintf(stderr,"Stepped too far\n");
          exit(1);
        }
        cal[0].LUT[i]=atof(temp2)/(cal[0].tran*cal[0].tran);  /*scale by filter transmission*/
        cal[1].LUT[i]=atof(temp3)/(cal[1].tran*cal[1].tran);
      }
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(cal);
}/*readCalibration*/


/*##########################################*/
/*set filter transmission*/

void setTransmission(calibration *cal,float filt)
{
  float tol=0;

  tol=0.0001;

  /*these were found with the ASD and are contained in*/
  /*/mnt/geodesy38/nsh103/ASD/filter/filterStrength.dat*/
  /* or */
  /*/mnt/geodesy38/nsh103/ASD/filter/rachelFilter.dat*/
  if(fabs(filt-0.0)<tol){
    cal[0].tran=cal[1].tran=1.0;
  }else if(fabs(filt-0.6)<tol){
    cal[0].tran=0.47958;   /*Rachel's*/
    cal[1].tran=0.29816;
    /*cal[0].tran=0.500667;
    cal[1].tran=0.316649;*/
  }else if(fabs(filt-1.0)<tol){
    cal[0].tran=0.217715;   /*Rachel's*/
    cal[1].tran=0.090333;
    /*cal[0].tran=0.333872;
    cal[1].tran=0.124325;*/
  }else if(fabs(filt-1.6)<tol){
    cal[0].tran=0.104413;   /*Rachel's*/
    cal[1].tran=0.026934;
    /*cal[0].tran=0.167159;
    cal[1].tran=0.0393674;*/
  }else{
    fprintf(stderr,"Don't know transmission for %f\n",filt);
    exit(1);
  }
  return;
}/*setTransmission*/


/*##########################################*/
/*tidy calibration arrays*/

void tidyCal(calibration *cal)
{
  int i=0;

  for(i=0;i<2;i++){
    TIDY(cal[i].LUT);
  }
  TIDY(cal);
  return;
}/*tidyCal*/


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

  //if(zM<0.0){  /*to keep consistent with SALCA input*/
  //  *cZen*=-1.0;
  //  *cAz+=M_PI;
  //}  /*to keep consistent with SALCA input*/

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
}/*rotateY*/


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
}/*rotateZ*/


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


/*the end*/
/*##########################################*/

