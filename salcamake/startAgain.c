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


/*function to optimise to*/
void gaussErr(float,float *,float *,float *,int);

#define MAX_ITER 30
#define TOL 0.00001

int sOffset;   /*wavestart offset*/
float backSig;  /*background signal*/
char thresh;   /*noise threshold*/
char satThresh; /*saturation threshold*/

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
  char calibrate;   /*calibration switch*/
  int nAz;
  int nZen;
  float azStep;
  float zStep;
  float maxZen;
  float azStart;
  int zenOffset[2]; /*zenith offset*/
  float azSquint;
  float zenSquint;
  float mSquint;    /*motor squint angle*/
  float *zen;       /*true zenith*/
  float *azOff;     /*azimuth offset*/
  float omega;      /*mirror slope angle*/
  float res;        /*range resolution*/
  char joy;
  float *smoother;    /*smoothing function quad*/
  float sWidth;       /*smoothing function width quad*/
  int nSmoo;          /*length of smoothing function quad*/
  float *smooLM;      /*smoothing function LM*/
  float lmWidth;      /*smoothing function width LM*/
  int nSmooLM;        /*length of smoothing function LM*/
  char doLM;          /*do LM Gaussian fitting flag*/
  char doQuad;        /*do Jupp's quadratic fitting*/
}control;


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
  void pointOut(char *,int,int,FILE *,float,int,int,int,char,calibration *,control *,int,int);
  void closeFiles(FILE **,char *);
  void tidyCal(calibration *);
  calibration *readCalibration(char *,float,float,float);
  calibration *cal=NULL;
  control *options=NULL;
  FILE **openOutput(char *);
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
  backSig=-116.0;
  thresh=-110;
  satThresh=125;      /*ringing thresholdf*/

  if(!(options=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  options->coarsen=1;   /*use at native resolution*/
  options->calibrate=0;/*output DN rather than reflectance*/
  options->azStep=0.06;   /*in degrees*/
  options->zenOffset[0]=options->zenOffset[1]=0;
  options->maxZen=190.0;   /*full hemispheric scan*/
  options->azSquint=1.6*M_PI/180.0;
  options->zenSquint=0.0;
  options->omega=M_PI/4.0;  /*45 degrees*/
  options->mSquint=0.0;
  options->azStart=0.0;
  options->nZen=3200;
  options->nAz=666;
  options->joy=1;        /*print out status*/
  options->res=0.15;
  options->smoother=NULL;
  options->sWidth=1.1;
  options->nSmoo=0;
  options->smooLM=NULL;
  options->lmWidth=0.0;
  options->nSmooLM=0;
  options->doLM=1;    /*do LM Gaussian fitting*/
  options->doQuad=1;  /*do Jupp's quadratic fitting*/

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
      }else if(!strncasecmp(argv[i],"-coarsen",8)){
        checkArguments(1,i,argc,"-coarsen");
        options->coarsen=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-zenOffset",10)){
        options->zenOffset[1]=1;
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
      }else if(!strncasecmp(argv[i],"-mSquint",8)){
        checkArguments(1,i,argc,"-mSquint");
        options->mSquint=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-background",11)){
        checkArguments(1,i,argc,"-background");
        backSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-joyless",8)){
        options->joy=0;
      }else if(!strncasecmp(argv[i],"-sWidth",7)){
        checkArguments(1,i,argc,"-sWidth");
        options->sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lmWidth",8)){
        checkArguments(1,i,argc,"-lmWidth");
        options->lmWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noLM",5)){
        options->doLM=0;
      }else if(!strncasecmp(argv[i],"-noQuad",7)){
        options->doQuad=0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-outRoot root;  output filename root\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in degrees\n-maxR range;    maximum recorded range, metres\n-calibrate;     calibrate to reflectance, need calibration file\n-calFile name;  calibrate to reflectance, need calibration file\n-filt n;        filter optical depth (0, 0.6, 1 or 1.6)\n-nFilt t1 t2;   specify transmissions rather than use defaults\n-coarsen n;     coarsen by a factor. CAUTIION, intensities will not scale properly yet\n-zenOffset;     offset 1064nm zenith by 1\n-maxZen zen;    maximum zenith angle, degrees. 190 by default\n-azSquint angle;    azimuth squint angle in degrees\n-zenSquint angle;   zenith squint angle, degrees\n-mSquint angle;     mirror squint angle in degrees\n-omega angle;       mirror angle in degrees\n-azStart angle;     azimuth start, degrees\n-background DN;     background DN value\n-joyless;           don't print out status\n-sWidth width;     smoothing width for quadratic fitting\n-lmWidth width;       smoothing for LM\n-noLM;     don't to LM Gaussian fitting\n-noQuad;       don't do quadratic fitting\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/



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
  opoo=openOutput(outRoot);


  /*if coarsened*/
  if(options->coarsen<1)options->coarsen=1;  /*so zero means uncoarsened*/
  options->nAz/=options->coarsen;
  numb/=options->coarsen;
  options->azStep*=(float)options->coarsen*M_PI/180.0; /*convert to radians as well*/
  options->nZen/=options->coarsen;
  options->zStep=0.001047198*(float)options->coarsen;  /*fixed for SALCA*/

  /*set up squint angles*/
  translateSquint(options);
  setSquint(options,numb);


  /*loop through azimuth files and write out points*/
  for(i=0;i<=options->nAz;i++){ /*az loop*/
    data=readData(inRoot,i,&numb,&nBins,&length,options); /*read binary data into an array*/

    if(data){ /*does the file exist?*/

      for(band=0;band<2;band++){   /*loop through wavebands*/
        pointOut(data,numb,nBins,opoo[band],(float)i*options->azStep+options->azStart,start[band],end[band],i,options->calibrate,&(cal[band]),options,options->zenOffset[band],band);
      }/*band loop*/
      TIDY(data);   /*clean up arrays as we go along*/
    }/*does file exist?*/
  }/*az loop*/

  closeFiles(opoo,outRoot);   /*close output files*/

  if(options){
    if(options->calibrate)tidyCal(cal);
    TIDY(options->azOff);
    TIDY(options->zen);
    TIDY(options->smoother);
    TIDY(options->smooLM);
    TIDY(options);
  }
  return(0);
}/*main*/


/*##########################################*/
/*output a point cloud*/

void pointOut(char *data,int numb,int nBins,FILE *opoo,float az,int start,int end,int azInd,char calibrate,calibration *cal,control *options,int zOff,int band)
{
  int i=0,j=0;  /*loop variables*/
  int jS=0;
  int place=0;  /*array place*/
  int maxPlace=0; /*for tracking the smoothed wave*/
  int nIn=0;    /*number of features*/
  int waveStart=0; /*outgoing pulse bin*/
  int width=0;      
  int findStart(int,int,char *,int);
  float findQstart(int,int,char *,int,control *);
  float findSstart(int,int,char *,int,float,float *);
  float range=0;         /*peak range*/
  float sRange=0;         /*mean range weighted by intensity*/
  float qInt=0,qRange=0;  /*quadratic fit parameters*/
  float qIntegral=0;      /*quadratic integral*/
  float lmInt=0,lmRange=0;  /*L-M fit parameters*/
  float lmSig=0;            /*gaussian width*/
  float integral=0.0;     /*return energy integral*/
  float rho=0;            /*target reflectance*/
  float res=0;
  float error=0;
  float qStart=0,sStart=0;   /*signal starts found from quadrativ and sum methods*/
  char max=0;    /*intensity threshold and max value*/
  char satTest=0;         /*saturation test indicator*/
  char ringTest=0;        /*test for ringing*/
  char checkRinging(int,int,char *,int);
  char brEak=0;           /*break indicator*/
  void markPoint(float,float,FILE *,char,int,int,int,int,float,char,float,int,control *,float,float,float,float,float,float,float,float,float);
  void fitFunction(char *,int,float *,float *,float *,float,int,int,float *,control *);
  void fitQuad(char *,int,float *,float *,float,control *,float *);

  float sCount=0;

  res=options->res;
  rho=-1.0;           /*nonesense value*/

  qRange=qInt=qIntegral=-1.0; /*set to nonesense*/
  error=lmRange=lmInt=-1.0;   /*set to nonesense*/

  for(i=0;i<numb;i++){ /*zenith loop*/
    nIn=1;                      /*number of points per waveform*/

    /*find the waveform start and check for saturation*/
    /*peak*/
    waveStart=findStart(start,end,data,i*nBins);
    if((band==0)&&(waveStart<0))waveStart=0; /*as the 1545 pulse is often lost off the end*/
    /*quadratic*/
    /*if(options->doQuad)qStart=findQstart(start,end,data,i*nBins,options);*/
    /*sum*/
    sStart=findSstart(start,end,data,i*nBins,options->res,&sCount);

    ringTest=checkRinging(start,end,data,i*nBins);

    if(waveStart>=0){
      for(j=start;j<end;j++){     /*loop along waveform*/
        place=i*nBins+j;            /*array place*/
        if(data[place]>thresh){   /*signal above noise level*/
          satTest=0;
          sRange=0.0;
          width=0;
          max=thresh;          /*reset max*/
          integral=0;

          jS=j;                   /*mark data for function fitting*/
          for(;j<end;j++){        /*step to end of feature*/
            place=i*nBins+j;        /*array place*/
            integral+=(float)data[place]-backSig;
            sRange+=(((float)(j-start)+1.5)*options->res-sStart)*((float)data[place]-backSig);
            if(data[place]>max){  /*record maximum intensity and position*/
              max=data[place];      /*set max*/
              maxPlace=place;       /*record array point for later tracking*/
              range=((float)(j-waveStart)+1.5)*0.15;  /*use outgoing pulse to set start*/  /*why +1.5?*/
            }                     /*max test*/

            if(data[place]==127)satTest=1;

            if(data[place]<=thresh){ /*we have left feature*/
              if(integral>0.0)sRange/=integral;
              else            sRange=range;

              /*fit quadratic*/
              if(options->doQuad)fitQuad(&(data[i*nBins+jS-1]),j-jS+1,&qRange,&qInt,(float)(jS-waveStart)*res,options,&qIntegral);

              /*fit Gaussian*/
              if(options->doLM)fitFunction(&(data[i*nBins+jS-1]),j-jS+1,&lmRange,&lmInt,&lmSig,(float)(jS-waveStart)*res,band,azInd,&error,options);


              /*if(calibrate){
                if(((int)integral>=0)&&((int)integral<cal->numb))rho=cal->LUT[(int)integral]*range*range;
                else{
                  if((int)integral<=0)rho=0.0;
                  else if((int)integral>=cal->numb)rho=1.0;
                }
              }*/

              markPoint(az,range,opoo,max,nIn,azInd,i+zOff,waveStart,integral,satTest,rho,width,options,qRange,qInt,lmRange,lmInt,error,qIntegral,sRange,sStart,sCount);
              width=0;
              if(ringTest&&(max>=satThresh))brEak=1;   /*if saturated only take the first return*/
              nIn++;    /*record number of points per waveform*/
              break;
            }/*left feature*/
            width++;
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
/*check for ringing*/

char checkRinging(int start,int end,char *data,int offset)
{
  int i=0,place=0;
  char satTest=0;

  /*test for saturation*/
  satTest=0;
  for(i=start;i<end;i++){  /*loop through full waveform to test for saturation*/
    place=offset+i;
    if(data[place]>=satThresh){   /*saturated*/
      satTest=1;
      break;
    }
  }/*saturation test loop*/

  return(satTest);
}/*checkRinging*/


/*##########################################*/
/*Prepare for David Jupp's quadratic fitting*/

void fitQuad(char *data,int length,float *qRange,float *qInt,float sRange,control *options,float *qIntegral)
{
  int i=0;
  int nLength=0;       /*smoothed length*/
  float *temp=NULL;    /*smoothed function to fit to*/
  float *setSmoother(float,float,int *);
  float *smooth(char *,int,float *,int,float *);
  void fitQuadratic(float *,float *,int,float *,float *,float *);
  float *x=NULL;


  /*if not already set up, set up the smoothing function*/
  if(!options->smoother)options->smoother=setSmoother(options->res,options->sWidth,&(options->nSmoo));

  /*smooth relevant bit of waveform*/
  nLength=length+options->nSmoo;
  x=falloc(nLength,"range",0);
  for(i=0;i<nLength;i++)x[i]=(float)((i-options->nSmoo/2)+0.5)*options->res+sRange;

  temp=smooth(&(data[-1*options->nSmoo/2]),nLength,options->smoother,options->nSmoo,x);

  fitQuadratic(x,temp,nLength,qRange,qInt,qIntegral);

  TIDY(x);
  TIDY(temp);
  return;
}/*fitQuad*/


/*###############################################*/
/*fit a quadratic using David Jupp's method*/

void fitQuadratic(float *x,float *y,int width,float *qRange,float *qInt,float *qIntegral)
{
  int i=0;
  int mid=0;
  float max=0;
  float a=0,b=0,c=0;  /*quadratic coefficients*/
  float d0=0,d1=0;    /*gradients*/

  if(width<1){
    fprintf(stderr,"Not enough points, need to pad\n");
    exit(1);
  }

  /*find three brightest points to fit to*/
  max=-100.0;
  for(i=0;i<width;i++){
    if(y[i]>max){
      max=y[i];
      mid=i;
    }
  }

  d0=(y[mid]-y[mid-1])/(x[mid]-x[mid-1]);
  d1=(y[mid+1]-y[mid])/(x[mid+1]-x[mid]);

  c=(d1-d0)/(x[mid+1]-x[mid-1]);
  b=d0-c*(x[mid]+x[mid-1]);
  a=y[mid]-b*x[mid]-c*x[mid]*x[mid];

  *qRange=-1.0*b/(2.0*c);
  /*enforce bounds. Do not let the range stray beyond these three points*/
  if((*qRange)<x[mid-1])     *qRange=x[mid-1];
  else if((*qRange)>x[mid+1])*qRange=x[mid+1];

  *qInt=a+b*(*qRange)+c*(*qRange)*(*qRange);
  *qIntegral=2.0*sqrt(b*b-4.0*c*(a-backSig));   /*integral above backSig*/

  return;
}/*fitQuadratic*/


/*###############################################*/
/*smooth a waveform*/

float *smooth(char *y,int numb,float *smoother,int nSmoo,float *x)
{
  int i=0,j=0;
  int place=0;
  float contN=0;
  float *smoothY=NULL;

  smoothY=falloc(numb,"smoothed",0);

  if(nSmoo>0){
    for(i=0;i<numb;i++){
      smoothY[i]=0.0;
      contN=0.0;
      for(j=0;j<nSmoo;j++){
        place=i+j-nSmoo/2;
        if((place>=0)&&(place<numb)){  /*check we're within the arrays*/
          if(x[place]>0.0){            /*check we've not stepped off the signal*/
            if(y[place]>thresh){
              smoothY[i]+=((float)y[place]-backSig)*smoother[j];
            }
            contN+=smoother[j];
          }
        }
      }
      if(contN>0.0)smoothY[i]/=contN;
    }
  }else{ /*no smoothing, copy raw*/
    for(i=0;i<numb;i++){
      smoothY[i]=(float)y[i]-backSig;
    }
  }
  return(smoothY);
}/*smooth*/


/*###############################################*/
/*set up a smoothing function*/

float *setSmoother(float sRes,float sWidth,int *nSmoo)
{
  int i=0;
  float tol=0,mid=0,y=0;
  float *smoother=NULL;

  tol=0.0001;

  if(sWidth<tol)*nSmoo=0;   /*no smoothing*/
  else{
    /*determine width*/
    i=0;
    do{
      y=gaussian((float)i*sRes,sWidth,0.0);
      i++;
    }while(y>=tol);
    *nSmoo=i*2;

    mid=(float)i*sRes;
    smoother=falloc(*nSmoo,"smoother",0);
    for(i=0;i<*nSmoo;i++){
      smoother[i]=gaussian((float)i*sRes,sWidth,mid);
    }
  }

  return(smoother);
}/*setSmoother*/


/*##########################################*/
/*fit a a Gaussian by Levenberg-Marquardt*/

void fitFunction(char *data,int numb,float *range,float *integral,float *width,float rOffset,int band,int azInd,float *error,control *options)
{
  int i=0;
  int nGauss=0,nParams=0;
  int *pSwitch=NULL;   /*switch to optimise each parameter*/
  float *x=NULL;
  float *temp=NULL;
  float *filterData(char *,int);
  float chisq=0;
  float minErr=0;
  float alambda=0;
  float **covar=NULL;
  float **alpha=NULL;
  float *sig=NULL;
  float *params=NULL;
  float maxErr=0;
  float *initialGuess(float *,float *,int,int);
  float *smooth(char *,int,float *,int,float *);
  int mrqminMine(float *,float *,float *,int,float *,int *,int,float **,float **,float *,void (*)(float,float *,float *,float *,int),float *);
  int indicator=0;
  char found=1;

  /*tolerances*/
  minErr=0.001;
  maxErr=10.0;


  /*if not already set up, set up the smoothing function*/
  if(!options->smooLM)options->smooLM=setSmoother(options->res,options->lmWidth,&(options->nSmooLM));

  /*smooth relevant bit of waveform*/
  numb+=options->nSmooLM;
  x=falloc(numb,"range",0);
  for(i=0;i<numb;i++)x[i]=(float)(i-options->nSmooLM/2)*options->res+rOffset;

  temp=smooth(&(data[-1*options->nSmooLM/2]),numb,options->smooLM,options->nSmooLM,x);


  /*number of Gaussians to fit*/
  nGauss=1;
  nParams=3*nGauss;
  params=initialGuess(x,temp,numb,nGauss);

  chisq=-1.0;   /*nonesense value*/
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


  /*numb is the number of return bins above noise in this feature*/
  if(numb<=3){    /*one point above noise, fix max and position*/
    pSwitch[1]=pSwitch[2]=0;
    /*leave params[1] and params[2] as the initial guesses*/
  }else if(numb<=4){  /*two points above noise, fix width*/
    params[3]=0.11;
    pSwitch[3]=0;
  }/*number of points check*/

  /*initialise*/
  alambda=-1.0;
  indicator=mrqminMine(x,temp,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
  if(indicator==0){
    /*do the iterations*/
    alambda=0.001;
    i=0;
    do{
      indicator=mrqminMine(x,temp,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
      if(indicator)break;
      i++;
      if(i>MAX_ITER){
        minErr*=10.0;
        if(minErr>=maxErr){
          found=0;
          break;
        }
      }
    }while(chisq>minErr);

    if(indicator==0){
      /*final call to get bits and pieces*/
      alambda=0.0;
      indicator=mrqminMine(x,temp,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
    }

    /*force widths to be positive*/
    for(i=0;i<nGauss;i++)params[3*i+3]=fabs(params[3*i+3]);
  }/*gaussj work check*/
  TTIDY((void **)covar,nParams+1);
  TTIDY((void **)alpha,nParams+1);
  TIDY(pSwitch);
  TIDY(sig);

  if(indicator==0){
    *range=params[1];
    *integral=params[2]*sqrt(2.0*M_PI)*params[3];
    *width=params[3];

  }else{
    *range=*integral=*width=-1.0;
  }

  *error=chisq;

  TIDY(x);
  TIDY(temp);
  return;
}/*fitFunction*/


/*###############################################*/
/*error function*/

void gaussErr(float x,float *params,float *yfit,float *dyda,int nParams)
{
  int i=0;
  int nGauss=0;
  float A=0,sig=0,mu=0;
  float arg=0,sqrt2pi=0;

  nGauss=nParams/3;
  sqrt2pi=sqrt(2.0*M_PI);

  *yfit=0.0;
  for(i=1;i<=nParams;i++)dyda[i]=0.0;

  for(i=0;i<nGauss;i++){  /*loop over all Gaussians and add up*/
    mu=params[i*3+1];
    A=params[i*3+2];
    sig=params[i*3+3];
    arg=-1.0*(x-mu)*(x-mu)/(2.0*sig*sig);

    *yfit+=A*gaussian(x,sig,mu);
    dyda[1]+=A*(x-mu)/(sig*sig)*exp(arg);
    dyda[2]+=exp(arg);
    dyda[3]+=A*(x-mu)*(x-mu)/(sig*sig*sig)*exp(arg);
  }/*Gaussian loop*/

  return;
}/*gaussian fit error*/


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
/*find outgoing pulse start using quadratic*/

float findQstart(int start,int end,char *data,int offset,control *options)
{
  int i=0,place=0;
  int b=0,e=0;
  int iS=0;
  float qStart=0;
  float qIntegral=0,qInt=0;   /*dummy variables*/
  char brEak=0;
  void fitQuad(char *,int,float *,float *,float,control *,float *);


  /*find feature*/
  b=start-50;   /*4.5m, from histograms*/
  if(b<0)b=0;  /*truncate at 0*/
  e=start+sOffset;   /*sOffset is a global variable*/
  if(e>end)e=end;
  brEak=0;

  iS=-1;
  for(i=b;i<e;i++){  /*loop around where we think the pulse might be*/
    place=offset+i;
    if(data[place]>thresh){ /*feature start*/
      iS=i;
      for(;i<e;i++){/*step to end of feature*/
        place=offset+i;
        if(data[place]<thresh)break;
      }
    }/*feature start*/
    if(brEak)break;
  }/*find feature loop*/


  if(iS>0){
    /*perform quadratic fitting*/
    fitQuad(&(data[offset+iS-1]),i-iS+1,&qStart,&qInt,0.0,options,&qIntegral);
  }else qStart=0.0;

  return(qStart);
}/*findQstart*/


/*##########################################*/
/*find outgoing pulse start using weighted mean*/

float findSstart(int start,int end,char *data,int offset,float res,float *sCount)
{
  int i=0,place=0;
  int b=0,e=0;
  float sStart=0;
  float contN=0;
  char brEak=0;

  /*find feature*/
  b=start-50;   /*4.5m, from histograms*/
  if(b<0)b=0;  /*truncate at 0*/
  e=start+sOffset;   /*sOffset is a global variable*/
  if(e>end)e=end;
  brEak=0;

  sStart=contN=0.0;
  for(i=b;i<e;i++){  /*loop around where we think the pulse might be*/
    place=offset+i;
    if(data[place]>thresh){ /*feature start*/
      for(;i<e;i++){
        place=offset+i;
        if(data[place]<thresh)break;
        sStart+=(float)(i-start)*res*((float)data[place]-backSig);
        contN+=(float)data[place]-backSig;
      }/*loop along feature*/
      brEak=1;
    }/*feature found*/
    if(brEak)break;
  }/*range loop*/

  if(contN>0.0)sStart/=contN;
  else         sStart=-1.0*(float)sOffset*res;  /*note start has an offset on it*/

  *sCount=contN;

  return(sStart);
}/*findSstart*/


/*##########################################*/
/*find outgoing pulse and check saturation*/

int findStart(int start,int end,char *data,int offset)
{
  int i=0,b=0,e=0;  /*loop control and bounds*/
  int place=0;      /*array index*/
  int waveStart=0;  /*outgoing pulse bin*/
  char max=0;       /*max intensity*/

  waveStart=-1; /*start-sOffset;*/

  max=-125;

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

  return(waveStart);
}/*findStart*/


/*##########################################*/
/*mark a point down*/

void markPoint(float az,float r,FILE *opoo,char max,int nIn,int azInd,int zenInd,int waveStart,float integral,char satTest,float rho,int width,control *options,float qRange,float qInt,float lmRange,float lmInt,float error,float qIntegral,float sRange,float sStart,float qStart)
{
  float x=0,y=0,z=0;  /*cartesian coordinates*/
  float cAz=0,cZen=0; /*for the wonky beam correction*/

  cZen=options->zen[zenInd];
  cAz=options->azOff[zenInd]+az;

  x=r*sin(cZen)*cos(cAz);
  y=r*sin(cZen)*sin(cAz);
  z=r*cos(cZen);

  fprintf(opoo,"%g %g %g %d %g %g %g %d %d %d %d %f %d %f %d %f %f %f %f %f %f %f %f %f\n",x,y,z,max,cZen*180.0/M_PI,cAz*180.0/M_PI,r,nIn,azInd,zenInd,waveStart,integral,satTest,rho,width,sRange,qRange,qInt,lmRange,lmInt,error,qIntegral,sStart,qStart);

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
    sprintf(outNamen,"%s.band.%d.shape",outRoot,band);
    if((opoo[band]=fopen(outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",outNamen);
      exit(1);
    }
    fprintf(opoo[band],"# 1 x, 2 y, 3 z, 4 int, 5 zen, 6 az, 7 peak range, 8 number, 9 azInd, 10 zenInd, 11 wave bin start, 12 energy, 13 sat flag, 14 rho, 15 width, 16 sum range, 17 quadRange, 18 quadInt, 19 lmRange, 20 lmInt, 21 lmError, 22 quad integral, 23 sum start, 24 quad start\n");
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
    fprintf(stdout,"Written to %s.band.%d.shape\n",outRoot,band);
  }
  TIDY(opoo);

  return;
}/*closeFiles*/


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
    zen=((float)(options->nZen/2)-(float)i)*options->zStep+options->mSquint;
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

  vect=falloc(3,"vector",0);
  vect[0]=lX;
  vect[1]=lY;
  vect[2]=lZ;

  thetaX=-1.0*atan2(sqrt(rX*rX+rY*rY),rZ);
  thetaZ=-1.0*atan2(rX,rY);
  rotateZ(vect,thetaZ);
  rotateX(vect,thetaX);
  rotateZ(vect,-2.0*inc);
  rotateX(vect,-1.0*thetaX);
  rotateZ(vect,-1.0*thetaZ);

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


/*the end*/
/*##########################################*/

