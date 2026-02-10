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

/*function to optimise to*/
void gaussErr(float,float *,float *,float *,int);


int sOffset;   /*wavestart offset*/
char backSig;  /*background signal*/
char thresh;   /*noise threshold*/
float res;     /*range resolution*/

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
  int zenOffset[2]; /*zenith offset*/
  float azSquint;
  float zenSquint;
  float *zen;       /*true zenith*/
  float *azOff;     /*azimuth offset*/
  float omega;      /*mirror slope angle*/
  char func;        /*perform function fitting*/
  char oldSquint;
  char joy;
  float *smoother;    /*smoothing function*/
  float sWidth;       /*smoothing function width*/
  int nSmoo;          /*length of smoothing function*/
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
  void pointOut(char *,int,int,FILE *,float,int,int,int,char,calibration *,control *,int,int,char *);
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
  maxR=60.0;
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
  options->zenOffset[0]=options->zenOffset[1]=0;
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
  options->smoother=NULL;
  options->sWidth=1.0;

  res=0.15;    /*SALCA range resolution*/


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
      }else if(!strncasecmp(argv[i],"-background",11)){
        checkArguments(1,i,argc,"-background");
        backSig=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fit",4)){
        options->func=1;
      }else if(!strncasecmp(argv[i],"-oldSquint",10)){
        options->oldSquint=1;
      }else if(!strncasecmp(argv[i],"-joyless",8)){
        options->joy=0;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-outRoot root;  output filename root\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in degrees\n-maxR range;    maximum recorded range, metres\n-calibrate;     calibrate to reflectance, need calibration file\n-calFile name;  calibrate to reflectance, need calibration file\n-filt n;        filter optical depth (0, 0.6, 1 or 1.6)\n-nFilt t1 t2;   specify transmissions rather than use defaults\n-dualOut;       output a combined point cloud\n-coarsen n;     coarsen by a factor. CAUTIION, intensities will not scale properly yet\n-zenOffset;     offset 1064nm zenith by 1\n-maxZen zen;    maximum zenith angle, degrees. 190 by default\n-azSquint angle;    azimuth squint angle in degrees\n-zenSquint angle;   zenith squint angle, degrees\n-omega angle;       mirror angle in degrees\n-azStart angle;     azimuth start, degrees\n-background DN;     background DN value\n-fit;               fit a function\n-oldSquint;         use old squint function\n-joyless;           don't print out status\n\nsmoothing;\n\n");
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
  nBins=(int)((150.0+maxR)/res);  /*150m of 1550 nm plus maxR of 1040 nm*/

  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(maxR/res);
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
  if(options->oldSquint==0){
    translateSquint(options);
    setSquint(options,numb);
  }


  /*loop through azimuth files and write out points*/
  for(i=0;i<=options->nAz;i++){ /*az loop*/
    data=readData(inRoot,i,&numb,&nBins,&length,options); /*read binary data into an array*/

    if(data){ /*does the file exist?*/
      for(band=0;band<2;band++){   /*loop through wavebands*/
        pointOut(data,numb,nBins,opoo[band],(float)i*options->azStep+options->azStart,start[band],end[band],i,options->calibrate,&(cal[band]),options,options->zenOffset[band],band,outRoot);
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
    TIDY(options);
  }
  return(0);
}/*main*/


/*##########################################*/
/*output a point cloud*/

void pointOut(char *data,int numb,int nBins,FILE *opoo,float az,int start,int end,int azInd,char calibrate,calibration *cal,control *options,int zOff,int band,char *outRoot)
{
  int i=0,j=0;  /*loop variables*/
  int jS=0;
  int place=0;  /*array place*/
  int maxPlace=0; /*for tracking the smoothed wave*/
  int nIn=0;    /*number of features*/
  int findStart(int,int,char *,char *,int);
  int startPlace=0;
  int waveStart=0; /*outgoing pulse bin*/
  int sum=0;
  int length=0;           /*feature length*/
  int binMax=0;           /*bin with maximum bin*/
  float zen=0;            /*zenith and step in radians*/
  float range=0;          /*range index*/
  float qInt=0,qRange=0;  /*quadratic fit parameters*/
  float width=0;      
  char max=0;             /*intensity threshold and max value*/
  char satTest=0;         /*saturation test indicator*/
  char ringTest=0;        /*test for ringing*/
  char brEak=0;           /*break indicator*/
  char satThresh=0;       /*saturation threshold*/
  float integral=0.0;     /*return energy integral*/
  float rho=0;            /*target reflectance*/
  float error=0;
  void fitFunction(char *,int,float *,float *,float *,float,int,int,float *,char *,int);
  void writeFitResults(float,float,float,int,char,FILE *,float,int,float,float,int,int,float,float,float,int);
  void fitQuad(char *,int,float *,float *,float,control *);

  satThresh=120;      /*ringing thresholdf*/
  rho=-1.0;           /*nonesense value*/

  for(i=0;i<numb;i++){ /*zenith loop*/
    zen=((float)(options->nZen/2)-(float)(i+zOff))*options->zStep;
    nIn=1;                      /*number of points per waveform*/

    /*find the waveform start and check for saturation*/
    waveStart=findStart(start,end,&ringTest,data,i*nBins);
    if((band==0)&&(waveStart<0))waveStart=0; /*as the 1545 pulse is often lost off the end*/

    if(waveStart>=0){
      for(j=start;j<end;j++){     /*loop along waveform*/
        place=i*nBins+j;            /*array place*/
        if(data[place]>thresh){   /*signal above noise level*/
          satTest=0;
          max=thresh;          /*reset max*/

          jS=j;                   /*mark data for function fitting*/
          sum=0;               /*reset sum of bins within return*/
          for(;j<end;j++){        /*step to end of feature*/
            place=i*nBins+j;        /*array place*/
            if(data[place]>max){  /*record maximum intensity and position*/
              max=data[place];      /*set max*/
              maxPlace=place;       /*record array point for later tracking*/
              binMax=j;           /*bin of maximum intesity*/
            }                     /*max test*/

            if(data[place]==127)satTest=1;       /*check for saturation*/
            sum+=(int)data[place]-(int)backSig;  /*running total of bin intensites*/

            if(data[place]<=thresh){ /*left feature*/
              startPlace=i*nBins+jS-1;   /*start of feature*/
              length=j-jS+1;             /*end of feature*/
              if(startPlace<0){
                length+=startPlace;
                startPlace=0;
              }

              /*Levenberg-Marquardt fit*/
              fitFunction(&(data[startPlace]),length,&range,&integral,&width,(float)(jS-waveStart)*res,band,azInd,&error,outRoot,waveStart);

              /*quadratic fit*/
              fitQuad(&(data[startPlace]),length,&qRange,&qInt,(float)(jS-waveStart)*res,options);

              /*write out all results to a single file*/
              writeFitResults(range,integral,width,j-jS+1,max,opoo,error,sum,qRange,qInt,i,azInd,(float)(binMax-waveStart)*res,options->zen[i],(float)azInd*options->azStep+options->azStart+options->azOff[i],waveStart);

              width=0;    /*reset bin width counter*/
              if((ringTest)&&(max>=satThresh))brEak=1;   /*if saturated only take the first return*/
              nIn++;    /*record number of points per waveform*/
              break;
            }/*left feature*/

            width++;
          }/*feature end*/
        }/*feature start*/
        if(brEak){   /*for saturated only take first return*/
          brEak=0;   /*this avoids ringing*/
          break;
        }/*break if*/
      }/*bin loop*/
    }else{/*found wavestart check*/
fprintf(stdout,"Missed %d %d %d %d\n",i,azInd,band,waveStart);
    }
  }/*zenith loop*/
  return;
}/*pointOut*/


/*##########################################*/
/*Prepare for David Jupp's quadratic fitting*/

void fitQuad(char *data,int length,float *qRange,float *qInt,float sRange,control *options)
{
  int i=0;
  int nLength=0;       /*smoothed length*/
  float *temp=NULL;    /*smoothed function to fit to*/
  float *setSmoother(float,float,int *);
  float *smooth(char *,int,float *,int,float *);
  void fitQuadratic(float *,float *,int,float *,float *);
  float *x=NULL;


  /*if not already set up, set up the smoothing function*/
  if(!options->smoother)options->smoother=setSmoother(res,options->sWidth,&(options->nSmoo));

  /*smooth relevant bit of waveform*/
  nLength=length+options->nSmoo;
  x=falloc(nLength,"range",0);
  for(i=0;i<nLength;i++)x[i]=(float)(i-options->nSmoo/2)*res+sRange;

  temp=smooth(&(data[-1*options->nSmoo/2]),nLength,options->smoother,options->nSmoo,x);

  fitQuadratic(x,temp,nLength,qRange,qInt);

  TIDY(x);
  TIDY(temp);
  return;
}/*fitQuad*/


/*###############################################*/
/*fit a quadratic using David Jupp's method*/

void fitQuadratic(float *x,float *y,int width,float *qRange,float *qInt)
{
  int i=0;
  int mid=0;
  float max=0;
  float a=0,b=0,c=0;  /*quadratic coefficients*/
  float d0=0,d1=0;    /*gradients*/

  if(width<3){
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

  for(i=0;i<numb;i++){
    smoothY[i]=0.0;
    contN=0.0;
    for(j=0;j<nSmoo;j++){
      place=i+j-nSmoo/2;
      if((place>=0)&&(place<numb)){  /*check we're within the arrays*/
        if(x[place]>0.0){            /*check we've not stepped off the signal*/
          if(y[place]>thresh){
            smoothY[i]+=((float)y[place]-(float)backSig)*smoother[j];
          }
          contN+=smoother[j];
        }
      }
    }
    if(contN>0.0)smoothY[i]/=contN;
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

  return(smoother);
}/*setSmoother*/


/*##########################################*/
/*mark funtion fit features*/

void writeFitResults(float lmRange,float lmInt,float lmSig,int numb,char max,FILE *opoo,float lmErr,int sum,float qRange,float qInt,int zen,int az,float midR,float zenAng,float azAng,int waveStart)
{
  float x=0,y=0,z=0;

  x=qRange*cos(azAng)*sin(zenAng);
  y=qRange*sin(azAng)*sin(zenAng);
  z=qRange*cos(zenAng);

  fprintf(opoo,"%d %d %d %d %f %f %f %f %f %f %f %f %f %f %d\n",zen,az,max,sum,lmRange,lmInt,lmSig,lmErr,qRange,qInt,midR,x,y,z,waveStart);

  return;
}/*writeFitResults*/


/*##########################################*/
/*fit a a Gaussian by Levenberg-Marquardt*/

void fitFunction(char *data,int numb,float *range,float *integral,float *width,float rOffset,int band,int azInd,float *error,char *outRoot,int offset)
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
  int mrqminMine(float *,float *,float *,int,float *,int *,int,float **,float **,float *,void (*)(float,float *,float *,float *,int),float *);
  int indicator=0;
  char found=1;

  /*tolerances*/
  minErr=0.001;
  maxErr=10.0;

  /*set up arrays to pass to mrqmin*/
  temp=filterData(data,numb);
  x=falloc(numb+1,"range",0);
  for(i=1;i<=numb;i++)x[i]=(float)(i-1)*res+rOffset;

  /*number of Gaussians to fit*/
  nGauss=1;
  nParams=3*nGauss;
  params=initialGuess(x,temp,numb,nGauss);

  /*numb is the number of return bins above noise in this feature*/
  if(numb<=3){  /*one point above noise, fix max and width*/
    params[3]=0.11;
  }else{        /*there are enough points to fit to*/
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

    /*for(i=1;i<=numb;i++)fprintf(stdout,"iter %f %f\n",x[i],temp[i]);*/

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
  }/*number of points check*/

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


/*##########################################*/
/*find outgoing pulse and check saturation*/

int findStart(int start,int end,char *satTest,char *data,int offset)
{
  int i=0,b=0,e=0;  /*loop control and bounds*/
  int place=0;      /*array index*/
  int waveStart=0;  /*outgoing pulse bin*/
  char max=0;       /*max intensity*/
  char satThresh=0; /*saturation threshold*/

  waveStart=-1; /*start-sOffset;*/

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
    fprintf(opoo[band],"# 1 zenInd, 2 azInd, 3 max, 4 sum, 5 lmRange, 6 lmInt, 7 lmSig, 8 lmErr, 9 qRange, 10 qInt, 11 midR, 12 x, 13 y, 14 z, 15 waveStart\n");
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


/*the end*/
/*##########################################*/

