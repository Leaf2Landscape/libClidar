#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "mpfit.h"
#include "libLidVoxel.h"
#include "libLasProcess.h"
#include "libOptWaveform.h"


/*#####################*/
/*# Reads TLS points  #*/
/*# within ALS beam   #*/
/*# and optimises TLS #*/
/*# voxel parameters  #*/
/*# S Hancock 05 2015 #*/
/*#####################*/

/*#######################################*/
/*# Copyright 2014-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of the lidar voxel library, compareWaves.          #*/
/*#                                                                      #*/
/*# compareWaves is free software: you can redistribute it and/or modify #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# compareWaves is distributed in the hope that it will be useful,      #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#  You should have received a copy of the GNU General Public License   #*/
/*#  along with compareWaves.  If not, see <http://www.gnu.org/licenses/>#*/
/*########################################################################*/

#define TOL 0.0001


/*##########################################*/
/*control structure*/

typedef struct{
  /*input output*/
  char **inList;       /*input filename list*/
  char listNamen[200]; /*input list filename*/
  char pNamen[400];    /*pulse shape filename*/
  char outNamen[200];  /*output filename*/
  char statsNamen[200];/*stats output file*/

  /*options*/
  char optDenoise;     /*optimise denoising parameters*/
  char firstLast;      /*optimise first-last returns*/
  char readTLSpar;     /*read TLS params from input file*/
  char doDenoise;      /*perform denoising without optimising*/
  char runCheck;       /*write a mark if we're running*/
  char grabDir[200];   /*directory to write mark to*/

  /*denoising parameter initial guesses*/
  denPar iniPar;

  /*tls parameters if we are reading from a file*/
  float *appRefl;
  float *minGap;

  /*optimisation control*/
  optControl ofimage;
}control;


/*##########################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  overStruct *data=NULL;
  overStruct *readData(control *);
  overStruct *tidyOverStruct(overStruct *);
  void setPars(overStruct *);
  void pulseReadForTLS(char *,overStruct *);
  void fitTLSparams(overStruct *);
  void outputTLS(overStruct *,control *);
  void outputDenoise(overStruct *,control *); 
  void makeTLSwaves(overStruct *);
  void optimiseDenoise(overStruct *,control *);
  void optimiseFirstLast(overStruct *,control *);
  void setDenoise(overStruct *,control *);
  void copyTLSpar(overStruct *,control *);
  void readListFile(control *);
  void writeRunCheck(control *);
  void tidySMoothPulse() ;


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*if needed, mark that we're running*/
  if(dimage->runCheck)writeRunCheck(dimage);

  /*read input list*/
  readListFile(dimage);

  /*read TLS and ALS data*/
  data=readData(dimage);

  /*set default parameters*/
  setPars(data);

  /*read pulse and rearrange pulse for convolution*/
  pulseReadForTLS(dimage->pNamen,data);

  /*optimise TLS parameters*/
  if(!dimage->readTLSpar)fitTLSparams(data);
  else                   copyTLSpar(data,dimage);

  /*optimise ALS denoising*/
  if(dimage->optDenoise||dimage->doDenoise){
    /*produce TLS waveforms*/
    makeTLSwaves(data);

    /*set up denoising parameters*/
    setDenoise(data,dimage);

    /*optimise denoising*/
    if(dimage->optDenoise){
      /*optimise ALS denoising*/
      if(!dimage->firstLast)optimiseDenoise(data,dimage);
      else                  optimiseFirstLast(data,dimage);
    }/*optimise denoising*/

    /*output results*/
    outputDenoise(data,dimage);
  }else{  /*only fitted TLS*/
    /*print out fitted TLS wave*/
    outputTLS(data,dimage);
  }

  /*tidy up arrays*/
  data=tidyOverStruct(data);
  tidySMoothPulse();
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->ofimage.nFiles);
    TIDY(dimage->appRefl);
    TIDY(dimage->minGap);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##########################################*/
/*optimise first last returns*/

void optimiseFirstLast(overStruct *data,control *dimage)
{
  int nParams=0,totNumb=0;
  int firstLastError(int,int,double *,double *,double **,void *);
  double *doPars=NULL;
  mp_par *mparS=NULL;
  mp_par *denBounds(int,overStruct *);
  mp_result *result=NULL;
  mp_config *config=NULL;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int fitCheck=0;

  fprintf(stdout,"Optimising first last\n");
  /*set initial parameter guess*/
  nParams=11;
  doPars=dalloc(nParams,"pars",0);
  doPars[0]=data->den.thresh;      /*noise threshold*/
  doPars[1]=data->den.meanN;       /*mean noise level*/
  doPars[2]=data->den.tailThresh;  /*trailing edge threshold, meanN by default*/
  doPars[3]=data->den.sWidth;      /*smoothing length*/
  doPars[4]=data->den.psWidth;     /*pre-smoothing width*/
  doPars[5]=data->den.minWidth;    /*min width above noise to accept*/
  doPars[6]=data->den.medLen;      /*median filter width*/
  doPars[7]=data->den.pScale;      /*scale pulse length by*/
  doPars[8]=log(data->den.deChang)/log(10.0);  /*deChang; 10^    change between decon iterations to stop*/
  doPars[9]=data->den.threshScale; /*variable noise threshold scale*/
  doPars[9]=data->den.gWidth;      /*Gaussian fit smoothing*/


  /*two points per beam, first and last*/
  totNumb=dimage->ofimage.nFiles*2;

  /*option space*/
  if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  config->maxiter=4000;

  /*result space*/
  if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
    fprintf(stderr,"error in mpfit structure.\n");
    exit(1);
  }
  result->resid=dalloc(totNumb,"",0);
  result->xerror=dalloc(nParams,"",0);
  result->covar=dalloc(nParams*nParams,"",0);
  /*bounds and more options*/
  mparS=denBounds(nParams,data);

  /*optimise parameters*/
  fitCheck=mpfit(firstLastError,totNumb,nParams,doPars,mparS,config,(void *)data,result);
  if(fitCheck<0){
    fprintf(stderr,"Check %d\n",fitCheck);
    exit(1);
  }

  /*copy parameters back*/
  data->den.thresh=(float)doPars[0];      /*noise threshold*/
  data->den.meanN=(float)doPars[1];       /*mean noise level*/
  data->den.tailThresh=(float)doPars[2];  /*trailing edge threshold, meanN by default*/
  data->den.sWidth=(float)doPars[3];      /*smoothing length*/
  data->den.psWidth=(float)doPars[4];     /*pre-smoothing width*/
  data->den.minWidth=(int)doPars[5];      /*min width above noise to accept*/
  data->den.medLen=(int)doPars[6];        /*median filter width*/
  data->den.pScale=(float)doPars[7];      /*scale pulse length by*/
  data->den.deChang=pow(10.0,doPars[8]);  /*change between decon iterations to stop*/
  data->den.threshScale=(float)doPars[9]; /*variable threshold scale*/
  data->den.gWidth=(float)doPars[10];     /*Gaussian fit feature smoothing*/


  /*tidy arrays*/
  TIDY(doPars);
  TIDY(mparS);
  if(result){
    TIDY(result->resid);
    TIDY(result->xerror);
    TIDY(result->covar);
    TIDY(result);
  }
  TIDY(config);
  return;
}/*optimiseFirstLast*/


/*##########################################*/
/*optimise ALS denoising*/

void optimiseDenoise(overStruct *data,control *dimage)
{
  int numb=0;
  int nParams=0,totNumb=0;
  int denError(int,int,double *,double *,double **,void *);
  double *doPars=NULL;
  mp_par *mparS=NULL;
  mp_par *denBounds(int,overStruct *);
  mp_result *result=NULL;
  mp_config *config=NULL;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int fitCheck=0;

  fprintf(stdout,"Optimising denoising\n");

  /*set initial parameter guess*/
  nParams=11;
  doPars=dalloc(nParams,"pars",0);
  doPars[0]=data->den.thresh;      /*noise threshold*/
  doPars[1]=data->den.meanN;       /*mean noise level*/
  doPars[2]=data->den.tailThresh;  /*trailing edge threshold, meanN by default*/
  doPars[3]=data->den.sWidth;      /*smoothing length*/
  doPars[4]=data->den.psWidth;     /*pre-smoothing width*/
  doPars[5]=data->den.minWidth;    /*min width above noise to accept*/
  doPars[6]=data->den.medLen;      /*median filter width*/
  doPars[7]=data->den.pScale;      /*scale pulse length by*/
  doPars[8]=log(data->den.deChang)/log(10.0);  /*deChang; 10^    change between decon iterations to stop*/
  doPars[9]=data->den.threshScale; /*variable noise threshold scale*/
  doPars[10]=data->den.gWidth;      /*Gaussian fit smoothing*/


  /*count up the total number of data points*/
  totNumb=0;
  for(numb=0;numb<dimage->ofimage.nFiles;numb++)totNumb+=data->dArr[numb].nALS;

  /*option space*/
  if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  config->maxiter=1000;

  /*result space*/
  if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
    fprintf(stderr,"error in mpfit structure.\n");
    exit(1);
  }
  result->resid=dalloc(totNumb,"",0);
  result->xerror=dalloc(nParams,"",0);
  result->covar=dalloc(nParams*nParams,"",0);
  /*bounds and more options*/
  mparS=denBounds(nParams,data);

  /*optimise parameters*/
  fitCheck=mpfit(denError,totNumb,nParams,doPars,mparS,config,(void *)data,result);
  if(fitCheck<0){
    fprintf(stderr,"Check %d\n",fitCheck);
    exit(1);
  }

  /*copy parameters back*/
  data->den.thresh=(float)doPars[0];      /*noise threshold*/
  data->den.meanN=(float)doPars[1];       /*mean noise level*/
  data->den.tailThresh=(float)doPars[2];  /*trailing edge threshold, meanN by default*/
  data->den.sWidth=(float)doPars[3];      /*smoothing length*/
  data->den.psWidth=(float)doPars[4];     /*pre-smoothing width*/
  data->den.minWidth=(int)doPars[5];      /*min width above noise to accept*/
  data->den.medLen=(int)doPars[6];        /*median filter width*/
  data->den.pScale=(float)doPars[7];      /*scale pulse length by*/
  data->den.deChang=pow(10.0,doPars[8]);  /*change between decon iterations to stop*/
  data->den.threshScale=(float)doPars[9]; /*variable threshold scale*/
  data->den.gWidth=(float)doPars[10];     /*Gaussian fit feature smoothing*/


  /*tidy arrays*/
  TIDY(doPars);
  TIDY(mparS);
  if(result){
    TIDY(result->resid);
    TIDY(result->xerror);
    TIDY(result->covar);
    TIDY(result);
  }
  TIDY(config);
  return;
}/*optimiseDenoise*/


/*##########################################*/
/*output ALS denoising results*/

void outputDenoise(overStruct *data,control *dimage)
{
  int i=0,numb=0;
  float *denoised=NULL;
  float aE=0,tE=0,aEraw=0,tEraw=0;
  float tlsRMSE=0,denRMSE=0;
  float *smoothDen=NULL;
  float *alsArea=NULL,*alsUse=NULL;
  float *correctAttenuation(float *,int);
  float err=0,weight=0,subFrac=0;
  float meanDenRMSE=0,meanTlsRMSE=0;
  char subT=0,failed=0;;
  char namen[300];
  void readPulse(denPar *);
  FILE *opoo=NULL;

  if(data->den.pulse==NULL)readPulse(&data->den);

  /*loop over files*/
  subT=failed=0;
  meanDenRMSE=meanTlsRMSE=0.0;
  for(numb=0;numb<dimage->ofimage.nFiles;numb++){
    if(dimage->ofimage.nFiles==1)strcpy(namen,dimage->outNamen);
    else                 sprintf(namen,"%s.%d.fitted",dimage->outNamen,numb);
    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }

    /*output waveforms*/
    denoised=processWave(data->dArr[numb].rawALS,data->dArr[numb].nALS,&data->den,data->dArr[numb].gbic);
    alsArea=correctAttenuation(denoised,data->dArr[numb].nALS);

    TIDY(data->den.gPar);

    /*set pointer to relevant array for error*/
    if(dimage->ofimage.testAtten==0)alsUse=denoised;
    else                            alsUse=alsArea;

    /*Smooth ALS area*/
    smoothDen=smooth(dimage->ofimage.denWidth,data->dArr[numb].nALS,alsUse,dimage->ofimage.res);

    /*energy*/
    aE=tE=aEraw=tEraw=0.0;
    for(i=0;i<data->dArr[numb].nALS;i++){
      aE+=smoothDen[i];
      tE+=data->dArr[numb].smoothTLS[i];
      aEraw+=data->dArr[numb].als[i];
      tEraw+=data->dArr[numb].tlsWave[0][i];
    }
    if(aE<=0.0)aE=1.0;
    if(isnan(aE))failed=1;

    /*errors*/
    if(dimage->firstLast==0){  /*wave deviation*/
      tlsRMSE=denRMSE=0.0;
      for(i=0;i<data->dArr[numb].nALS;i++){
        tlsRMSE+=(data->dArr[numb].als[i]/aEraw-data->dArr[numb].tlsWave[0][i]/tEraw)*\
                 (data->dArr[numb].als[i]/aEraw-data->dArr[numb].tlsWave[0][i]/tEraw);
        err=smoothDen[i]/aE-data->dArr[numb].smoothTLS[i]/tE;
        if(dimage->ofimage.weight){
          if((i<data->dArr[numb].tStart)&&(err>0.0)){
            weight=(float)(data->dArr[numb].tStart-i)*2.6;
          }else if((i>data->dArr[numb].tEnd)&&(err>0.0)){ /*weight to avoid false positives outside envelope*/
            weight=(float)(i-data->dArr[numb].tEnd)*2.6;
            subFrac+=smoothDen[i]/aE;
            subT=1;   
          }else  weight=1.0;
        }else weight=1.0;
        if(dimage->firstLast==0)denRMSE+=weight*err*err;
      }
      if(data->dArr[numb].nALS>0)tlsRMSE=sqrt(tlsRMSE/(float)data->dArr[numb].nALS);
      if(data->dArr[numb].nALS>0)denRMSE=sqrt(denRMSE/(float)data->dArr[numb].nALS);
    }else{   /*first last deviation*/
      tlsRMSE=denRMSE=0.0;
      /*first return*/
      for(i=0;i<data->dArr[numb].nALS;i++){
        if(smoothDen[i]>TOL)break;
      }
      if(dimage->ofimage.weight){
        if(i<(data->dArr[numb].tStart-1))weight=dimage->ofimage.flWeight;
        else                             weight=1.0;
      }else weight=1.0;
      denRMSE+=weight*(float)(data->dArr[numb].tStart-i)*(float)(data->dArr[numb].tStart-i);

      /*last return*/
      for(i=data->dArr[numb].nALS-1;i>=0;i--){
        if(smoothDen[i]>TOL)break;
      }
      if(dimage->ofimage.weight){
        if(i>(data->dArr[numb].tEnd+1))weight=dimage->ofimage.flWeight;
        else                           weight=1.0;
      }else weight=1.0;
      denRMSE+=weight*(float)(data->dArr[numb].tEnd-i)*(float)(data->dArr[numb].tEnd-i);
      denRMSE/=2.0;
    }
    meanTlsRMSE+=tlsRMSE;
    meanDenRMSE+=denRMSE;
    TIDY(smoothDen);
    alsUse=NULL;

    fprintf(opoo,"# 1 bin, 2 TLS-als, 3 TLS visible, 4 TLS area, 5 raw ALS, 6 denoised ALS, 7 deconvolved ALS, 8 alsArea\n");
    for(i=0;i<data->dArr[numb].nALS;i++)fprintf(opoo,"%d %f %f %f %d %f %f %f\n",i,data->dArr[numb].tlsWave[0][i]*aE/data->dArr[numb].tE,\
                                         data->dArr[numb].tlsWave[1][i]*aE/data->dArr[numb].tE,
                                         data->dArr[numb].tlsWave[2][i]*aE/data->dArr[numb].tE,\
                                         data->dArr[numb].rawALS[i],data->dArr[numb].als[i],denoised[i],alsArea[i]);

    /*output parameters*/
    fprintf(opoo,"# denPar %f %f %f %f %f %d %d %f %0.14f %f %f\n",data->den.thresh,data->den.meanN,data->den.tailThresh,\
                      data->den.sWidth,data->den.psWidth,data->den.minWidth,data->den.medLen,data->den.pScale,data->den.deChang,data->den.threshScale,data->den.gWidth); 
    fprintf(opoo,"# appRefl %f minGap %f\n",data->dArr[numb].tlsPar.appRefl,data->dArr[numb].tlsPar.minGap);
    fprintf(opoo,"# errors TLS %f den %f energy %f\n",tlsRMSE,denRMSE,aE);

    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",namen);

    TIDY(denoised);
    TIDY(alsArea);
  }/*file loop*/

  /*tag stats on end of file*/
  meanDenRMSE/=(float)dimage->ofimage.nFiles;
  meanTlsRMSE/=(float)dimage->ofimage.nFiles;
  subFrac/=(float)dimage->ofimage.nFiles;
  if((opoo=fopen(dimage->statsNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->statsNamen);
    exit(1);
  }
  fprintf(opoo,"%f %f %d %f %d %f %f %f %f %f %d %d %f %0.14f %f %f %s\n",meanDenRMSE,meanTlsRMSE,subT,subFrac,failed,\
          data->den.thresh,data->den.meanN,data->den.tailThresh,data->den.sWidth,data->den.psWidth,\
          data->den.minWidth,data->den.medLen,data->den.pScale,data->den.deChang,data->den.threshScale,\
          data->den.gWidth,dimage->outNamen);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }

  return;
}/*outputDenoise*/


/*##########################################*/
/*output TLS results and waveform*/

void outputTLS(overStruct *data,control *dimage)
{
  int numb=0;
  int i=0,contN=0;
  double *doPar=NULL;
  float *tlsWave=NULL;
  float *makeOneTLSwaveform(overStruct *,double *,int);
  float tE=0;
  float RMSE=0;
  char namen[200];
  FILE *opoo=NULL;


  for(numb=0;numb<dimage->ofimage.nFiles;numb++){
    if(dimage->ofimage.nFiles==1)strcpy(namen,dimage->outNamen);
    else                 sprintf(namen,"%s.%d.fitted",dimage->outNamen,numb);

    if((opoo=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }

    doPar=dalloc(2,"",0);
    doPar[0]=(double)data->dArr[numb].tlsPar.appRefl;
    doPar[1]=(double)data->dArr[numb].tlsPar.minGap;
    tlsWave=makeOneTLSwaveform(data,doPar,numb);
    TIDY(doPar);

    /*match intensities*/
    tE=0.0;
    for(i=0;i<data->dArr[numb].nALS;i++)tE+=tlsWave[i];

    RMSE=0.0;
    for(i=0;i<data->dArr[numb].nALS;i++){
      fprintf(opoo,"%d %f %f\n",i,tlsWave[i]*data->dArr[numb].aE/tE,data->dArr[numb].als[i]);
      if(i<150){
        RMSE+=(tlsWave[i]*data->dArr[numb].aE/tE-data->dArr[numb].als[i])*(tlsWave[i]*data->dArr[numb].aE/tE-data->dArr[numb].als[i]);
        contN++;
      }
    }
    if(contN>0)RMSE=sqrt(RMSE/(float)contN);
    else       RMSE=-1.0;
    fprintf(opoo,"# appRefl %f minGap %f RMSE %f\n",data->dArr[numb].tlsPar.appRefl,data->dArr[numb].tlsPar.minGap,RMSE);
    fprintf(opoo,"# coord %.2f %.2f %.2f\n",data->dArr[numb].x0,data->dArr[numb].y0,data->dArr[numb].z0);

    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",namen);
  }/*file loop*/

  return;
}/*outputTLS*/


/*################################################*/
/*copy TLS parameters*/

void copyTLSpar(overStruct *data,control *dimage)
{
  int i=0;

  for(i=0;i<dimage->ofimage.nFiles;i++){
    data->dArr[i].tlsPar.appRefl=dimage->appRefl[i];
    data->dArr[i].tlsPar.minGap=dimage->minGap[i];
  }

  TIDY(dimage->minGap);
  TIDY(dimage->appRefl);
  return;
}/*copyTLSpar*/


/*##########################################*/
/*set bounds for denoising optimisation*/

mp_par *denBounds(int nParams,overStruct *data)
{
  mp_par *mparS=NULL;

  if(!(mparS=(mp_par *)calloc(nParams,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

  if(data->den.varNoise==0){  /*fixed noise thresholds*/
    /*thresh;      noise threshold*/
    mparS[0].fixed=0;
    mparS[0].side=2;
    mparS[0].limited[0]=1;
    mparS[0].limits[0]=12.0;
    mparS[0].limited[1]=1;
    mparS[0].limits[1]=24.0;
    mparS[0].step=4.0;
    /*meanN;       mean noise level*/
    mparS[1].fixed=1;
    /*tailThresh;  trailing edge threshold, meanN by default*/
    if(data->den.tailThresh>=0.0){
      mparS[2].fixed=0;
      mparS[2].side=2;
      mparS[2].limited[0]=1;
      mparS[2].limits[0]=12.0;
      mparS[2].limited[1]=1;
      mparS[2].limits[1]=30.0;
      mparS[2].step=4.0;
    }else mparS[2].fixed=1;
    /*threshScale;   variable noise thresold scale from mean*/
    mparS[9].fixed=1;
  }else{               /*variable noise threshold*/
    mparS[0].fixed=1; /*thresh;      noise threshold*/
    mparS[1].fixed=1; /*meanN;       mean noise level*/
    /*tailThresh;  trailing edge threshold offset from leading threshold*/
    if(data->den.tailThresh>=0.0){
      mparS[2].fixed=0;
      mparS[2].side=2;
      mparS[2].limited[0]=1;
      mparS[2].limits[0]=0.0;
      mparS[2].limited[1]=1;
      mparS[2].limits[1]=10.0;
      mparS[2].step=2.0;
    }else mparS[2].fixed=1;
    /*threshScale;   variable noise thresold scale from mean*/
    mparS[9].fixed=0;
    mparS[9].side=2;
    mparS[9].limited[0]=1;
    mparS[9].limits[0]=0.2;
    mparS[9].limited[1]=1;
    mparS[9].limits[1]=5.0;
    mparS[9].step=0.1;
  }/*noise threshold switch*/

  /*sWidth;      smoothing length*/
  if(data->dimage->optSmoo){
    mparS[3].fixed=0;
    mparS[3].side=2;
    mparS[3].limited[0]=1;
    mparS[3].limits[0]=0.0;
    mparS[3].limited[1]=1;
    mparS[3].limits[1]=1.5;
    mparS[3].step=0.5;
  }else mparS[3].fixed=1;
  /*psWidth;     pre-smoothing width*/
  if(data->dimage->optPsmoo){
    mparS[4].fixed=0;
    mparS[4].side=2;
    mparS[4].limited[0]=1;
    mparS[4].limits[0]=0.0;
    mparS[4].limited[1]=1;
    mparS[4].limits[1]=5.0;
    mparS[4].step=0.5;
  }else mparS[4].fixed=1;
  /*minWidth;    min width above noise to accept*/
  mparS[5].fixed=0;
  mparS[5].side=2;
  mparS[5].limited[0]=1;
  mparS[5].limits[0]=1.0;
  mparS[5].limited[1]=1;
  mparS[5].limits[1]=20.0;
  mparS[5].step=2.0;
  /*medLen;      median filter width*/
  mparS[6].fixed=1;
  //mparS[6].side=2;
  //mparS[6].limited[0]=1;
  //mparS[6].limits[0]=0.0;
  //mparS[6].limited[1]=1;
  //mparS[6].limits[1]=5.0;
  //mparS[6].step=0.5;


  /*deconvolution*/
  if(data->den.deconMeth>=0){
    /*pScale;      scale pulse length by*/
    if(data->dimage->fixPS==0){
      mparS[7].fixed=0;
      mparS[7].side=2;
      mparS[7].limited[0]=1;
      mparS[7].limits[0]=0.75;
      mparS[7].limited[1]=1;
      mparS[7].limits[1]=3.0;
      mparS[7].step=0.1;
    }else mparS[7].fixed=1;  /*fix, ideally to one*/
    /*deChang;     change between decon iterations to stop*/
    mparS[8].fixed=0;
    mparS[8].side=2;
    mparS[8].limited[0]=1;
    mparS[8].limits[0]=-14.0; //0.0000000001;
    mparS[8].limited[1]=1;
    mparS[8].limits[1]=-4.0; //0.0001;
    mparS[8].step=1.0;
  }else{  /*no deconvolution*/
    mparS[7].fixed=1; /*pScale;      scale pulse length by*/
    mparS[8].fixed=1; /*deChang;     change between decon iterations to stop*/
  }/*deconvolution switch*/

  /*Gaussian fitting*/
  if(data->den.fitGauss){
  /*gWidth;  Gaussian fitting*/
    mparS[10].fixed=0;
    mparS[10].side=2;
    mparS[10].limited[0]=1;
    mparS[10].limits[0]=0.0;
    mparS[10].limited[1]=1;
    mparS[10].limits[1]=5.0;
    mparS[10].step=0.2;
  }else{
  /*gWidth;  Gaussian fitting*/
    mparS[10].fixed=1;
  }

  return(mparS);
}/*denBounds*/


/*##########################################*/
/*set fixed denoising parameters*/

void setDenoise(overStruct *data,control *dimage)
{

  /*fixed parameters*/
  data->den.maxIter=2000;
  data->den.res=dimage->ofimage.res;

  /*pulse filename*/
  strcpy(data->den.pNamen,dimage->pNamen);

  /*initial guess for denoising*/
  data->den.thresh=dimage->iniPar.thresh;          /*noise threshold*/
  data->den.meanN=dimage->iniPar.meanN;            /*mean noise level*/
  data->den.tailThresh=dimage->iniPar.tailThresh;  /*trailing edge threshold, meanN by default*/
  data->den.sWidth=dimage->iniPar.sWidth;          /*smoothing length*/
  data->den.psWidth=dimage->iniPar.psWidth;        /*pre-smoothing width*/
  data->den.minWidth=dimage->iniPar.minWidth;      /*min width above noise to accept*/
  data->den.medLen=dimage->iniPar.medLen;          /*median filter width*/
  data->den.pScale=dimage->iniPar.pScale;          /*scale pulse length by*/
  data->den.deChang=dimage->iniPar.deChang;        /*change between decon iterations to stop*/
  data->den.varNoise=dimage->iniPar.varNoise;      /*variable noise switch*/
  data->den.medStats=dimage->iniPar.medStats;      /*use median rather than mean statistics (ARSF)*/
  data->den.threshScale=dimage->iniPar.threshScale;/*variable noise threshold scale*/
  data->den.gWidth=dimage->iniPar.gWidth;          /*Gaussian fit smoothing*/
  data->den.deconMeth=dimage->iniPar.deconMeth;    /*Gold or R-L method*/
  data->den.fitGauss=dimage->iniPar.fitGauss;      /*Gaussian fit switch*/
  data->den.gaussPulse=dimage->iniPar.gaussPulse;  /*pulse to Gaussian switch*/
  data->den.noiseTrack=dimage->iniPar.noiseTrack;  /*noise tracking switch*/
  data->den.matchHard=dimage->iniPar.matchHard;    /*match hard target switch*/
  data->den.hardThresh=dimage->iniPar.hardThresh;  /*hard target threshold*/
  data->den.bitRate=dimage->iniPar.bitRate;        /*bit rate for noise*/
  data->den.maxDN=dimage->iniPar.maxDN;            /*maximum DN for scaling*/

  return;
}/*setDenoise*/


/*##########################################*/
/*ALS first-last error function*/

int firstLastError(int numb,int npar,double *doPars,double *deviates,double **derivs,void *private)
{
  int i=0,fNumb=0;
  int globI=0;    /*global counter for multiple files*/
  float *denoised=NULL;
  float *smoothALS=NULL;
  float weight=0;
  overStruct *data=NULL;
  void readPulse(denPar *);
  float *correctAttenuation(float *,int);
  float *trueArea=NULL;

  /*set data structure from void pointer*/
  data=(overStruct *)private;

  /*copy parameters to denoising structure*/
  data->den.thresh=(float)doPars[0];      /*noise threshold*/
  data->den.meanN=(float)doPars[1];       /*mean noise level*/
  data->den.tailThresh=(float)doPars[2];  /*trailing edge threshold, meanN by default*/
  data->den.sWidth=(float)doPars[3];      /*smoothing length*/
  data->den.psWidth=(float)doPars[4];     /*pre-smoothing width*/
  data->den.minWidth=(int)doPars[5];      /*min width above noise to accept*/
  data->den.medLen=(int)doPars[6];        /*median filter width*/
  data->den.pScale=(float)doPars[7];      /*scale pulse length by*/
  data->den.deChang=pow(10.0,doPars[8]);  /*change between decon iterations to stop*/
  data->den.threshScale=(float)doPars[9]; /*variable threshold scale*/
  data->den.gWidth=(float)doPars[10];     /*Gaussian fit smoothing*/


  if(data->den.varNoise==0){
    /*if thresh<meanN noise tracking causes a blank denoised and a seg fault*/
    if(data->den.thresh<data->den.meanN)data->den.thresh=data->den.meanN;
    if(data->den.tailThresh<data->den.meanN)data->den.tailThresh=data->den.meanN;
  }


  /*read pulse here as it gets smoothed and scaled*/
  TTIDY((void **)data->den.pulse,2);
  data->den.pulse=NULL;
  readPulse(&data->den);

  /*loop over files*/
  globI=0;
  for(fNumb=0;fNumb<data->dimage->nFiles;fNumb++){
    /*denoise ALS*/
    denoised=processWave(data->dArr[fNumb].rawALS,data->dArr[fNumb].nALS,&data->den,data->dArr[fNumb].gbic);

    /*attenuation corrected waveforms*/
    trueArea=correctAttenuation(denoised,data->dArr[fNumb].nALS);
    TIDY(denoised);

    /*smooth both a little to account for small range differences*/
    smoothALS=smooth(data->dimage->denWidth,data->dArr[fNumb].nALS,trueArea,data->den.res);
    TIDY(denoised); 
    if(globI>numb){
      fprintf(stderr,"Point number issue %d of %d\n",globI,numb);
      exit(1);
    }


    /*first return*/
    for(i=0;i<data->dArr[fNumb].nALS;i++){
      if(smoothALS[i]>TOL)break;
    }
    if(data->dimage->weight){
      if(i<(data->dArr[fNumb].tStart-1))weight=data->dimage->flWeight;
      else                              weight=1.0;
    }else weight=1.0;
    deviates[globI++]=(double)(i-data->dArr[fNumb].tStart)*(double)weight;

    /*last return*/
    for(i=data->dArr[fNumb].nALS-1;i>=0;i--){
      if(smoothALS[i]>TOL)break;
    }
    if(data->dimage->weight){
      if(i>(data->dArr[numb].tEnd+1))weight=data->dimage->flWeight;
      else                           weight=1.0;
    }else weight=1.0;
    deviates[globI++]=(double)(data->dArr[numb].tEnd-i)*(double)weight;

    TIDY(data->den.gPar);
    TIDY(smoothALS);
  }/*file loop*/


  data=NULL;
  return(0);
}/*firstLastError*/


/*##########################################*/
/*ALS waveform error function*/

int denError(int numb,int npar,double *doPars,double *deviates,double **derivs,void *private)
{
  int i=0,fNumb=0;
  int globI=0;    /*global counter for multiple files*/
  float aE=0;     /*waveform energies*/
  float *denoised=NULL;
  float *smoothALS=NULL;
  float err=0,weight=0;
  overStruct *data=NULL;
  void readPulse(denPar *);
  float *correctAttenuation(float *,int);
  float *trueArea=NULL;


  /*set data structure from void pointer*/
  data=(overStruct *)private;

  /*copy parameters to denoising structure*/
  data->den.thresh=(float)doPars[0];      /*noise threshold*/
  data->den.meanN=(float)doPars[1];       /*mean noise level*/
  data->den.tailThresh=(float)doPars[2];  /*trailing edge threshold, meanN by default*/
  data->den.sWidth=(float)doPars[3];      /*smoothing length*/
  data->den.psWidth=(float)doPars[4];     /*pre-smoothing width*/
  data->den.minWidth=(int)doPars[5];      /*min width above noise to accept*/
  data->den.medLen=(int)doPars[6];        /*median filter width*/
  data->den.pScale=(float)doPars[7];      /*scale pulse length by*/
  data->den.deChang=pow(10.0,doPars[8]);  /*change between decon iterations to stop*/
  data->den.threshScale=(float)doPars[9]; /*variable threshold scale*/
  data->den.gWidth=(float)doPars[10];     /*Gaussian fit smoothing*/


  if(data->den.varNoise==0){
    /*if thresh<meanN noise tracking causes a blank denoised and a seg fault*/
    if(data->den.thresh<data->den.meanN)data->den.thresh=data->den.meanN;
    if(data->den.tailThresh<data->den.meanN)data->den.tailThresh=data->den.meanN;
  }

  if(data->den.pScale<=0.0){
    fprintf(stderr,"pScale error %f\n",data->den.pScale);
    exit(1);
  }

  /*fprintf(stdout,"t %d mN %d tail %d sW %f psW %f minW %d med %d pScale %f deCh %.10f\n",(int)doPars[0],(int)doPars[1],(int)doPars[2],doPars[3],doPars[4],(int)doPars[5],(int)doPars[6],doPars[7],doPars[8]);*/

  /*read pulse here as it gets smoothed and scaled*/
  TTIDY((void **)data->den.pulse,2);
  data->den.pulse=NULL;
  readPulse(&data->den);   

  /*loop over files*/
  globI=0;
  for(fNumb=0;fNumb<data->dimage->nFiles;fNumb++){
    /*denoise ALS*/
    denoised=processWave(data->dArr[fNumb].rawALS,data->dArr[fNumb].nALS,&data->den,data->dArr[fNumb].gbic);

    /*if required test attenuation corrected waveforms*/
    if(data->dimage->testAtten==1){
      trueArea=correctAttenuation(denoised,data->dArr[fNumb].nALS);
      TIDY(denoised);
      denoised=trueArea;
      trueArea=NULL;
    }/*attenuation correction*/

    /*smooth both a little to account for small range differences*/
    smoothALS=smooth(data->dimage->denWidth,data->dArr[fNumb].nALS,denoised,data->den.res);
    TIDY(denoised);

    /*calculate deviates*/
    aE=0.0;
    for(i=0;i<data->dArr[fNumb].nALS;i++)aE+=smoothALS[i];   /*integral to balance energy*/
    if(isnan(aE)){  /*check for nan*/
     for(i=0;i<data->dArr[fNumb].nALS;i++)deviates[globI++]*=1.5;  /*add some error to get away*/
    }else{
      for(i=0;i<data->dArr[fNumb].nALS;i++){
        err=smoothALS[i]/aE-data->dArr[fNumb].smoothTLS[i]/data->dArr[fNumb].tE;
        if(data->dimage->weight){
          if((i<data->dArr[fNumb].tStart)&&(err>0.0)){
            weight=(float)(data->dArr[fNumb].tStart-i)*2.6;
          }else if((i>data->dArr[fNumb].tEnd)&&(err>0.0)){ /*weight to avoid false positives outside envelope*/
            weight=(float)(i-data->dArr[fNumb].tEnd)*2.6; 
          }else  weight=1.0; 
        }else weight=1.0;
        deviates[globI++]=(double)(weight*err);
      }
    }
    if(globI>numb){
      fprintf(stderr,"Point number issue %d of %d\n",globI,numb);
      exit(1);
    }

    TIDY(data->den.gPar);
    TIDY(smoothALS);
  }/*file loop*/
  TTIDY((void **)data->den.pulse,2);
  TIDY(data->den.matchPulse);
  data->den.pulse=NULL;
  return(0);
}/*denError*/


/*##########################################*/
/*read data*/

overStruct *readData(control *dimage)
{
  int i=0,bin=0,numb=0;
  overStruct *data=NULL;
  char line[300],temp1[100],temp2[100],temp3[100];
  char temp4[100],temp5[100],temp6[100],temp7[100];
  FILE *ipoo=NULL;


  if(!(data=(overStruct *)calloc(1,sizeof(overStruct)))){
    fprintf(stderr,"error in data allocation.\n");
    exit(1);
  }
  if(!(data->dArr=(dataStruct *)calloc(dimage->ofimage.nFiles,sizeof(dataStruct)))){
    fprintf(stderr,"error in data allocation.\n");
    exit(1);
  }
  data->dimage=&dimage->ofimage;   /*share control pointer with data*/

  /*loop over multiple input files*/
  for(numb=0;numb<dimage->ofimage.nFiles;numb++){
    if((ipoo=fopen(dimage->inList[numb],"r"))==NULL){
      fprintf(stderr,"Error opening input file \"%s\"\n",dimage->inList[numb]);
      exit(1);
    }

    /*count up arrays*/
    data->dArr[numb].nALS=data->dArr[numb].nTLS=0;
    while(fgets(line,300,ipoo)!=NULL){
      if(strncasecmp(line,"#",1)){
        if(strncasecmp(line,"ALS",3))data->dArr[numb].nTLS++;
        else                         data->dArr[numb].nALS++;
      }
    }/*data reading*/

    /*allocate space*/
    data->dArr[numb].als=falloc(data->dArr[numb].nALS,"denoised als data",0);
    data->dArr[numb].rawALS=uchalloc(data->dArr[numb].nALS,"raw als data",0);
    data->dArr[numb].aE=0.0;
    if(!(data->dArr[numb].tls=(tlsPoint *)calloc(data->dArr[numb].nTLS,sizeof(tlsPoint)))){
      fprintf(stderr,"error in tls data allocation.\n");
      exit(1);
    }

    /*rewind*/
    if(fseek(ipoo,0,SEEK_SET)){
      fprintf(stderr,"fseek error\n");
      exit(1);
    }

    /*read data*/
    i=0;
    while(fgets(line,300,ipoo)!=NULL){
      if(!strncasecmp(line,"#",1)){  /*read centre coordinate*/
        if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
          data->dArr[numb].x0=atof(temp2);
          data->dArr[numb].y0=atof(temp3);
          data->dArr[numb].z0=atof(temp4);
        }
      }else{   /*read data*/
        if(!strncasecmp(line,"ALS",3)){  /*read ALS*/
          if(sscanf(line,"%s %s %s %s %s",temp1,temp2,temp3,temp4,temp5)==5){
            bin=atoi(temp2);
            if((bin<0)||(bin>data->dArr[numb].nALS)){
              fprintf(stderr,"Error\n");
              exit(1);
            }
            data->dArr[numb].rawALS[bin]=(unsigned char)atoi(temp3);  /*raw*/
            data->dArr[numb].als[bin]=atof(temp4);                    /*denoised*/
            data->dArr[numb].aE+=data->dArr[numb].als[bin];           /*total energy for normalisation*/
            if(dimage->ofimage.doGBIC)data->dArr[numb].gbic=atof(temp5);                        /*GBIC*/
            else                      data->dArr[numb].gbic=1.0;
          }
        }else{  /*read TLS*/
          if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
            if(i>data->dArr[numb].nTLS){
              fprintf(stderr,"Counting error\n");
              exit(1);
            }
            data->dArr[numb].tls[i].bin=atoi(temp1);
            data->dArr[numb].tls[i].x=atof(temp2);
            data->dArr[numb].tls[i].y=atof(temp3);
            data->dArr[numb].tls[i].z=atof(temp4);
            data->dArr[numb].tls[i].gap=atof(temp5);
            data->dArr[numb].tls[i].r=atof(temp6);
            data->dArr[numb].tls[i].refl=(uint16_t)atoi(temp7);
            i++;
          }
        }
      }
    }/*file reading*/

    if(ipoo){
      fclose(ipoo);
      ipoo=NULL;
    }
  }/*file loop*/

  return(data);
}/*readData*/


/*##########################################*/
/*read a file with list of filenames*/

void readListFile(control *dimage)
{
  int i=0;
  char line[300],temp[300];
  char temp1[200],temp2[200];
  FILE *ipoo=NULL;

  /*if we are only reading a single file then skip this function*/
  if(dimage->ofimage.nFiles==1)return;


  if((ipoo=fopen(dimage->listNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input list file %s\n",dimage->listNamen);
    exit(1);
  }

  /*count number of files*/
  dimage->ofimage.nFiles=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))dimage->ofimage.nFiles++;
  }/*file counting*/

  /*allocate and rewind*/
  dimage->inList=chChalloc(dimage->ofimage.nFiles,"input list",0);
  if(dimage->readTLSpar){
    dimage->appRefl=falloc(dimage->ofimage.nFiles,"appRefl",0);
    dimage->minGap=falloc(dimage->ofimage.nFiles,"minGap",0);
  }
  if(fseek(ipoo,0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  i=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(!dimage->readTLSpar){  /*only read filename*/
        if(sscanf(line,"%s",temp)==1){
          dimage->inList[i]=challoc(strlen(temp)+1,"input list",i+1);
          strcpy(dimage->inList[i],temp);
          i++;
        }
      }else{   /*read filename and tls pars*/
        if(sscanf(line,"%s %s %s",temp,temp1,temp2)==3){
          dimage->inList[i]=challoc(strlen(temp)+1,"input list",i+1);
          strcpy(dimage->inList[i],temp);
          dimage->appRefl[i]=atof(temp1);
          dimage->minGap[i]=atof(temp2);
          i++;
        }
      }
    }
  }/*file counting*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*readListFile*/


/*##########################################*/
/*write a run check file*/

void writeRunCheck(control *dimage)
{
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s/%s.check",dimage->grabDir,dimage->outNamen);

  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening rucn check file %s\n",namen);
    exit(1);
  }

  fprintf(opoo,"Running in this file\n");

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  return;
}/*writeRunCheck*/


/*##########################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);


  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in control allocation.\n");
    exit(1);
  }

  /*defaults*/
  dimage->ofimage.nFiles=1;
  dimage->inList=chChalloc(dimage->ofimage.nFiles,"in list",0);
  dimage->inList[i]=challoc(200,"",0);
  strcpy(dimage->inList[0],"/media/hda/shancock/data/bess/tls_tune/wBits/res.0.5/tlsTune.0.5.103.wBits");

  strcpy(dimage->pNamen,"/Users/stevenhancock/data/bess/leica_shape/leicaPulseFilt.dat");
  strcpy(dimage->outNamen,"resultWave.dat");
  strcpy(dimage->statsNamen,"teastStats.dat");
  dimage->optDenoise=0;
  dimage->ofimage.denWidth=0.0;
  dimage->ofimage.res=0.15;
  dimage->ofimage.testAtten=0;
  dimage->ofimage.doGBIC=0;
  dimage->ofimage.fixGap=0;
  dimage->readTLSpar=0;
  dimage->ofimage.fixPS=0;
  dimage->ofimage.optSmoo=1;
  dimage->ofimage.optPsmoo=0;
  dimage->ofimage.weight=0;
  dimage->firstLast=0;
  dimage->doDenoise=0;
  dimage->ofimage.flWeight=3.6;
  dimage->runCheck=0;


  /*initial denoising guesses*/
  setDenoiseDefault(&dimage->iniPar);
  dimage->iniPar.deconMeth=0;      /*Gold's method*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->ofimage.nFiles);
        dimage->inList=NULL;
        dimage->ofimage.nFiles=1;
        dimage->inList=chChalloc(dimage->ofimage.nFiles,"in list",0);
        dimage->inList[0]=challoc(200,"",0);
        strcpy(dimage->inList[0],argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",6)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->ofimage.nFiles);
        dimage->inList=NULL;
        dimage->ofimage.nFiles=0;  /*set to nonesense value for later*/
        strcpy(dimage->listNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-pFile",6)){
        checkArguments(1,i,argc,"-pFile");
        strcpy(dimage->pNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-statsOut",9)){
        checkArguments(1,i,argc,"-statsOut");
        strcpy(dimage->statsNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-optDenoise",11)){
        dimage->optDenoise=1;
      }else if(!strncasecmp(argv[i],"-denWidth",9)){
        checkArguments(1,i,argc,"-denWidth");
        dimage->ofimage.denWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-atten",6)){
        dimage->ofimage.testAtten=1;
        dimage->optDenoise=1;
      }else if(!strncasecmp(argv[i],"-thresh",7)){
        checkArguments(1,i,argc,"-thresh");
        dimage->iniPar.thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-meanN",6)){
        checkArguments(1,i,argc,"-meanN");
        dimage->iniPar.meanN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-tailThresh",11)){
        checkArguments(1,i,argc,"-tailThresh");
        dimage->iniPar.tailThresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sWidth",7)){
        checkArguments(1,i,argc,"-sWidth");
        dimage->iniPar.sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-psWidth",8)){
        checkArguments(1,i,argc,"-psWidth");
        dimage->iniPar.psWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minWidth",9)){
        checkArguments(1,i,argc,"-minWidth");
        dimage->iniPar.minWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-medLen",7)){
        checkArguments(1,i,argc,"-medLen");
        dimage->iniPar.medLen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pScale",7)){
        checkArguments(1,i,argc,"-pScale");
        dimage->iniPar.pScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-deChang",8)){
        checkArguments(1,i,argc,"-deChang");
        dimage->iniPar.deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-varScale",9)){
        checkArguments(1,i,argc,"-varScale");
        dimage->iniPar.threshScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-varNoise",9)){
        dimage->iniPar.varNoise=1;
        dimage->iniPar.medStats=1;
      }else if(!strncasecmp(argv[i],"-GBIC",5)){
        dimage->ofimage.doGBIC=1;
      }else if(!strncasecmp(argv[i],"-richLucy",9)){
        dimage->iniPar.deconMeth=1;
      }else if(!strncasecmp(argv[i],"-gauss",6)){
        dimage->iniPar.fitGauss=1;
        dimage->iniPar.gaussPulse=1;
      }else if(!strncasecmp(argv[i],"-noDecon",7)){
        dimage->iniPar.deconMeth=-1;
      }else if(!strncasecmp(argv[i],"-gWidth",6)){
        checkArguments(1,i,argc,"-gWidth");
        dimage->iniPar.gWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-readTLSpar",11)){
        dimage->readTLSpar=1;
      }else if(!strncasecmp(argv[i],"-fixPS",6)){
        dimage->ofimage.fixPS=1;
      }else if(!strncasecmp(argv[i],"-noWeight",9)){
        dimage->ofimage.weight=0;
      }else if(!strncasecmp(argv[i],"-fixSmoo",8)){
        dimage->ofimage.optSmoo=0;
      }else if(!strncasecmp(argv[i],"-optPsmoo",9)){
        dimage->ofimage.optPsmoo=1;
      }else if(!strncasecmp(argv[i],"-firstLast",10)){
        dimage->firstLast=1;
        dimage->optDenoise=1;
      }else if(!strncasecmp(argv[i],"-noTrack",8)){
        dimage->iniPar.noiseTrack=0;
      }else if(!strncasecmp(argv[i],"-flWeight",9)){
        checkArguments(1,i,argc,"-flWeight");
        dimage->ofimage.flWeight=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lockGap",8)){
        dimage->ofimage.fixGap=1;
      }else if(!strncasecmp(argv[i],"-findHard",9)){
        dimage->iniPar.matchHard=1;
      }else if(!strncasecmp(argv[i],"-denoise",8)){
        dimage->doDenoise=1;
      }else if(!strncasecmp(argv[i],"-runCheck",9)){
        checkArguments(1,i,argc,"-runCheck");
        dimage->runCheck=1;
        strcpy(dimage->grabDir,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;     input filename\n-inList name;    list of input files to solve simultaneously\n-output name;    output filename\n-pFile name;     pulse filename\n-statsOut name;  stats for multi-fit output\n-optDenoise;     optimise ALS denoising\n-denWidth s;     smoothing width when testing denoising accuracy\n-atten;          optimise with attenuation correction\n-GBIC;           apply GBIC\n-readTLSpar;     read TLS par from file rather than optimise\n-denoise;        perform denoising without optimisation\n-runCheck dir;     write a run check file to a directory\n\nDenoising options\n-varNoise;       use variable noise levels\n-richLucy;       Richardson-Lucy deconvolution rather than Gold's\n\nInitial parameter guess\n-thresh x;\n-meanN x;\n-tailThresh x;\n-sWidth x;\n-psWidth x;\n-minWidth x;\n-medLen x;\n-pScale x;\n-deChang x;\n-varScale x;     variable noise threshold scale\n-gauss;          fit a Gaussian\n-gWidth sig;     smooth to choose Gaussian features\n-fixPS;          fix pulse scaling\n-noWeight;       do not weight fit to avoid subterranean\n-fixSmoo;        fix smoothing width\n-optPsmoo;       optimise pre-smoothing\n-firstLast;      optimise to find first and last return\n-noTrack;        do not use noise tracking\n-flWeight w;     scale factor for first-last weighting\n-lockGap;       set min gap fraction to 1\n-findHard;      find hard targets\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*arg loop*/

  return(dimage);
}/*readCommands*/

/*the end*/
/*##########################################*/

