#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidVoxel.h"
#include "libOptWaveform.h"
#include "mpfit.h"



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

/*####################################*/
/*clean up a data arrays*/

overStruct *tidyOverStruct(overStruct *data)
{
  int i=0;

  if(data){
    if(data->dArr){
      for(i=0;i<data->dimage->nFiles;i++){
        TIDY(data->dArr[i].tls);
        TIDY(data->dArr[i].als);
        TIDY(data->dArr[i].rawALS);
        TTIDY((void **)data->dArr[i].tlsWave,3);
        TIDY(data->dArr[i].smoothTLS);
      }
    }
    if(data->dimage->errOut){
      fclose(data->dimage->errOut);
      data->dimage->errOut=NULL;
    }
    TTIDY((void **)data->den.pulse,2);
    TIDY(data->den.matchPulse);
    TIDY(data->pulse);
    data->dimage=NULL;
    TIDY(data->dArr);
    TIDY(data->den.gPar);
    TIDY(data);
  }

  return(NULL);
}/*tidyOverStruct*/


/*##########################################*/
/*read and rearrange ALS pulse*/

void pulseReadForTLS(char *pNamen,overStruct *data)
{
  int i=0,place=0;
  int pStart=0,pEnd=0;
  int *contN=NULL;
  float max=0,thresh=0;
  float inRes=0;
  denPar tempPars;   /*temporary denoising structure*/
  void readPulse(denPar *);
  void setDenoiseDefault(denPar *);

  /*read pulse*/
  setDenoiseDefault(&tempPars);
  tempPars.pScale=1.0;
  tempPars.sWidth=0.0;
  strcpy(tempPars.pNamen,pNamen);
  tempPars.gaussPulse=0;
  readPulse(&tempPars);

  data->res=0.15;

  /*find pulse start*/
  thresh=0.0002;  /*chosen from a graph of the pulse*/
  pStart=-1;
  for(i=0;i<tempPars.pBins;i++){
    if((pStart<0)&&(tempPars.pulse[1][i]>thresh)){
      pStart=i-1;
      break;
    }
  }/*find maximum and pulse start*/
  
  /*trace back to find end*/
  pEnd=-1;
  for(i=tempPars.pBins-1;i>=0;i--){
    if((pEnd<0)&&(tempPars.pulse[1][i]>thresh)){
      pEnd=i+1;
      break;
    }
  }/*trace back to find end*/

  inRes=tempPars.pulse[0][1]-tempPars.pulse[0][0];
  data->pBins=(int)((tempPars.pulse[0][pEnd]-tempPars.pulse[0][pStart])/data->res)+1;
  data->pulse=falloc(data->pBins,"pulse",0);

  /*tot up pulse to use*/
  contN=ialloc(data->pBins,"count",0);
  for(i=pStart;i<pEnd;i++){
    place=(int)((float)(i-pStart)*inRes/data->res);
    if((place>=0)&&(place<data->pBins)){
      data->pulse[place]+=tempPars.pulse[1][i];
      contN[place]++;
    }
  }

  /*normalise*/
  max=-1000.0;
  for(i=0;i<data->pBins;i++){
    if(contN[i]>0){
      data->pulse[i]/=(float)contN[i];
      if(data->pulse[i]>max){
        max=data->pulse[i];
        data->maxBin=i;
      }
    }else fprintf(stderr,"Whit?\n");
  }/*normalise*/

  TIDY(contN);
  TTIDY((void **)tempPars.pulse,2);
  return;
}/*pulseReadForTLS*/


/*################################################*/
/*fit with minpack*/

void fitTLSparams(overStruct *data)
{
  int numb=0;
  int nParams=0;
  int tlsError(int,int,double *,double *,double **,void *);
  double *doPars=NULL;
  mp_par *mparS=NULL;
  mp_par *tlsBounds(int,char);
  mp_result *result=NULL;
  mp_config *config=NULL;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int fitCheck=0;
  /*variable starting points*/
  float minRefl=0,maxRefl=0,reflRes=0;
  float minGap=0,maxGap=0,gapRes=0;
  float minErr=0;
  int nRefl=0,nGap=0;
  int i=0,j=0;


  minRefl=0.2;
  maxRefl=0.7;
  reflRes=0.1;
  nRefl=(int)((maxRefl-minRefl)/reflRes);
  minGap=0.01;
  maxGap=1.0;
  gapRes=0.2;
  nGap=(int)((maxGap-minGap)/gapRes);
  if(data->dimage->fixGap==1)nGap=1;

  for(numb=0;numb<data->dimage->nFiles;numb++){
    data->thisTLS=numb;    /*TLS index to pass through to optimisation function*/
    fprintf(stdout,"Optimising TLS %d of %d\n",numb+1,data->dimage->nFiles);

    /*set initial parameter guess*/
    nParams=2;
    doPars=dalloc(nParams,"pars",0);
    minErr=0.0;
    for(i=0;i<nRefl;i++){
      doPars[0]=(double)(minRefl+(float)i*reflRes);
      for(j=0;j<nGap;j++){
        if(data->dimage->fixGap==0){
          doPars[1]=(double)(minGap+(float)j*gapRes);
        }else{
          doPars[1]=1.0;
        }

        /*option space*/
        if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
          fprintf(stderr,"error in control structure.\n");
          exit(1);
        }
        //config->nofinitecheck=1;
        //config->stepfactor=10.0;
        config->maxiter=1000;
        //config->douserscale=0.1;

        /*result space*/
        if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
          fprintf(stderr,"error in mpfit structure.\n");
          exit(1);
        }
        result->resid=dalloc(data->dArr[numb].nALS,"",0);
        result->xerror=dalloc(nParams,"",0);
        result->covar=dalloc(nParams*nParams,"",0);
        /*bounds and more options*/
        mparS=tlsBounds(nParams,data->dimage->fixGap);


        /*optimise parameters*/
        fitCheck=mpfit(tlsError,data->dArr[numb].nALS,nParams,doPars,mparS,config,(void *)data,result);
        if(fitCheck<0){
          fprintf(stderr,"Check %d\n",fitCheck);
          exit(1);
        }

        /*copy parameters back*/
        if((float)result->bestnorm<minErr){
          data->dArr[numb].tlsPar.appRefl=(float)doPars[0];
          data->dArr[numb].tlsPar.minGap=(float)doPars[1];
          minErr=(float)result->bestnorm;
        }
        /*tidy arrays*/
        TIDY(mparS);
        if(result){
          TIDY(result->resid);
          TIDY(result->xerror);
          TIDY(result->covar);
          TIDY(result);
        }
        TIDY(config);

      }/*gap start loop*/
    }/*refl start loop*/
    TIDY(doPars);
  }/*file loop*/

  return;
}/*fitTLSparams*/


/*##########################################*/
/*TLS error function*/

int tlsError(int numb,int npar,double *p,double *deviates,double **derivs,void *private)
{
  int i=0,n=0;
  float tE=0;     /*waveform energies*/
  float *tlsWave=NULL;
  float *makeOneTLSwaveform(overStruct *,double *,int);
  overStruct *data=NULL;

  data=(overStruct *)private;
  n=data->thisTLS;    /*copy to a local variable to save screen space*/

  /*make a waveform*/
  tlsWave=makeOneTLSwaveform(data,p,n);

  /*count up total energies*/
  tE=0.0;
  for(i=0;i<data->dArr[n].nALS;i++)tE+=tlsWave[i];
  if(tE<=0.0)tE=0.000000001;   /*put large error to avoid NaN*/

  /*fprintf(stdout,"Estimates %.10f %.10f\n",p[0],p[1]);*/
  /*assess deviates*/
  for(i=0;i<data->dArr[n].nALS;i++)deviates[i]=(double)(data->dArr[n].als[i]-tlsWave[i]*data->dArr[n].aE/tE);

  data=NULL;
  TIDY(tlsWave);
  return(0);
}/*tlsError*/


/*##########################################*/
/*set bounds for TLS optimisation*/

mp_par *tlsBounds(int nParams,char fixGap)
{
  mp_par *mparS=NULL;

  if(!(mparS=(mp_par *)calloc(nParams,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

  /*appRefl*/
  mparS[0].fixed=0;
  mparS[0].side=2;
  mparS[0].limited[0]=1;
  mparS[0].limits[0]=0.00000001;
  mparS[0].limited[1]=1;
  mparS[0].limits[1]=2.0;
  mparS[0].step=0.001;
  /*minGap*/
  if(fixGap==0){
    mparS[1].fixed=0;
    mparS[1].side=2;
    mparS[1].limited[0]=1;
    mparS[1].limits[0]=0.0001;
    mparS[1].limited[1]=1;
    mparS[1].limits[1]=1.0;
    mparS[1].step=0.001;
  }else{
    mparS[1].fixed=1;
  }

  return(mparS);
}/*tlsBounds*/


/*##########################################*/
/*smooth a wave with these structures*/

float *smoothOptWave(float *wave,overStruct *data,int numb)
{
  int i=0,j=0;
  int place=0;
  float *smoothed=NULL;
  float weight=0;

  smoothed=falloc(data->dArr[numb].nALS,"",0);

  for(i=0;i<data->dArr[numb].nALS;i++){
    smoothed[i]=weight=0.0;
    for(j=0;j<data->pBins;j++){
      place=(data->maxBin-j)+i;
      if((place>=0)&&(place<data->dArr[numb].nALS)){
        smoothed[i]+=wave[place]*data->pulse[j];
        weight+=data->pulse[j];
      }
    }
    if(weight>0.0)smoothed[i]/=weight;
  }/*bin loop*/

  return(smoothed);
}/*smoothOptWave*/


/*################################################*/
/*make TLS waveforms for denoising optimisation*/

void makeTLSwaves(overStruct *data)
{
  int numb=0;
  int i=0,bin=0;
  int nX=0,nY=0;
  float beamRad=0,rimRes=0;
  float *smoothOptWave(float *,overStruct *,int);
  float area=0,count=0;
  double vect[3];
  char **rImage=NULL;
  void makeBinImage(double *,float *,float *,char *,int,float,double,float,float,int,int,float,lidVoxPar *);
  void waveFromImageOld(char **,float **,int,int,int);

  beamRad=0.165;
  rimRes=0.001;   /*1 mm*/
  nX=nY=(int)(2.0*beamRad/rimRes);

  /*loop over files*/
  for(numb=0;numb<data->dimage->nFiles;numb++){
    /*first make range images*/
    rImage=chChalloc(data->dArr[numb].nALS,"",0);
    for(i=0;i<data->dArr[numb].nALS;i++)rImage[i]=challoc(nX*nY,"",0);

    for(i=0;i<data->dArr[numb].nTLS;i++){
      bin=data->dArr[numb].tls[i].bin;
      if((bin<0)||(bin>=data->dArr[numb].nALS)){
        fprintf(stderr,"BValls %d of %d\n",bin,data->dArr[numb].nALS);
        exit(1);
      }
      vect[0]=data->dArr[numb].tls[i].x;
      vect[1]=data->dArr[numb].tls[i].y;
      vect[2]=data->dArr[numb].tls[i].z;
      makeBinImage(&(vect[0]),&area,&count,&(rImage[bin][0]),data->dArr[numb].nALS,\
                   data->dArr[numb].tls[i].gap,(double)data->dArr[numb].tls[i].r,\
                   data->dArr[numb].tls[i].refl,beamRad,nX,nY,rimRes,&data->dArr[numb].tlsPar);
    }/*tls point loop*/

    /*trace through to make waveform*/
    data->dArr[numb].tlsWave=fFalloc(3,"tls wave",0);
    for(i=1;i<3;i++)data->dArr[numb].tlsWave[i]=falloc(data->dArr[numb].nALS,"tls wave",i+1);
    waveFromImageOld(rImage,&(data->dArr[numb].tlsWave[1]),data->dArr[numb].nALS,nX,nY);
    TTIDY((void **)rImage,data->dArr[numb].nALS);
    /*smooth to make ALS wave*/
    data->dArr[numb].tlsWave[0]=smoothOptWave(&(data->dArr[numb].tlsWave[1][0]),data,numb);

    /*smoothing for ALS optimisation*/
    if(data->dimage->testAtten==0)data->dArr[numb].smoothTLS=smooth(data->dimage->denWidth,\
                   data->dArr[numb].nALS,&(data->dArr[numb].tlsWave[1][0]),data->dimage->res);
    else                          data->dArr[numb].smoothTLS=smooth(data->dimage->denWidth,\
                   data->dArr[numb].nALS,&(data->dArr[numb].tlsWave[2][0]),data->dimage->res);

    /*determine wave start and end*/
    data->dArr[numb].tStart=-1;
    for(i=0;i<data->dArr[numb].nALS;i++){
      if(data->dArr[numb].smoothTLS[i]>TOL){
        if(data->dArr[numb].tStart<0)data->dArr[numb].tStart=i;
        data->dArr[numb].tEnd=i;
      }
    }/*tls wave bound determination*/

    /*TLS energy for denoising optimisation*/
    data->dArr[numb].tE=0.0;
    for(i=0;i<data->dArr[numb].nALS;i++)data->dArr[numb].tE+=data->dArr[numb].smoothTLS[i];
  }/*file loop*/

  return;
}/*makeTLSwaves*/


/*##########################################*/
/*make a TLS waveform*/

float *makeOneTLSwaveform(overStruct *data,double *params,int numb)
{
  int i=0,bin=0;
  int nX=0,nY=0;
  double vect[3];
  float **wave=NULL;
  float count=0,area=0;
  float rimRes=0;
  float beamRad=0;
  float *smoothed=NULL;
  float *smoothOptWave(float *,overStruct *,int);
  char **rImage=NULL;
  void makeBinImage(double *,float *,float *,char *,int,float,double,float,float,int,int,float,lidVoxPar *);
  void waveFromImageOld(char **,float **,int,int,int);


  /*set params*/
  data->dArr[numb].tlsPar.appRefl=(float)params[0];
  data->dArr[numb].tlsPar.minGap=(float)params[1];

  beamRad=0.165;
  rimRes=0.001;   /*1 mm*/
  nX=nY=(int)(2.0*beamRad/rimRes);


  /*first make range images*/
  rImage=chChalloc(data->dArr[numb].nALS,"",0);
  for(i=0;i<data->dArr[numb].nALS;i++)rImage[i]=challoc(nX*nY,"",0);

  for(i=0;i<data->dArr[numb].nTLS;i++){
    bin=data->dArr[numb].tls[i].bin;
    if((bin<0)||(bin>=data->dArr[numb].nALS)){
      fprintf(stderr,"BValls %d of %d\n",bin,data->dArr[numb].nALS);
      exit(1);
    }
    vect[0]=data->dArr[numb].tls[i].x;
    vect[1]=data->dArr[numb].tls[i].y;
    vect[2]=data->dArr[numb].tls[i].z;
    makeBinImage(&(vect[0]),&area,&count,&(rImage[bin][0]),data->dArr[numb].nALS,\
                 data->dArr[numb].tls[i].gap,(double)data->dArr[numb].tls[i].r,\
                 data->dArr[numb].tls[i].refl,beamRad,nX,nY,rimRes,&data->dArr[numb].tlsPar);
  }/*tls point loop*/
  /*trace through to make waveform*/
  wave=fFalloc(2,"tls wave",0);
  for(i=0;i<2;i++)wave[i]=falloc(data->dArr[numb].nALS,"tls wave",i+1);
  waveFromImageOld(rImage,wave,data->dArr[numb].nALS,nX,nY);
  TTIDY((void **)rImage,data->dArr[numb].nALS);


  /*convolve with ALS pulse*/
  smoothed=smoothOptWave(&(wave[0][0]),data,numb);

  TTIDY((void **)wave,2);
  return(smoothed);
}/*makeOneTLSwaveform*/


/*##########################################*/
/*set parameters*/

void setPars(overStruct *data)
{
  int i=0;

  /*loop over multiple files*/
  for(i=0;i<data->dimage->nFiles;i++){
    /*to optimise*/
    data->dArr[i].tlsPar.appRefl=0.7;      /*apparent TLS reflectance*/
    data->dArr[i].tlsPar.minGap=0.2;       /*minimum gap fraction to apply*/

    /*fixed for now*/
    data->dArr[i].tlsPar.beamTanDiv=sin(0.35/(2.0*1000.0))/cos(0.35/(2.0*1000.0));  /*tan half of 0.35 mrad*/
    data->dArr[i].tlsPar.beamRad=0.035;     /*7 mm diameter beam start*/
    data->dArr[i].tlsPar.minRefl=14.0;     /*minimum refletance value to scale between 0 and 1*/
    data->dArr[i].tlsPar.maxRefl=33.57;    /*maximum refletance value to scale between 0 and 1*/
  }/*file loop*/

  return;
}/*setPars*/


/*###############################################*/
/*make a waveform from a point cloud image*/

void waveFromImageOld(char **rImage,float **wave,int numb,int rNx,int rNy)
{
  int i=0,j=0;
  char doneIt=0;

  for(i=0;i<numb;i++)wave[0][i]=wave[1][i]=0.0;

  /*turn range images into a waveforms*/
  for(i=rNx*rNy-1;i>=0;i--){ /*image x-y loop*/
    doneIt=0;
    for(j=0;j<numb;j++){ /*image bin loop*/
      if(rImage[j][i]>0){
        if(doneIt==0){
          wave[0][j]+=1.0; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/
          doneIt=1;
        }/*first hit only*/
        wave[1][j]+=1.0; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/  /*hits per bin*/
      }
    }/*range image loop*/
  }/*range image bin loop*/
  return;
}/*waveFromImageOld*/


/*the end*/
/*####################################*/

