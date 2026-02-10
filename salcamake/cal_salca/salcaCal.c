#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "mpfit.h"


/*##############################################*/
/*control structure*/

typedef struct{
  /*input output*/
  char inNamen[200];
  char outRoot[200];

  /*switches*/
  char fixPar;    /*optimise for laser power only*/
  char fixP0;     /*fix laser power*/
  char useFix;    /*use fixed range only*/
  char fixRange;  /*fix range parameter switch*/
  char fixShortR; /*fix short range*/
  char useAll;

  /*initial guesses*/
  float m;
  float c;
  float k;
  float a;
  float k2;
  float *P0;

  /*resolutions to set weightings by*/
  float rRes;   /*range resolution*/
  float rhoRes; /*reflectance resolution*/
}control;


/*##############################################*/
/*data structure*/

typedef struct{
  float *mR;    /*movable range*/
  float *mRho;  /*movable rho*/
  float *mDN;   /*movable DN*/
  float *mWeight; /*weights for fitting*/
  int nMove;   /*number of movable data points*/
  float *fR;   /*fixed range*/
  float *fRho; /*fixed rho*/
  float *fDN;  /*fixed DN*/
  float *fWeight; /*weights for fitting*/
  int nFix;    /*number of fixed data points*/

  /*file labels*/
  int nFiles;  /*number of scan files*/
  int *fFile;  /*fixed file label*/
  int *mFile;  /*movable file label*/

  /*pointers for data to fit*/
  int nUse;
  int *uFile;
  float *uDN;
  float *uR;
  float *uRho;
  float *uWeight;

  /*fitted parameters*/
  float m;
  float c;
  float k;
  float a;
  float k2;
  float *P0;   /*outgoing laser power*/

  /*pointer to control*/
  control *dimage;
}dataStruct;


/*##############################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;
  dataStruct *readData(control *);
  void fitCalParams(dataStruct *,control *);
  void testFit(dataStruct *,control *);
  void writeParams(dataStruct *,control *);
  void initialGuess(dataStruct *,control *);

  /*read commands*/
  dimage=readCommands(argc,argv);

  /*read data*/
  data=readData(dimage);

  /*set initial guess*/
  initialGuess(data,dimage);

  /*perform optimsation*/
  dimage->useFix=0;
  dimage->fixP0=1;
  dimage->fixPar=0;
  fitCalParams(data,dimage);

  /*output parameters*/
  writeParams(data,dimage);

  /*test fit*/
  testFit(data,dimage);

  /*tidy up arrays*/
  if(data){
    TIDY(data->mR);
    TIDY(data->mRho);
    TIDY(data->mDN);
    TIDY(data->mWeight);
    TIDY(data->fWeight);
    TIDY(data->fR);
    TIDY(data->P0);
    TIDY(data->fRho);
    TIDY(data->fDN);
    TIDY(data->fFile);
    TIDY(data->mFile);
    data->uFile=NULL;
    data->uDN=NULL;
    data->uWeight=NULL;
    data->uR=NULL;
    data->uRho=NULL;
    data->dimage=NULL;
    TIDY(data);
  }
  if(dimage){
    TIDY(dimage->P0);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##############################################*/
/*initial parameter guess*/

void initialGuess(dataStruct *data,control *dimage)
{
  int i=0;

  data->m=dimage->m;
  data->c=dimage->c;
  data->k=dimage->k;
  data->a=dimage->a;
  data->k2=dimage->k2;
  data->P0=falloc(data->nFiles,"P0",0);
  for(i=0;i<data->nFiles;i++){  /*set laser power by filter*/
    data->P0[i]=1.0; //dimage->P0[i];
  }

  return;
}/*initialGuess*/


/*##############################################*/
/*write parameters*/

void writeParams(dataStruct *data,control *dimage)
{
  int i=0;
  int maxDN=0;
  float DN=0,RMSE=0;
  float rhoRMSE=0;
  float calibrate(float,float,float,float,float,float,float,float);
  float *makeLUT(dataStruct *,int);
  float rhoEff=0,rho=0;
  char namen[200];
  float *LUT=NULL;   /*look up table for calibration*/
  FILE *opoo=NULL;

  maxDN=2000;
  //LUT=makeLUT(data,maxDN);


  sprintf(namen,"%s.par",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"m %f\nc %f\nk %f\na %f\nk2 %f\n",data->m,data->c,data->k,data->a,data->k2);
  for(i=0;i<data->nFiles;i++)fprintf(opoo,"P0 file %d %f\n",i,data->P0[i]);
  RMSE=rhoRMSE=0.0;
  for(i=0;i<data->nUse;i++){
    //rhoEff=LUT[(int)data->uDN[i]];
    rho=(rhoEff/data->P0[data->uFile[i]])*(pow(data->uR[i],data->a)/(1.0-exp(-1.0*data->k2*data->uR[i])));
    rhoRMSE+=(rho-data->uRho[i])*(rho-data->uRho[i]);

    DN=calibrate(data->uR[i],data->uRho[i],data->m,data->c,data->k,data->a,data->k2,data->P0[data->uFile[i]]);
    RMSE+=((DN-data->uDN[i])*(DN-data->uDN[i])*data->uWeight[i]);
  }
  RMSE=sqrt(RMSE/(float)data->nUse);
  rhoRMSE=sqrt(rhoRMSE/(float)data->nUse);
  fprintf(opoo,"RMSE %f\n",RMSE);
  fprintf(opoo,"rhoRMSE %f\n",rhoRMSE);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Parameters written to %s\n",namen);
  TIDY(LUT);

  return;
}/*writeParams*/


/*##############################################*/
/*make LUT*/

float *makeLUT(dataStruct *data,int maxDN)
{
  int i=0;
  int *nIn=NULL;
  float *LUT=NULL;
  float rho=0,res=0;
  float DN=0,lastDN=0;
  float expo=0;

  LUT=falloc(maxDN,"LUT",0);
  nIn=ialloc(maxDN,"nIn",0);

  res=0.000003;
  rho=0.0;
  lastDN=0.0;
  do{
    if((data->k*rho)<700.0)expo=exp(-1.0*data->k*rho);
    else                   expo=0.0;
    DN=(data->m*rho+data->c)*(1.0-expo);
    if((DN-lastDN)>1.0){
      fprintf(stderr,"Increased resolution needed %f %f\n",DN,lastDN);
      exit(1);
    }
    lastDN=DN;
    if(DN>=maxDN)break;
    LUT[(int)DN]+=rho;
    nIn[(int)DN]++;

    rho+=res;
  }while(DN<maxDN);

  for(i=0;i<maxDN;i++){
    if(nIn[i]>0)LUT[i]/=(float)nIn[i];
  }

  TIDY(nIn);
  return(LUT);
}/*makeLUT*/


/*##############################################*/
/*test fit*/

void testFit(dataStruct *data,control *dimage)
{
  int i=0;
  float DN=0;
  float calibrate(float,float,float,float,float,float,float,float);
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.dat",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"# 1 cal DN, 2 true DN, 3 range, 4 rho\n");
  for(i=0;i<data->nUse;i++){
    DN=calibrate(data->uR[i],data->uRho[i],data->m,data->c,data->k,data->a,data->k2,data->P0[data->uFile[i]]);
    fprintf(opoo,"%f %f %f %f %f\n",DN,data->uDN[i],data->uR[i],data->uRho[i],data->P0[data->uFile[i]]);
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Test fit written to %s\n",namen);

  return;
}/*testFit*/


/*##############################################*/
/*do the calibration*/

void fitCalParams(dataStruct *data,control *dimage)
{
  int i=0;
  int nPars=0;
  int fitError(int,int,double *,double *,double **,void *);
  double *pars=NULL;
  mp_par *mparS=NULL;
  mp_par *fitBounds(int,control *);
  mp_result *result=NULL;
  mp_config *config=NULL;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int fitCheck=0;

  if(dimage->useFix){
    data->nUse=data->nFix;
    data->uFile=data->fFile;
    data->uDN=data->fDN;
    data->uR=data->fR;
    data->uRho=data->fRho;
    data->uWeight=data->fWeight;
  }else{
    data->nUse=data->nMove;
    data->uFile=data->mFile;
    data->uDN=data->mDN;
    data->uR=data->mR;
    data->uRho=data->mRho;
    data->uWeight=data->mWeight;
  }


  /*copy initial estimates*/
  nPars=5+data->nFiles;
  pars=dalloc(nPars,"parameters",0);
  pars[0]=data->m;  /*m*/
  pars[1]=data->c;  /*c*/
  pars[2]=data->k;   /*k*/
  pars[3]=data->a;    /*a*/
  pars[4]=data->k2;    /*k2*/
  for(i=5;i<nPars;i++)pars[i]=data->P0[i-5];  /*assume all lasers balanced*/


  /*option space*/
  if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  //config->nofinitecheck=1;
  //config->stepfactor=10.0;
  config->maxiter=1000;

  /*result space*/
  if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
    fprintf(stderr,"error in mpfit structure.\n");
    exit(1);
  }
  result->resid=dalloc(data->nUse,"",0);
  result->xerror=dalloc(nPars,"",0);
  result->covar=dalloc(nPars*nPars,"",0);

  /*set bounds*/
  mparS=fitBounds(nPars,dimage);

  /*optimise parameters*/
  fitCheck=mpfit(fitError,data->nUse,nPars,pars,mparS,config,(void *)data,result);
  if(fitCheck<0){
    fprintf(stderr,"Check %d\n",fitCheck);
    exit(1);
  }

  /*copy data to structure*/
  data->m=(float)pars[0];
  data->c=(float)pars[1];
  data->k=(float)pars[2];
  data->a=(float)pars[3];
  data->k2=(float)pars[4];
  if(data->dimage->fixP0==0)for(i=5;i<nPars;i++)data->P0[i-5]=(float)pars[i];

  /*fprintf(stdout,"%f %f %f %f %f\n",data->m,data->c,data->k,data->a,data->k2);
  for(i=0;i<data->nFiles;i++)fprintf(stdout,"%f ",data->P0[i]);
  fprintf(stdout,"\n");*/

  if(result){
    TIDY(result->resid);
    TIDY(result->xerror);
    TIDY(result->covar);
    TIDY(result);
  }
  TIDY(mparS);
  TIDY(config);
  TIDY(pars);
  return;
}/*fitCalParams*/


/*##############################################*/
/*set bounds*/

mp_par *fitBounds(int nPars,control *dimage)
{
  int i=0;
  mp_par *mparS=NULL;

  if(!(mparS=(mp_par *)calloc(nPars,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

  if(dimage->fixPar==0){
    /*m*/
    mparS[0].fixed=0;
    mparS[0].side=3;
    mparS[0].limited[0]=1;
    mparS[0].limits[0]=0.0;
    mparS[0].limited[1]=1;
    mparS[0].limits[1]=100000.0;
    mparS[0].step=10.0;
    /*c*/
    mparS[1].fixed=0;
    mparS[1].side=3;
    mparS[1].limited[0]=1;
    mparS[1].limits[0]=0.0;
    mparS[1].limited[1]=1;
    mparS[1].limits[1]=10000.0;
    mparS[1].step=10.0;
    /*k*/
    mparS[2].fixed=0;
    mparS[2].side=3;
    mparS[2].limited[0]=1;
    mparS[2].limits[0]=0.01;
    mparS[2].limited[1]=1;
    mparS[2].limits[1]=2000.0;
    mparS[2].step=4.0;
    /*range parameters*/
    if(dimage->fixRange==0){
      /*a*/
      mparS[3].fixed=0;
      mparS[3].side=3;
      mparS[3].limited[0]=1;
      mparS[3].limits[0]=0.2;
      mparS[3].limited[1]=1;
      mparS[3].limits[1]=4.0;
      mparS[3].step=0.5;
      /*k2*/
      if(dimage->fixShortR==0){
        mparS[4].fixed=0;
        mparS[4].side=3;
        mparS[4].limited[0]=1;
        mparS[4].limits[0]=0.02;
        mparS[4].limited[1]=1;
        mparS[4].limits[1]=40.0;
        mparS[4].step=0.02;
      }else mparS[4].fixed=1;
    }else{
      mparS[3].fixed=1;
      mparS[4].fixed=1;
    }
  }else{   /*lock all but laser power*/
    for(i=0;i<5;i++)mparS[i].fixed=1;
  }
  /*laser powers*/
  if(dimage->fixP0==0){
    for(i=5;i<nPars;i++){
      mparS[i].fixed=0;
      mparS[i].side=3;
      mparS[i].limited[0]=1;
      mparS[i].limits[0]=0.001;
      mparS[i].limited[1]=1;
      mparS[i].limits[1]=40.0;
      mparS[i].step=0.2;
    }
    //mparS[5+7].fixed=1;  /*set all relative to first. 7th is brightest*/
  }else{
    for(i=5;i<nPars;i++)mparS[i].fixed=1;
  }

  return(mparS);
}/*fitBounds*/


/*##############################################*/
/*error funtion*/

int fitError(int numb,int npar,double *pars,double *deviates,double **derivs,void *private)
{
  int i=0;
  float DN=0;
  float m=0,c=0,k=0,a=0,k2=0,P0=0;
  float calibrate(float,float,float,float,float,float,float,float);
  float x=0;
  control dimage,*tempDimage=NULL;
  dataStruct *data=NULL;
  void fitCalParams(dataStruct *,control *);

  data=(dataStruct *)private;

  m=(float)pars[0];
  c=(float)pars[1];
  k=(float)pars[2];
  a=(float)pars[3];
  k2=(float)pars[4];

  /*determine laser power*/
  if(data->dimage->fixP0==1){
    /*fprintf(stdout,"Laser power\n");*/
    /*new conntrol variable*/
    tempDimage=data->dimage;
    data->dimage=&dimage;
    dimage.useFix=1;
    dimage.fixP0=0;
    dimage.fixPar=1;
    dimage.fixRange=0;
    dimage.useAll=0;
    data->m=m;
    data->c=c;
    data->k=k;
    data->a=a;
    data->k2=k2;
    /*fprintf(stdout,"before %f %f\n",data->k2,data->P0[0]);*/
    fitCalParams(data,&dimage);
    /*fprintf(stdout," after %f %f\n",data->k2,data->P0[0]);*/
    /*fprintf(stdout,"Parameters\n");*/
    data->nUse=data->nMove;
    data->uFile=data->mFile;
    data->uDN=data->mDN;
    data->uR=data->mR;
    data->uRho=data->mRho;
    data->uWeight=data->mWeight;
    data->dimage=tempDimage;
    tempDimage=NULL;
    /*copy laser powers back*/
    for(i=6;i<npar;i++)pars[i]=data->P0[i-6];
  }


  for(i=0;i<numb;i++){
    P0=pars[5+data->uFile[i]]; //data->P0[data->uFile[i]];

    DN=calibrate(data->uR[i],data->uRho[i],m,c,k,a,k2,P0);

    if(data->uWeight==data->mWeight){
      if(i>=data->nMove){
        fprintf(stderr,"Error %d of %d of %d\n",i,data->nMove,numb);
        exit(1);
      }
    }

    deviates[i]=(data->uDN[i]-DN)*data->uWeight[i];;


    if(derivs){ /*all negative as it's DN(truth)-DN(model)*/
      x=data->uRho[i]*P0*(1.0-exp(-1.0*k2*data->uR[i]))/pow(data->uR[i],a);
      if(derivs[0])derivs[0][i]=-1.0*x*(1.0-exp(-1.0*k*x));
      if(derivs[1])derivs[1][i]=-1.0*(1.0-exp(-1.0*k*x));
      if(derivs[2])derivs[2][i]=-1.0*x*(m*x+c)*exp(-1.0*k*x);
      if(derivs[3])derivs[3][i]=-1.0*(x*log(data->uR[i])*(m*(1.0-exp(-1.0*k*x))+k*exp(-1.0*k*x)*(m*x+c)));
      if(derivs[4])derivs[4][i]=-1.0*(data->uRho[i]*P0*(exp(-1.0*k2*data->uR[i])/pow(data->uR[i],a-1.0))*(m*(1.0-exp(-1.0*k*x))+k*exp(-1.0*k*x)*(m*x+c)));
      if(derivs[5+data->uFile[i]])derivs[5+data->uFile[i]][i]=-1.0*(x/P0)*(m*(1.0-exp(-1.0*k*x))+k*exp(-1.0*k*x)*(m*x+c));
    }
  }/*point loop*/

  data=NULL;
  return(0);
}/*fitError*/


/*##############################################*/
/*generate a DN*/

float calibrate(float r,float rho,float m,float c,float k,float a,float k2,float P0)
{
  float DN=0,rhoEff=0;
  float shortR=0,nonLin=0;
  float arg=0;

  arg=k2*r;
  if(arg<730.0)shortR=1.0-exp(-1.0*arg);
  else         shortR=1.0;
  rhoEff=P0*rho*shortR/pow(r,a);
  arg=k*rhoEff;
  if(arg<730.0)nonLin=1.0-exp(-1.0*arg);
  else         nonLin=1.0;
  DN=(m*rhoEff+c)*nonLin;

  return(DN);
}/*calibrate*/


/*##############################################*/
/*read data*/

dataStruct *readData(control *dimage)
{
  int i=0,j=0,k=0;
  int maxIn=0;
  int place=0,rBin=0,rhoBin=0;
  int *nIn=NULL,nR=0,nRho=0,total=0;
  float minR=0,maxR=0,minRho=0,maxRho=0;
  float minRf=0,maxRf=0,minRhoF=0,maxRhoF=0;
  char line[300],temp1[100],temp2[100];
  char temp3[300],temp4[100],temp5[100];
  char temp6[300],temp7[100],temp8[100];
  char temp9[100],temp10[100];
  char newFile=0;
  char **fileList=NULL;
  dataStruct *data=NULL;
  FILE *ipoo=NULL;

  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error in data structure allocation.\n");
    exit(1);
  }

  if((ipoo=fopen(dimage->inNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",dimage->inNamen);
    exit(1);
  }

  /*count number of points*/
  data->nFix=data->nMove=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
        if(!strncasecmp(temp4,"m",1))data->nMove++;
        else if(!strncasecmp(temp4,"f",1))data->nFix++;
      }
    }
  }/*data reading*/

  dimage->P0=falloc(data->nFix,"P0",0);

  /*movable targets*/
  data->mR=falloc(data->nMove,"range",0);
  data->mRho=falloc(data->nMove,"rho",0);
  data->mDN=falloc(data->nMove,"DN",0);
  data->mFile=ialloc(data->nMove,"mFile",0);
  /*fixed targets*/
  data->fR=falloc(data->nFix,"range",0);
  data->fRho=falloc(data->nFix,"rho",0);
  data->fDN=falloc(data->nFix,"DN",0);
  data->fFile=ialloc(data->nFix,"fFile",0);

  /*rewind*/
  if(fseek(ipoo,0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  minR=minRho=minRf=minRhoF=1000.0;
  maxR=maxRho=maxRf=maxRhoF=-100.0;
  i=j=0;
  data->nFiles=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){  /*read data*/
      if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
        if(!strncasecmp(temp4,"m",1)){  /*movable target*/
          data->mR[i]=atof(temp1);
          data->mRho[i]=atof(temp5);
          data->mDN[i]=atof(temp2);

          if(data->mR[i]<minR)minR=data->mR[i];
          if(data->mR[i]>maxR)maxR=data->mR[i];
          if(data->mRho[i]<minRho)minRho=data->mRho[i];
          if(data->mRho[i]>maxRho)maxRho=data->mRho[i];

          /*see if file is new*/
          newFile=1;
          for(k=0;k<data->nFiles;k++){
            if(!strncasecmp(temp7,fileList[k],strlen(fileList[k]))){
              newFile=0;
              break;
            }
          }
          data->mFile[i]=k;

          if(newFile){
            if(!(fileList=(char **)realloc(fileList,(data->nFiles+1)*sizeof(char *)))){
              fprintf(stderr,"Balls\n");
              exit(1);
            }
            fileList[data->nFiles]=challoc(strlen(temp7)+1,"",0);
            strcpy(fileList[data->nFiles],temp7);
            dimage->P0[data->nFiles]=1.0; //atof(temp10);
            data->nFiles++;
          }

          i++;
        }else if(!strncasecmp(temp4,"f",1)){  /*fixed target*/
          data->fR[j]=atof(temp1);
          data->fRho[j]=atof(temp5);
          data->fDN[j]=atof(temp2);

          if(data->fR[j]<minRf)minRf=data->fR[j];
          if(data->fR[j]>maxRf)maxRf=data->fR[j];
          if(data->fRho[j]<minRhoF)minRhoF=data->fRho[j];
          if(data->fRho[j]>maxRhoF)maxRhoF=data->fRho[j];

          /*see if file is new*/
          newFile=1;
          for(k=0;k<data->nFiles;k++){
            if(!strncasecmp(temp7,fileList[k],strlen(fileList[k]))){
              newFile=0;
              break;
            }
          }
          data->fFile[j]=k;

          if(newFile){
            if(!(fileList=(char **)realloc(fileList,(data->nFiles+1)*sizeof(char *)))){
              fprintf(stderr,"Balls\n");
              exit(1);
            }
            fileList[data->nFiles]=challoc(strlen(temp7)+1,"",0);
            strcpy(fileList[data->nFiles],temp7);
            dimage->P0[data->nFiles]=1.0; //atof(temp10);
            data->nFiles++;
          }

          j++;
        }
      }
    }
  }


  /*sort through moving and determine fitting weights*/
  /*weights for movable targets*/
  nR=(int)((maxR-minR)/dimage->rRes)+2;
  nRho=(int)((maxRho-minRho)/dimage->rhoRes)+2;
  nIn=ialloc(nR*nRho,"",0);

  for(i=nR*nRho-1;i>=0;i--)nIn[i]=0;
  total=0;
  maxIn=0;
  for(i=0;i<data->nMove;i++){
    rBin=(int)((data->mR[i]-minR)/dimage->rRes+0.5);
    rhoBin=(int)((data->mRho[i]-minRho)/dimage->rhoRes+0.5);
    place=rhoBin*nR+rBin;
if((place<0)||(place>=(nR*nRho))){
  fprintf(stderr,"Error %d of %d and %d of %d\n",rBin,nR,rhoBin,nRho);
  exit(1);
}
    nIn[place]++;
    if(nIn[place]>maxIn)maxIn=nIn[place];
    total++;
  }

  data->mWeight=falloc(data->nMove,"weight",0);
  for(i=0;i<data->nMove;i++){
    rBin=(int)((data->mR[i]-minR)/dimage->rRes+0.5);
    rhoBin=(int)((data->mRho[i]-minRho)/dimage->rhoRes+0.5);
    place=rhoBin*nR+rBin;
    data->mWeight[i]=(float)maxIn/(float)nIn[place];
  }
  TIDY(nIn);
  /*weights for fixed targets*/
  nR=(int)((maxRf-minRf)/dimage->rRes)+2;
  nRho=(int)((maxRhoF-minRhoF)/dimage->rhoRes)+2;
  nIn=ialloc(nR*nRho,"",0);

  for(i=nR*nRho-1;i>=0;i--)nIn[i]=0;
  total=0;
  maxIn=0;
  for(i=0;i<data->nFix;i++){
    rBin=(int)((data->fR[i]-minRf)/dimage->rRes+0.5);
    rhoBin=(int)((data->fRho[i]-minRhoF)/dimage->rhoRes+0.5);
    place=rhoBin*nR+rBin;
    nIn[place]++;
    if(nIn[place]>maxIn)maxIn=nIn[place];
    total++;
  }

  data->fWeight=falloc(data->nFix,"weight",0);
  for(i=0;i<data->nFix;i++){
    rBin=(int)((data->fR[i]-minRf)/dimage->rRes+0.5);
    rhoBin=(int)((data->fRho[i]-minRhoF)/dimage->rhoRes+0.5);
    place=rhoBin*nR+rBin;
    data->fWeight[i]=(float)maxIn/(float)nIn[place];
  }
  TIDY(nIn);



  TTIDY((void **)fileList,data->nFiles);
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  if(dimage->useAll){  /*combine datasets*/
    if(!(data->mR=(float *)realloc(data->mR,(data->nFix+data->nMove)*sizeof(float)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
    if(!(data->mFile=(int *)realloc(data->mFile,(data->nFix+data->nMove)*sizeof(int)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
    if(!(data->mRho=(float *)realloc(data->mRho,(data->nFix+data->nMove)*sizeof(float)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
    for(i=data->nMove;i<data->nFix+data->nMove;i++){
      data->mR[i]=data->fR[i-data->nMove];
      data->mRho[i]=data->fRho[i-data->nMove];
      data->mFile[i]=data->fFile[i-data->nMove];
    }
    data->nMove+=data->nFix;
  }/*combine datasets*/

  data->dimage=dimage;

  return(data);
}/*readData*/


/*##############################################*/
/*read commands*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in control allocation.\n");
    exit(1);
  }
  strcpy(dimage->inNamen,"/mnt/urban-bess/shancock_work/data/SALCA/raw/brisbane/brisbaneTest.1.dat");
  strcpy(dimage->outRoot,"fittedSALCA");
  dimage->fixPar=0;  /*optimise for all*/
  dimage->useFix=0;  /*use moving range only*/
  dimage->fixP0=0;   /*optiimse laser power*/
  dimage->fixRange=0;
  dimage->fixShortR=0;
  dimage->useAll=0;

  /*initial guesses*/
  dimage->m=15000.0;  /*m*/
  dimage->c=500.0;  /*c*/
  dimage->k=200.0;   /*k*/
  dimage->a=1.5;    /*a*/
  dimage->k2=1.4;    /*k2*/
  dimage->P0=NULL;

  /*weighting resolutions*/
  dimage->rRes=100.0;   /*range resolution*/
  dimage->rhoRes=100.0; /*reflectance resolution*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-fixPar",7)){
        dimage->fixPar=1;
      }else if(!strncasecmp(argv[i],"-useFix",7)){
        dimage->useFix=1;
      }else if(!strncasecmp(argv[i],"-fixP0",6)){
        dimage->fixP0=1;
      }else if(!strncasecmp(argv[i],"-m",2)){
        checkArguments(1,i,argc,"-m");
        dimage->m=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-k2",3)){
        checkArguments(1,i,argc,"-k2");
        dimage->k2=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-a",2)){
        checkArguments(1,i,argc,"-a");
        dimage->a=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-k",2)){
        checkArguments(1,i,argc,"-k");
        dimage->k=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-c",2)){
        checkArguments(1,i,argc,"-c");
        dimage->c=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fixRange",9)){
        dimage->fixRange=1;
      }else if(!strncasecmp(argv[i],"-useAll",7)){
        dimage->useAll=1;
      }else if(!strncasecmp(argv[i],"-rRes",5)){
        checkArguments(1,i,argc,"-rRes");
        dimage->rRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rhoRes",7)){
        checkArguments(1,i,argc,"-rhoRes");
        dimage->rhoRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fixShortR",10)){
        dimage->k2=1000.0;
        dimage->fixShortR=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;    input filename\n-outRoot root;  output name root\n-useFix;        use fixed range\n-fixPar;        lock all parameters except laser power\n-fixP0;     fix laser power\n-fixRange;      fix range parameters\n-useAll;      use all points in fitting\n\nInitial guesses:\n-m m;\n-c c;\n-k k;\n-a a;\n-k2 k2;\n\n-rRes res;     range weighting res\n-rhoRes res;     reflectance weighting resolution\n-fixShortR;     fix short range effect\n\n");
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
/*##############################################*/
