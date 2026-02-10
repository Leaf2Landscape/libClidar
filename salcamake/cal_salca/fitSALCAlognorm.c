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

  /*calibration parameters*/
  float m;
  float c;
  float k;

  /*initial guesses*/
  float a;
  float mu;
  float sigma;

  /*resolutions to set weightings by*/
  float rRes;   /*range resolution*/
  float rhoRes; /*reflectance resolution*/
}control;


/*##############################################*/
/*data structure*/

typedef struct{
  float *r;    /*movable range*/
  float *rho;  /*movable rho*/
  float *DN;   /*movable DN*/
  float *weight;  /*fitting weights*/
  int numb;   /*number of movable data points*/

  /*fitted parameters*/
  float m;
  float c;
  float k;
  float a;
  float mu;
  float sigma;
  float *P0;   /*outgoing laser power*/

  /*specificed parameters*/
  float *P0r;

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
  fitCalParams(data,dimage);

  /*output parameters*/
  writeParams(data,dimage);

  /*test fit*/
  testFit(data,dimage);

  /*tidy up arrays*/
  if(data){
    TIDY(data->r);
    TIDY(data->rho);
    TIDY(data->DN);
    TIDY(data->P0);
    TIDY(data->P0r);
    TIDY(data->weight);
    data->dimage=NULL;
    TIDY(data);
  }
  TIDY(dimage);
  return(0);
}/*main*/


/*##############################################*/
/*initial parameter guess*/

void initialGuess(dataStruct *data,control *dimage)
{

  data->m=dimage->m;
  data->c=dimage->c;
  data->k=dimage->k;
  data->a=dimage->a;
  data->mu=dimage->mu;
  data->sigma=dimage->sigma;

  return;
}/*initialGuess*/


/*##############################################*/
/*write parameters*/

void writeParams(dataStruct *data,control *dimage)
{
  int i=0;
  float DN=0,RMSE=0;
  float calibrate(float,float,float,float,float,float,float,float,float);
  char namen[200];
  FILE *opoo=NULL;


  sprintf(namen,"%s.par",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"m %f\nc %f\nk %f\na %f\nmu %f\nsigma %f\n",data->m,data->c,data->k,data->a,data->mu,data->sigma);
  RMSE=0.0;
  for(i=0;i<data->numb;i++){
    DN=calibrate(data->r[i],data->rho[i],data->m,data->c,data->k,data->a,data->mu,data->sigma,data->P0[i]);
    RMSE+=(DN-data->DN[i])*(DN-data->DN[i])*data->weight[i];
  }
  RMSE=sqrt(RMSE/(float)data->numb);
  fprintf(opoo,"RMSE %f\n",RMSE);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Parameters written to %s\n",namen);

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
  float calibrate(float,float,float,float,float,float,float,float,float);
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.dat",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"# 1 cal DN, 2 true DN, 3 range, 4 rho\n");
  for(i=0;i<data->numb;i++){
    DN=calibrate(data->r[i],data->rho[i],data->m,data->c,data->k,data->a,data->mu,data->sigma,data->P0[i]);
    fprintf(opoo,"%f %f %f %f %f\n",DN,data->DN[i],data->r[i],data->rho[i],data->P0[i]);
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
  int nPars=0;
  int fitError(int,int,double *,double *,double **,void *);
  double *pars=NULL;
  mp_par *mparS=NULL;
  mp_par *fitBounds(int,control *);
  mp_result *result=NULL;
  mp_config *config=NULL;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int fitCheck=0;


  /*copy initial estimates*/
  nPars=6;
  pars=dalloc(nPars,"parameters",0);
  pars[0]=data->m;  /*m*/
  pars[1]=data->c;  /*c*/
  pars[2]=data->k;   /*k*/
  pars[3]=data->a;    /*a*/
  pars[4]=data->mu;    /*mu*/
  pars[5]=data->sigma;    /*sigma*/


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
  result->resid=dalloc(data->numb,"",0);
  result->xerror=dalloc(nPars,"",0);
  result->covar=dalloc(nPars*nPars,"",0);

  /*set bounds*/
  mparS=fitBounds(nPars,dimage);

  /*optimise parameters*/
  fitCheck=mpfit(fitError,data->numb,nPars,pars,mparS,config,(void *)data,result);
  if(fitCheck<0){
    fprintf(stderr,"Check %d\n",fitCheck);
    exit(1);
  }

  /*copy data to structure*/
  data->m=(float)pars[0];
  data->c=(float)pars[1];
  data->k=(float)pars[2];
  data->a=(float)pars[3];
  data->mu=(float)pars[4];
  data->sigma=(float)pars[5];

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
  mp_par *mparS=NULL;

  if(!(mparS=(mp_par *)calloc(nPars,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

  /*m*/
  mparS[0].fixed=0;
  mparS[3].side=3;
  mparS[3].limited[0]=1;
  mparS[3].limits[0]=50.0;
  mparS[3].limited[1]=1;
  mparS[3].limits[1]=1000.0;
  mparS[3].step=10.0;
  /*c*/
  mparS[1].fixed=0;
  mparS[3].side=3;
  mparS[3].limited[0]=1;
  mparS[3].limits[0]=15.0;
  mparS[3].limited[1]=1;
  mparS[3].limits[1]=500.0;
  mparS[3].step=5.0;
  /*k*/
  mparS[2].fixed=0;
  mparS[3].side=3;
  mparS[3].limited[0]=1;
  mparS[3].limits[0]=0.05;
  mparS[3].limited[1]=1;
  mparS[3].limits[1]=100.0;
  mparS[3].step=1.0;
  /*a*/
  mparS[3].fixed=0;
  mparS[3].side=3;
  mparS[3].limited[0]=1;
  mparS[3].limits[0]=1.5;
  mparS[3].limited[1]=1;
  mparS[3].limits[1]=10.0;
  mparS[3].step=0.5;
  /*mu*/
  mparS[4].fixed=0;
  mparS[4].side=3;
  mparS[4].limited[0]=1;
  mparS[4].limits[0]=0.0;
  mparS[4].limited[1]=1;
  mparS[4].limits[1]=100.0;
  mparS[4].step=0.5;
  /*sigma*/
  mparS[5].fixed=0;
  mparS[5].side=3;
  mparS[5].limited[0]=1;
  mparS[5].limits[0]=0.1;
  mparS[5].limited[1]=1;
  mparS[5].limits[1]=10.0;
  mparS[5].step=0.25;

  return(mparS);
}/*fitBounds*/


/*##############################################*/
/*error funtion*/

int fitError(int numb,int npar,double *pars,double *deviates,double **derivs,void *private)
{
  int i=0;
  float DN=0;
  float m=0,c=0,k=0,a=0,mu=0,sigma=0,P0=0;
  float calibrate(float,float,float,float,float,float,float,float,float);
  float x=0,arg=0,expo=0;
  float arg2=0,argRoot=0,expo2=0;
  float salcaRange(float,float,float,float);
  float dfR=0,dDNdfR=0;
  dataStruct *data=NULL;
  void fitCalParams(dataStruct *,control *);

  data=(dataStruct *)private;

  m=(float)pars[0];
  c=(float)pars[1];
  k=(float)pars[2];
  a=(float)pars[3];
  mu=(float)pars[4];
  sigma=(float)pars[5];

  for(i=0;i<numb;i++){
    P0=data->P0r[i]/salcaRange(a,mu,sigma,8.0);
    data->P0[i]=P0;

    DN=calibrate(data->r[i],data->rho[i],m,c,k,a,mu,sigma,P0);

    deviates[i]=(data->DN[i]-DN)*data->weight[i];


    if(derivs){ /*all negative as it's DN(truth)-DN(model)*/
      x=data->rho[i]*P0*salcaRange(a,mu,sigma,data->r[i]);
      arg=k*x;
      if(arg<700.0)expo=exp(-1.0*arg);
      else         expo=0.0; 
      if(derivs[0])derivs[0][i]=-1.0*x*(1.0-expo)*data->weight[i];    /*m*/
      if(derivs[1])derivs[1][i]=-1.0*(1.0-expo)*data->weight[i];      /*c*/
      if(derivs[2])derivs[2][i]=-1.0*x*(m*x+c)*expo*data->weight[i];  /*k*/

      /*range parameters*/
      argRoot=(log(data->r[i])-mu)/sigma;
      arg2=argRoot*argRoot/2.0;
      if(arg2<730.0)expo2=exp(-1.0*arg2);
      else          expo2=0.0;
      dDNdfR=(m*x+c)*k*data->rho[i]*P0*expo+m*data->rho[i]*P0*(1.0-expo);

      if(derivs[3]){  /*a*/
        dfR=log(data->r[i])*pow(data->r[i],-1.0*a)*expo2;
        derivs[3][i]=-1.0*dfR*dDNdfR*data->weight[i];
      }
      if(derivs[4]){  /*mu*/
        dfR=expo2*argRoot/(sigma*pow(data->r[i],a));
        derivs[4][i]=-1.0*dfR*dDNdfR*data->weight[i];
      }

      if(derivs[5]){  /*sigma*/
        dfR=arg2*2.0/(sigma*pow(data->r[i],a))*expo2;
        derivs[5][i]=-1.0*dfR*dDNdfR*data->weight[i];
      }
    }
  }/*point loop*/

  data=NULL;
  return(0);
}/*fitError*/


/*##############################################*/
/*generate a DN*/

float calibrate(float r,float rho,float m,float c,float k,float a,float mu,float sigma,float P0)
{
  float DN=0,rhoEff=0;
  float nonLin=0;
  float arg=0;
  float salcaRange(float,float,float,float);

  rhoEff=P0*rho*salcaRange(a,mu,sigma,r);
  arg=k*rhoEff;
  if(arg<730.0)nonLin=1.0-exp(-1.0*arg);
  else         nonLin=1.0;
  DN=(m*rhoEff+c)*nonLin;

  return(DN);
}/*calibrate*/


/*##############################################*/
/*salca range function*/

float salcaRange(float a,float mu,float sigma,float r)
{
  float fR=0;
  float arg=0,shortR=0;

  arg=(log(r)-mu)/sigma;
  arg*=arg/2.0;  /*square and divide by 2*/
  if(arg<730.0)shortR=exp(-1.0*arg);
  else         shortR=0.0;
  fR=shortR/pow(r,a);

  return(fR);
}/*salcaRange*/


/*##############################################*/
/*read data*/

dataStruct *readData(control *dimage)
{
  int i=0;
  int *nIn=NULL,maxIn=0;
  int place=0,nR=0,nRho=0;
  int rhoBin=0,rBin=0;
  float minR=0,maxR=0,minRho=0,maxRho=0;
  char line[300],temp1[100],temp2[100];
  char temp3[300],temp4[100],temp5[100];
  char temp6[300],temp7[100];
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
  data->numb=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
        data->numb++;
      }
    }
  }/*data reading*/


  /*allocate*/
  data->r=falloc(data->numb,"range",0);
  data->rho=falloc(data->numb,"rho",0);
  data->DN=falloc(data->numb,"DN",0);
  data->P0r=falloc(data->numb,"P0",0);
  data->P0=falloc(data->numb,"P0",0);

  /*rewind*/
  if(fseek(ipoo,0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  minR=minRho=1000.0;
  maxR=maxRho=-100.0;
  i=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){  /*read data*/
      if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
        data->r[i]=atof(temp1);
        data->rho[i]=atof(temp5);
        data->DN[i]=atof(temp2);
        data->P0r[i]=atof(temp7);

        if(data->r[i]<minR)minR=data->r[i];
        if(data->r[i]>maxR)maxR=data->r[i];
        if(data->rho[i]<minRho)minRho=data->rho[i];
        if(data->rho[i]>maxRho)maxRho=data->rho[i];

        i++;
      }
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  data->dimage=dimage;

  /*determine fitting weights*/
  nR=(int)((maxR-minR)/dimage->rRes)+2;
  nRho=(int)((maxRho-minRho)/dimage->rhoRes)+2;
  nIn=ialloc(nR*nRho,"",0);
  for(i=nR*nRho-1;i>=0;i--)nIn[i]=0;
  maxIn=0;
  for(i=0;i<data->numb;i++){
    rBin=(int)((data->r[i]-minR)/dimage->rRes+0.5);
    rhoBin=(int)((data->rho[i]-minRho)/dimage->rhoRes+0.5);
    place=rhoBin*nR+rBin;
    nIn[place]++;
    if(nIn[place]>maxIn)maxIn=nIn[place];
  }

  data->weight=falloc(data->numb,"weight",0);
  for(i=0;i<data->numb;i++){
    rBin=(int)((data->r[i]-minR)/dimage->rRes+0.5);
    rhoBin=(int)((data->rho[i]-minRho)/dimage->rhoRes+0.5);
    place=rhoBin*nR+rBin;
    data->weight[i]=(float)maxIn/(float)nIn[place];
  }

  TIDY(nIn);
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

  /*initial guesses*/
  dimage->m=15000.0;  /*m*/
  dimage->c=500.0;  /*c*/
  dimage->k=200.0;   /*k*/
  dimage->a=2.0;    /*a*/
  dimage->mu=4.0;    /*mu*/
  dimage->sigma=2.0;  /*sigma*/

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
      }else if(!strncasecmp(argv[i],"-mu",3)){
        checkArguments(1,i,argc,"-mu");
        dimage->mu=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-m",2)){
        checkArguments(1,i,argc,"-m");
        dimage->m=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sigma",6)){
        checkArguments(1,i,argc,"-sigma");
        dimage->sigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-a",2)){
        checkArguments(1,i,argc,"-a");
        dimage->a=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-k",2)){
        checkArguments(1,i,argc,"-k");
        dimage->k=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-c",2)){
        checkArguments(1,i,argc,"-c");
        dimage->c=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rRes",5)){
        checkArguments(1,i,argc,"-rRes");
        dimage->rRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rhoRes",7)){
        checkArguments(1,i,argc,"-rhoRes");
        dimage->rhoRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;    input filename\n-outRoot root;  output name root\n\nInitial guesses:\n-m m;\n-c c;\n-k k;\n-a a;\n-mu mu;\n-sigma sigma;\n\n-rRes res;     range weighting res\n-rhoRes res;     reflectance weighting resolution\n\n");
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
