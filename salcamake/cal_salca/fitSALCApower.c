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

  /*preset calibration parameters*/
  float m;
  float c;
  float k;
  float a;
  float k2;
  float k3;

  /*parameter to fit*/
  float P0;   /*product of power and range function*/
}control;


/*##############################################*/
/*data structure*/

typedef struct{
  float *r;    /*movable range*/
  float *rho;  /*movable rho*/
  float *DN;   /*movable DN*/
  int numb;   /*number of movable data points*/

  /*fitted parameters*/
  float m;
  float c;
  float k;
  float a;
  float k2;
  float k3;
  float P0;   /*outgoing laser power*/

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
  data->k2=dimage->k2;
  data->k3=dimage->k3;
  data->P0=dimage->P0;

  return;
}/*initialGuess*/


/*##############################################*/
/*write parameters*/

void writeParams(dataStruct *data,control *dimage)
{
  int i=0;
  int maxDN=0;
  float DN=0,RMSE=0;
  float calibrate(float,float,float,float,float,float,float,float,float);
  char namen[200];
  FILE *opoo=NULL;

  maxDN=2000;
  //LUT=makeLUT(data,maxDN);


  sprintf(namen,"%s.par",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"m %f\nc %f\nk %f\n",data->m,data->c,data->k);
  fprintf(opoo,"a %f\nk2 %f\nk3 %f\n",data->a,data->k2,data->k3);
  fprintf(opoo,"P0 %f\n",data->P0);
  RMSE=0.0;
  for(i=0;i<data->numb;i++){
    DN=calibrate(data->r[i],data->rho[i],data->m,data->c,data->k,data->a,data->k2,data->k3,data->P0);
    RMSE+=(DN-data->DN[i])*(DN-data->DN[i]);
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
    DN=calibrate(data->r[i],data->rho[i],data->m,data->c,data->k,data->a,data->k2,data->k3,data->P0);
    fprintf(opoo,"%f %f %f %f %f\n",DN,data->DN[i],data->r[i],data->rho[i],data->P0);
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
  nPars=1;
  pars=dalloc(nPars,"parameters",0);
  pars[0]=data->P0;

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
  data->P0=(float)pars[0];

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

  /*P0*/
  mparS[0].fixed=0;
  mparS[0].side=3;
  mparS[0].limited[0]=1;
  mparS[0].limits[0]=0.0;
  mparS[0].limited[1]=1;
  mparS[0].limits[1]=100000.0;
  mparS[0].step=0.2;

  return(mparS);
}/*fitBounds*/


/*##############################################*/
/*error funtion*/

int fitError(int numb,int npar,double *pars,double *deviates,double **derivs,void *private)
{
  int i=0;
  float DN=0;
  float m=0,c=0,k=0;
  float a=0,k2=0,k3=0;
  float P0=0;
  float arg=0,expo=0;
  float rhoEff=0;
  float calibrate(float,float,float,float,float,float,float,float,float);
  float salcaRange(float,float,float,float);
  dataStruct *data=NULL;

  data=(dataStruct *)private;

  m=data->m;
  c=data->c;
  k=data->k;
  a=data->a;
  k2=data->k2;
  k3=data->k3;
  P0=(float)pars[0];

  for(i=0;i<numb;i++){
    DN=calibrate(data->r[i],data->rho[i],m,c,k,a,k2,k3,P0);

    deviates[i]=data->DN[i]-DN;

    if(derivs){ /*all negative as it's DN(truth)-DN(model)*/
      rhoEff=P0*data->rho[i]*salcaRange(a,k2,k3,data->r[i]);
      arg=k*rhoEff;
      if(arg<700.0)expo=exp(-1.0*arg);
      else         expo=0.0;

      if(derivs[0])derivs[0][i]=-1.0*data->rho[i]*salcaRange(a,k2,k3,data->r[i])*((m*rhoEff+c)*k*expo+m*(1.0-expo));
    }
  }/*point loop*/

  data=NULL;
  return(0);
}/*fitError*/


/*##############################################*/
/*generate a DN*/

float calibrate(float r,float rho,float m,float c,float k,float a,float k2,float k3,float P0)
{
  float DN=0,rhoEff=0;
  float nonLin=0;
  float arg=0;
  float salcaRange(float,float,float,float);

  rhoEff=P0*rho*salcaRange(a,k2,k3,r);
  arg=k*rhoEff;
  if(arg<730.0)nonLin=1.0-exp(-1.0*arg);
  else         nonLin=1.0;
  DN=(m*rhoEff+c)*nonLin;

  return(DN);
}/*calibrate*/


/*##############################################*/
/*salca range function*/

float salcaRange(float a,float k2,float k3,float r)
{
  float fR=0;
  float arg=0,shortR=0;

  arg=k2*pow(r,k3);
  if(arg<730.0)shortR=1.0-exp(-1.0*arg);
  else         shortR=1.0;
  fR=shortR/pow(r,a);

  return(fR);
}/*salcaRange*/


/*##############################################*/
/*read data*/

dataStruct *readData(control *dimage)
{
  int i=0;
  char line[300],temp1[100],temp2[100];
  char temp3[300],temp4[100],temp5[100];
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

fprintf(stdout,"There are %d\n",data->numb);

  /*allocate*/
  data->r=falloc(data->numb,"range",0);
  data->rho=falloc(data->numb,"rho",0);
  data->DN=falloc(data->numb,"DN",0);

  /*rewind*/
  if(fseek(ipoo,0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  i=0;
  while(fgets(line,300,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){  /*read data*/
      if(sscanf(line,"%s %s %s %s %s",temp1,temp2,temp3,temp4,temp5)==5){
        data->r[i]=atof(temp1);
        data->rho[i]=atof(temp5);
        data->DN[i]=atof(temp2);

        i++;
      }
    }
  }


  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

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

  /*initial guesses*/
  dimage->m=15000.0;
  dimage->c=500.0;
  dimage->k=200.0;
  dimage->a=2.0;  
  dimage->k2=2.0;
  dimage->k3=2.0;
  dimage->P0=1.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-m",2)){
        checkArguments(1,i,argc,"-m");
        dimage->m=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-k2",3)){
        checkArguments(1,i,argc,"-k2");
        dimage->k2=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-k3",3)){
        checkArguments(1,i,argc,"-k3");
        dimage->k3=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-a",2)){
        checkArguments(1,i,argc,"-a");
        dimage->a=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-k",2)){
        checkArguments(1,i,argc,"-k");
        dimage->k=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-c",2)){
        checkArguments(1,i,argc,"-c");
        dimage->c=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;    input filename\n-outRoot root;  output name root\n\nInitial guesses:\n-m m;\n-c c;\n-k k;\n-a a;\n-k2 k2;\n-k3 k3;\n\n");
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
