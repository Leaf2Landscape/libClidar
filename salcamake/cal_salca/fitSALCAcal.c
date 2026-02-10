#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "mpfit.h"


/*#########################*/
/*# Stripped down version #*/
/*# of salcaCal.c         #*/
/*# September 2015        #*/
/*#########################*/



/*##############################################*/
/*control structure*/

typedef struct{
  /*input output*/
  char inNamen[200];
  char outRoot[200];

  /*initial guesses*/
  float m;
  float c;
  float k;
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
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##############################################*/
/*initial parameter guess*/

void initialGuess(dataStruct *data,control *dimage)
{
  data->m=dimage->m;
  data->c=dimage->c;
  data->k=dimage->k;

  return;
}/*initialGuess*/


/*##############################################*/
/*write parameters*/

void writeParams(dataStruct *data,control *dimage)
{
  int i=0;
  float DN=0,RMSE=0;
  float calibrate(float,float,float,float,float,float,float,float);
  float meanDN=0;
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.par",dimage->outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"m %f\nc %f\nk %f\n",data->m,data->c,data->k);
  RMSE=0.0;
  for(i=0;i<data->numb;i++){
    DN=calibrate(data->r[i],data->rho[i],data->m,data->c,data->k,0.1,10000.0,1.0);
    meanDN+=data->DN[i];
    RMSE+=(DN-data->DN[i])*(DN-data->DN[i]);
  }
  RMSE=sqrt(RMSE/(float)data->numb);
  meanDN/=(float)data->numb;
  fprintf(opoo,"RMSE %f\nmeanDN %f\n",RMSE,meanDN);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Parameters written to %s\n",namen);

  return;
}/*writeParams*/


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
  for(i=0;i<data->numb;i++){
    DN=calibrate(data->r[i],data->rho[i],data->m,data->c,data->k,0.0,100000.0,1.0);
    fprintf(opoo,"%f %f %f %f %f\n",DN,data->DN[i],data->r[i],data->rho[i],1.0);
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
  nPars=3;
  pars=dalloc(nPars,"parameters",0);
  pars[0]=data->m;  /*m*/
  pars[1]=data->c;  /*c*/
  pars[2]=data->k;   /*k*/

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
  mp_par *mparS=NULL;

  if(!(mparS=(mp_par *)calloc(nPars,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

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

  return(mparS);
}/*fitBounds*/


/*##############################################*/
/*error funtion*/

int fitError(int numb,int npar,double *pars,double *deviates,double **derivs,void *private)
{
  int i=0;
  float DN=0;
  float m=0,c=0,k=0,P0=0;
  float calibrate(float,float,float,float,float,float,float,float);
  float x=0;
  dataStruct *data=NULL;
  void fitCalParams(dataStruct *,control *);

  data=(dataStruct *)private;

  m=(float)pars[0];
  c=(float)pars[1];
  k=(float)pars[2];

  for(i=0;i<numb;i++){
    P0=1.0;
    DN=calibrate(data->r[i],data->rho[i],m,c,k,0.0,1000000.0,P0);
    deviates[i]=(data->DN[i]-DN);

    if(derivs){ /*all negative as it's DN(truth)-DN(model)*/
      x=data->rho[i]*P0;
      if(derivs[0])derivs[0][i]=-1.0*x*(1.0-exp(-1.0*k*x));
      if(derivs[1])derivs[1][i]=-1.0*(1.0-exp(-1.0*k*x));
      if(derivs[2])derivs[2][i]=-1.0*x*(m*x+c)*exp(-1.0*k*x);
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

fprintf(stdout,"%d targets\n",data->numb);

  /*targets*/
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
fprintf(stdout,"%d %f %f %f\n",i,data->r[i],data->rho[i],data->DN[i]); 
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
  dimage->m=15000.0;  /*m*/
  dimage->c=500.0;  /*c*/
  dimage->k=200.0;   /*k*/

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
      }else if(!strncasecmp(argv[i],"-k",2)){
        checkArguments(1,i,argc,"-k");
        dimage->k=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-c",2)){
        checkArguments(1,i,argc,"-c");
        dimage->c=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;    input filename\n-outRoot root;  output name root\n\n\nInitial guesses:\n-m m;\n-c c;\n-k k;\n\n");
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
