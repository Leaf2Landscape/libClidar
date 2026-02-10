#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"


/*####################################*/
/*# Fits function to a single return #*/
/*# read from an ASCII file. The     #*/
/*# files are extracted elsewhere.   #*/
/*####################################*/


#define MAX_ITER 30

typedef struct{
  float *x;   /*range*/
  char *data; /*intensity*/
  int numb;   /*number of bins*/
}feats;


void gaussErr(float,float *,float *,float *,int);
float gauss(float,float,float);

int main(int argc,char **argv)
{
  int i=0;
  int sBin=0,width=0;
  float *setData(int,float *);
  float *params=NULL;
  float *temp=NULL;
  float *fitGauss(float *,float *,int,float);
  float *filterData(char *,float *,int, int *,int *);
  float minErr=0;
  void writeFitted(char *,float *,float *,int,float *,int);
  feats *feat=NULL;
  feats *readData(char *);
  char namen[200];  /*input filename*/
  char outNamen[200];


  strcpy(namen,"/mnt/geodesy38/nsh103/SALCA/analysis/return_shape/return.fail.zen.1087.az.0.bin.12.band.0.dat");
  strcpy(outNamen,"teast.dat");

  minErr=0.001;


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(namen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;    input filename\n-output name;   output filename\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*read the command line*/



  /*load data*/
  feat=readData(namen);

  /*filter data and determine width*/
  temp=filterData(feat->data,feat->x,feat->numb,&width,&sBin);

  /*optimise*/
  params=fitGauss(&(feat->x[sBin-1]),temp,width,minErr);


  /*test the fit*/
  writeFitted(outNamen,feat->x,temp,width,params,1);


  /*test output*/
  fprintf(stdout,"Retrieved mu %f sig %f A %f\n",params[1],params[3],params[2]);

  TIDY(params);
  TIDY(temp);
  if(feat){
    TIDY(feat->data);
    TIDY(feat->x);
    TIDY(feat);
  }

  return(0);
}/*main*/


/*###############################################*/
/*read data*/

feats *readData(char *namen)
{
  int i=0;
  feats *feat=NULL;
  char line[200];
  char temp1[100],temp2[100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  if(!(feat=(feats *)calloc(1,sizeof(feats)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }

  feat->numb=0;
  /*count number of lines*/
  while(fgets(line,200,ipoo)!=NULL){
    if(sscanf(line,"%s",temp1)==1){
      if(strncasecmp(temp1,"#",1)){
        feat->numb++;
      }
    }
  }/*line counting*/

  feat->data=challoc(feat->numb+1,"data",0);
  feat->x=falloc(feat->numb+1,"range",0);

  /*rewind*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read the data*/
  i=1;
  while(fgets(line,200,ipoo)!=NULL){
    if(sscanf(line,"%s %s",temp1,temp2)==2){
      if(strncasecmp(temp1,"#",1)){
        if(i>feat->numb){
          fprintf(stderr,"Error\n");
          exit(1);
        }
        feat->x[i]=atof(temp1);
        feat->data[i]=atoi(temp2);
        i++;
      }
    }
  }/*data reading*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(feat);
}/*readData*/


/*###############################################*/
/*extract relevant data and determine width*/

float *filterData(char *y,float *x,int numb, int *width,int *sBin)
{
  int i=0;
  int start=0,end=0;
  int count=0;
  char thresh=0;
  char backSig=0;
  char saturated=0;
  char *useMap=NULL;
  float *temp=NULL;
  float *newX=NULL;

  (*width)=0;
  thresh=-110;
  backSig=-115;
  start=-1;
  useMap=challoc(numb+1,"useable map",0);

  /*determine width and check for saturation*/
  saturated=0;
  for(i=1;i<=numb;i++){
    if(y[i]>thresh){
      if(start<0)start=i;
      end=i;
      if(saturated&&(y[i]<127)){  /*in a saturated feature*/
        if((y[i]>y[i-1])&&(y[i]>y[i+1]))saturated=0;  /*a local maximum is the end*/
      }
      if(!saturated){
        useMap[i]=1;
        (*width)++;
      }else useMap[i]=0;
      if(y[i]>=127)saturated=1;
    }
  }
  start--;  /*step to 0 signal*/
  end++;
  (*width)+=2; /*=end-start+1;*/   /*add zero point before and after*/
  (*sBin)=start;

  /*allocate space and copy data*/
  temp=falloc((*width)+1,"temporary data",0);
  newX=falloc((*width)+1,"new range",0);
  count=1;
  for(i=start;i<=end;i++){
    if(y[i]<=thresh){
       temp[count]=0.0;
       newX[count]=x[i];
       count++;
    }else if(useMap[i]){
      temp[count]=(float)y[i]-(float)backSig;
      newX[count]=x[i];
      count++;
    }
  }
  for(i=1;i<=*width;i++){
    x[i]=newX[i];
fprintf(stdout,"data %f %f\n",x[i],temp[i]);
  }


  TIDY(newX);
  TIDY(useMap);
  return(temp);
}/*filterData*/


/*###############################################*/
/*fit Gaussians*/

float *fitGauss(float *x,float *y,int numb,float minErr)
{
  int i=0;
  int nGauss=0,nParams=0;
  int *pSwitch=NULL;   /*switch to optimise each parameter*/
  int indicator=0;
  float chisq=0;
  float alambda=0;
  float **covar=NULL;
  float **alpha=NULL;
  float *sig=NULL;
  float *params=NULL;
  float *initialGuess(float *,float *,int,int);
  int mrqminMine(float *,float *,float *,int,float *,int *,int,float **,float **,float *,void (*)(float,float *,float *,float *,int),float *);

  nGauss=1;
  nParams=3*nGauss;
  params=initialGuess(x,y,numb,nGauss);


  minErr=(float)numb*20.0;   /*base minimum acceptable error on number of bins*/


  if(numb==3){  /*one point above noise, fix max and width*/
    params[3]=0.17;

    return(params);
  //}else if(numb==4){  /*two points above noise, fix width*/
  //  pSwitch[3]=0;
  }
  fprintf(stdout,"Width %d\n",numb);

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

  pSwitch[3]=0;
  params[3]=0.102;

  /*initialise*/
  alambda=-1.0;
  indicator=mrqminMine(x,y,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
  fprintf(stdout,"Initialiased, status %d\n",indicator);

  if(indicator==0){
    /*do the iterations*/
    alambda=0.001;
    i=0;
    do{
      indicator=mrqminMine(x,y,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
      if(indicator)break;
      fprintf(stdout,"Iteration %d error %f aiming for %f\n",i,chisq,minErr);
      i++;
      if(i>=MAX_ITER){
        fprintf(stderr,"Maximum iterations exceeded, %d\n",i);
        break;
      }
    }while(chisq>minErr);

    if(indicator==0){
      /*final call to get bits and pieces*/
      alambda=0.0;
      indicator=mrqminMine(x,y,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
      fprintf(stdout,"Finalised, status %d\n",indicator);
    }
  }

  /*force widths to be positive*/
  for(i=0;i<nGauss;i++)params[3*i+3]=fabs(params[3*i+3]);

  if(indicator){
    for(i=1;i<=nParams;i++)params[i]=-1.0;
  }


  TTIDY((void **)covar,nParams+1);
  TTIDY((void **)alpha,nParams+1);
  TIDY(pSwitch);
  TIDY(sig);
  return(params);
}/*fitGauss*/


/*###############################################*/
/*write out fitted waveform*/

void writeFitted(char *outNamen,float *x,float *y,int numb,float *params,int nGauss)
{
  int i=0,j=0;
  float f=0;
  float backSig=0;
  FILE *opoo=NULL;

  backSig=-115.0;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }

  for(i=1;i<=numb;i++){
    f=0.0;
    for(j=0;j<nGauss;j++){
      f+=params[3*j+2]*gauss(x[i],params[3*j+3],params[3*j+1]);
    }
    fprintf(opoo,"%f %f %f\n",x[i],y[i],f+backSig);
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);

  return;
}/*writeFitted*/


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

    *yfit+=A*gauss(x,sig,mu);
    dyda[1]+=A*(x-mu)/(sig*sig)*exp(arg);
    dyda[2]+=exp(arg);
    dyda[3]+=A*(x-mu)*(x-mu)/(sig*sig*sig)*exp(arg);
  }/*Gaussian loop*/

  return;
}/*gaussErr*/


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


  fprintf(stdout,"Initial guess mu %f A %f sig %f\n",params[1],params[2],params[3]);

  return(params);
}/*initialGuess*/


/*###############################################*/
/*set data*/

float *setData(int numb,float *x)
{
  int i=0;
  float *y=NULL;
  float mu=0,sig=0;
  float A=0;

  mu=6.623;
  sig=0.015;
  A=0.5;

  y=falloc(numb+1,"data",0);
  for(i=1;i<=numb;i++){
    y[i]=A*gauss(x[i],sig,mu);
  }
  fprintf(stdout,"Truth mu %f sig %f A %f\n",mu,sig,A);

  return(y);
}/*setData*/


/*###############################################*/
/*non normalised Gaussian*/

float gauss(float x,float sigma,float offset)
{
  float y=0;

  y=exp(-1.0*(x-offset)*(x-offset)/(2.0*sigma*sigma));

  return(y);
}/*gauss*/


/*the end*/
/*###############################################*/

