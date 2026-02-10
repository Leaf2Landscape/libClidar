#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"


/*###############################*/
/*# Creates a Gaussian return,  #*/
/*# resamples then fits to that #*/
/*# by Levenberg-Marquardt      #*/
/*###############################*/


#define TOL 0.00001

void gaussErr(float,float *,float *,float *,int);
float gauss(float,float,float);

int main(int argc,char **argv)
{
  int numb=0,i=0;
  int sBin=0,width=0;
  float mu=0,sig=0,A=0;
  float *y=NULL,*x=NULL;
  float *setData(int,float *,float,float,float);
  float *params=NULL;
  float *temp=NULL;
  float *fitGauss(float *,float *,int,float);
  float *filterData(float *,int, int *,int *);
  float *fitQuadratic(float *,float *,int);
  float minErr=0;
  float *quad=NULL;
  float sWidth=0,*smoothY=NULL;
  float *smooth(float *,int,float);
  void writeFitted(char *,float *,float *,int,float *,int);
  char smooBoth=0;

  minErr=0.001;

  /*return characteristics*/
  mu=6.623;
  sig=0.11;
  A=0.5;

  /*smoothing width*/
  sWidth=1.1;
  smooBoth=0;   /*don't smooth for Levenberg-Marquardt*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-A",2)){
        checkArguments(1,i,argc,"-A");
        A=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sig",4)){
        checkArguments(1,i,argc,"-sig");
        sig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-mu",3)){
        checkArguments(1,i,argc,"-mu");
        mu=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sWidth",7)){
        checkArguments(1,i,argc,"-sWidth");
        sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-smooBoth",9)){
        smooBoth=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-A A;    amplitude\n-sig sig;    width\n-mu x;    centre range\n-sWidth sig;    smoothing sidth\n-smooBoth;    smooth for L-M\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/


  /*load data*/
  numb=100;
  x=falloc(numb+1,"range",0);
  for(i=1;i<=numb;i++)x[i]=0.15*(float)i;
  y=setData(numb,x,mu,sig,A);

  /*smooth for quadratic fitting*/
  if(sWidth>TOL)smoothY=smooth(y,numb,sWidth);
  else          smoothY=&(y[0]);


  /*Extract relevant bins. filter data and determine width*/
  if(!smooBoth)temp=filterData(y,numb,&width,&sBin);
  else         temp=filterData(smoothY,numb,&width,&sBin);

  /*optimise*/
  params=fitGauss(&(x[sBin-1]),temp,width,minErr);
  if(!smooBoth)TIDY(temp);

  /*use D Jupp's method*/
  if(!smooBoth)temp=filterData(smoothY,numb,&width,&sBin);
  quad=fitQuadratic(&(x[sBin-1]),temp,width);

  /*test the fit*/
  writeFitted("teast.dat",x,y,numb,params,1);


  /*test output*/
  fprintf(stdout,"Retrieved mu %f sig %f A %f\n",params[1],params[3],params[2]);
  fprintf(stdout,"Quadratic numbers x %f y %f\n",quad[0],quad[1]);


  TIDY(params);
  TIDY(temp);
  TIDY(quad);
  TIDY(y);
  TIDY(x);

  return(0);
}/*main*/


/*###############################################*/
/*smooth a waveform*/

float *smooth(float *y,int numb,float sWidth)
{
  int i=0,j=0;
  int nSmoo=0;
  int place=0;
  float contN=0;
  float *smoothY=NULL;
  float *smoother=NULL;
  float *setSmoother(float,float,int *);

  smoothY=falloc(numb,"smoothed",0);
  smoother=setSmoother(0.15,sWidth,&nSmoo);

  for(i=0;i<numb;i++){
    smoothY[i]=0.0;
    contN=0.0;
    for(j=0;j<nSmoo;j++){
      place=i-j+nSmoo/2;
      if((place>=0)&&(place<numb)){
        smoothY[i]+=y[place]*smoother[j];
        contN+=smoother[j];
      }
    }
    if(contN>0.0)smoothY[i]/=contN;
  }

  TIDY(smoother);
  return(smoothY);
}/*smooth*/


/*###############################################*/
/*set up a smoothing function*/

float *setSmoother(float res,float sWidth,int *nSmoo)
{
  int i=0;
  float tol=0,mid=0,y=0;
  float *smoother=NULL;

  tol=0.0001;

  /*determine width*/
  i=0;
  do{
    y=gauss((float)i*res,sWidth,0.0);
    i++;
  }while(y>=tol);
  *nSmoo=i*2;

  mid=(float)i*res;
  smoother=falloc(*nSmoo,"smoother",0);
  for(i=0;i<*nSmoo;i++){
    smoother[i]=gauss((float)i*res,sWidth,mid);
  }

  return(smoother);
}/*setSmoother*/


/*###############################################*/
/*fit a quadratic using David Jupp's method*/

float *fitQuadratic(float *x,float *y,int width)
{
  int i=0;
  int mid=0;
  float max=0;
  float a=0,b=0,c=0;  /*quadratic coefficients*/
  float *quad=NULL;   /*x and y of maximum of quadratic*/
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

  quad=falloc(2,"quadratic",0);
  d0=(y[mid]-y[mid-1])/(x[mid]-x[mid-1]);
  d1=(y[mid+1]-y[mid])/(x[mid+1]-x[mid]);

  c=(d1-d0)/(x[mid+1]-x[mid-1]);
  b=d0-c*(x[mid]+x[mid-1]);
  a=y[mid]-b*x[mid]-c*x[mid]*x[mid];

  quad[0]=-1.0*b/(2.0*c);
  /*enforce bounds. Do not let the range stray beyond these three points*/
  if((quad[0])<x[mid-1])     quad[0]=x[mid-1];
  else if((quad[0])>x[mid+1])quad[0]=x[mid+1];

  quad[1]=a+b*quad[0]+c*quad[0]*quad[0];

  return(quad);
}/*fitQuadratic*/


/*###############################################*/
/*extract relevant data and determine width*/

float *filterData(float *y,int numb, int *width,int *sBin)
{
  int i=0;
  int start=0,end=0;
  float thresh=0;
  float *temp=NULL;

  (*width)=0.0;
  thresh=0.0001;
  start=-1;

  /*determine width*/
  for(i=1;i<=numb;i++){
    if(y[i]>thresh){
      if(start<0)start=i;
      end=i;
    }
  }
  start--;  /*step to 0 signal*/
  end++;
  (*width)=end-start+1;
  (*sBin)=start;

  /*allocate space and copy data*/
  temp=falloc((*width)+1,"temporary data",0);
  for(i=start;i<=end;i++){
    temp[i-start+1]=y[i];
  }

  return(temp);
}/*filterData*/


/*###############################################*/
/*fit Gaussians*/

float *fitGauss(float *x,float *y,int numb,float minErr)
{
  int i=0;
  int nGauss=0,nParams=0;
  int *pSwitch=NULL;   /*switch to optimise each parameter*/
  float chisq=0;
  float alambda=0;
  float **covar=NULL;
  float **alpha=NULL;
  float *sig=NULL;
  float *params=NULL;
  float *initialGuess(float *,float *,int,int);
  void mrqminMine(float *,float *,float *,int,float *,int *,int,float **,float **,float *,void (*)(float,float *,float *,float *,int),float *);

  nGauss=1;
  nParams=3*nGauss;
  params=initialGuess(x,y,numb,nGauss);


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


  /*initialise*/
  alambda=-1.0;
  mrqminMine(x,y,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);


  /*do the iterations*/
  alambda=0.001;
  i=0;
  do{
    mrqminMine(x,y,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);
    //fprintf(stdout,"Iteration %d error %f\n",i,chisq);
    i++;
  }while(chisq>minErr);


  /*final call to get bits and pieces*/
  alambda=0.0;
  mrqminMine(x,y,sig,numb,params,pSwitch,nParams,covar,alpha,&chisq,gaussErr,&alambda);

  /*force widths to be positive*/
  for(i=0;i<nGauss;i++)params[3*i+3]=fabs(params[3*i+3]);


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
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }

  for(i=1;i<=numb;i++){
    f=0.0;
    for(j=0;j<nGauss;j++){
      f+=params[3*j+2]*gauss(x[i],params[3*j+3],params[3*j+1]);
    }
    fprintf(opoo,"%f %f %f\n",x[i],y[i],f);
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
}/*gaussian fit error*/


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

float *setData(int numb,float *x,float mu,float sig,float A)
{
  int i=0;
  float *y=NULL;

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

