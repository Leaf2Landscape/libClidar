#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"

/*##########################*/
/*# Reads raw SALCA binary #*/
/*# and outputs an ascii   #*/
/*# point cloud            #*/
/*##########################*/


/*function to optimise to*/
//void gaussErr(float,float *,float *,float *,int);
double *rotateVect(double *,double **);


#define MAX_ITER 30
#define TOL 0.00001

int sOffset;   /*wavestart offset*/
float backSig;  /*background signal*/
char thresh;   /*noise threshold*/
char satThresh; /*saturation threshold*/

typedef struct{ /*to hold calibration data*/
  float *LUT;   /*DN to reflectance array*/
  float a;      /*range decay exponent*/
  float k2;     /*short range parameter*/
  float k3;     /*short range parameter*/
  float minDN;  /*min DN in LUT*/
  float maxDN;  /*max DN in LUT*/
  float tran;   /*filter transmission*/
  int numb;     /*number of elements in LUT*/
}calibration;


/*######################################################################*/

typedef struct{ /*to hold control options*/
  /*options*/
  int coarsen;     /*aggragate beams*/
  char calibrate;  /*calibration switch*/
  char writeBin;   /*make a binary file rather than ASCII*/
  char joy;
  /*input outputs*/
  char inRoot[200];  /*input filename root*/
  char outRoot[200]; /*output filename root*/
  char calFile[100]; /*DN to refl calibration file*/
  char matName[200]; /*geolocation matrix name*/
  uint32_t shotN[2]; /*total number of shots written per file*/
  char transName[200];/*collocation matrix*/
  char globName[200]; /*geolocation matrix*/
  /*scan characteristics*/
  int nAz;
  int nZen;
  float azStep;
  float zStep;
  float maxZen;
  float azStart;
  float azSquint[2];  /*laser squint angles per beam*/
  float zenSquint[2]; /*laser squint angles per beam*/
  float omega;      /*mirror slope angle*/
  float mSquint;    /*motor squint angle*/
  float **zen;       /*true zenith per beam*/
  float **azOff;     /*azimuth offset per beam*/
  float res;        /*range resolution*/
  float *smoother;    /*smoothing function quad*/
  float sWidth;       /*smoothing function width quad*/
  float I0[2];        /*laser power and filter scalar*/
  int nSmoo;          /*length of smoothing function quad*/
  char doQuad;        /*do Jupp's quadratic fitting*/
  char calOldStyle;   /*for backwards compatability*/
  float maxR;     /*maximum range*/
  double scanCent[3];  /*scan centre coordinate*/
  /*filter strengths*/
  float filt;      /*filter strength (0, 0.6, 1.0 or 1.6)*/
  float filt1;     /*filter strength (0, 0.6, 1.0 or 1.6)*/
  float filt2;     /*filter strength (0, 0.6, 1.0 or 1.6)*/
  /*collocation matrices*/
  double **trans;  /*scan collocation matrix*/
  double **glob;   /*global collocation matrix*/
  double dZen;     /*change in zenith*/
  double dAz;      /*change in azimuth*/
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
  char *data=NULL;  /*waveform data*/
  char *readData(char *,int,int *,int *,int *,control *);
  double **setMatrix(char *);
  double **combineGeolocation(double **,char *);
  void pointOut(char *,int,int,FILE *,float,int,int,int,char,calibration *,control *,int);
  void closeFiles(FILE **,char *,char);
  void tidyCal(calibration *);
  calibration *readCalibration(char *,float,float,float,control *);
  calibration *cal=NULL;
  control *options=NULL;
  control *readCommands(int,char **);
  FILE **openOutput(char *,char);
  FILE **opoo=NULL;              /*pointer to output files*/
  void translateSquint(control *);
  void setSquint(control *,int);
  void writeBinEnding(uint32_t *,double *,FILE **);


  /*the defaults are for the hedge test*/
  nBins=1200;

  sOffset=7;   /*1.05 m*/
  backSig=-116.0;
  thresh=-110;

  /*read the command line*/
  options=readCommands(argc,argv);


  /*set number of shots per line*/
  numb=(int)(options->maxZen/0.059375);

  /*read calibration if needed*/
  if(options->calibrate)cal=readCalibration(options->calFile,options->filt,options->filt1,options->filt2,options); /*set up calibration bits*/

  /*set scan array dimensions*/
  nBins=(int)((150.0+options->maxR)/0.15);  /*150m of 1550 nm plus maxR of 1040 nm*/

  /*band start and end bins for sampling*/
  start[0]=sOffset;   /*1.05 m to avoid the outgoing pulse*/
  end[0]=(int)(options->maxR/0.15);
  start[1]=1000+sOffset;   /*7 to avoid outgoing pulses*/
  end[1]=nBins;

  /*open the two output files, one for each band*/
  opoo=openOutput(options->outRoot,options->writeBin);


  /*if coarsened*/
  if(options->coarsen<1)options->coarsen=1;  /*so zero means uncoarsened*/
  options->nAz/=options->coarsen;
  numb/=options->coarsen;
  /*set the scan offsets*/
  options->azStep*=(float)options->coarsen*M_PI/180.0; /*convert to radians as well*/
  options->nZen/=options->coarsen;
  options->zStep=0.001047198*(float)options->coarsen;  /*fixed for SALCA*/

  /*set up squint angles*/
  translateSquint(options);
  setSquint(options,numb);

  /*read geolocation matrix if needed*/
  options->trans=setMatrix(options->transName);
  options->trans=combineGeolocation(options->trans,options->globName);
  /*use as scan centre*/
  for(i=0;i<3;i++)options->scanCent[i]=options->trans[i][3];

  /*loop through azimuth files and write out points*/
  for(i=0;i<=options->nAz;i++){ /*az loop*/
    data=readData(options->inRoot,i,&numb,&nBins,&length,options); /*read binary data into an array*/

    if(data){ /*does the file exist?*/
      for(band=0;band<2;band++){   /*loop through wavebands*/
        pointOut(data,numb,nBins,opoo[band],(float)i*options->azStep+options->azStart+(float)options->dAz,start[band],end[band],i,options->calibrate,&(cal[band]),options,band);
      }/*band loop*/
      TIDY(data);   /*clean up arrays as we go along*/
    }/*does file exist?*/
  }/*az loop*/

  /*write file end if binary*/
  if(options->writeBin)writeBinEnding(options->shotN,options->scanCent,opoo);
  closeFiles(opoo,options->outRoot,options->writeBin);   /*close output files*/

  if(options){
    if(options->calibrate)tidyCal(cal);
    TTIDY((void **)options->azOff,2);
    TTIDY((void **)options->zen,2);
    TTIDY((void **)options->trans,4);
    TTIDY((void **)options->glob,4);
    TIDY(options->smoother);
    TIDY(options);
  }
  return(0);
}/*main*/


/*##########################################*/
/*combine geolocation and collocation matrices*/

double **combineGeolocation(double **trans,char *globName)
{
  int i=0,j=0,k=0;
  double **tempMat=NULL,**newMat=NULL;
  void readMatrix(char *,double **);

  /*allocate space*/
  newMat=dDalloc(4,"new matrix",0);
  tempMat=dDalloc(4,"new matrix",0);
  for(i=0;i<4;i++){
    tempMat[i]=dalloc(4,"new matrix",i+1);
    newMat[i]=dalloc(4,"new matrix",i+1);
  }

  /*do we need to do anything?*/
  if(strncasecmp(globName,"none",4)){   /*read a file*/
    readMatrix(globName,tempMat);
    /*apply global matrix*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        newMat[i][j]=0.0;
        for(k=0;k<4;k++){
          newMat[i][j]+=tempMat[i][k]*trans[k][j];
        }
      }
    }
  }else{
    for(i=0;i<4;i++){
      for(j=0;j<4;j++)newMat[i][j]=trans[i][j];
    }
  }

  TTIDY((void **)trans,4);
  TTIDY((void **)tempMat,4);
  return(newMat);
}/*combineGeolocation*/


/*##########################################*/
/*read a collocation matrix*/

void readMatrix(char *namen,double **mat)
{
  int i=0,j=0;
  FILE *ipoo=NULL;
  char line[400];
  char *token=NULL;

  /*open input*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file list \"%s\"\n",namen);
    exit(1);
  }

  /*read file*/
  i=0;
  while(fgets(line,399,ipoo)!=NULL){
    if(i>=4){
      fprintf(stderr,"Badly formatted matrix file %s %d \"%s\"\n",namen,i,line);
      exit(1);
    }
    token=strtok(line," ");
    j=0;
    while(token!=NULL) {
      if(j<4)mat[i][j]=atof(token);
      token=strtok(NULL," ");
      j++;
    }
    i++;
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*readMatrix*/


/*##########################################*/
/*set a collocation matrix*/

double **setMatrix(char *transName)
{
  int i=0,j=0;
  double **mat=NULL;
  void readMatrix(char *,double **);

  /*allocate space*/
  mat=dDalloc(4,"matrix",0);
  for(i=0;i<4;i++)mat[i]=dalloc(4,"matrix",0);

  /*read or unity?*/
  if(strncasecmp(transName,"none",4)){   /*read a file*/
    readMatrix(transName,mat);
  }else{   /*make unity matrix*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        if(i!=j)mat[i][j]=0.0;
        else    mat[i][j]=1.0;
      }
    }
  }

  return(mat);
}/*setMatrix*/


/*##########################################*/
/*output a point cloud*/

void pointOut(char *data,int numb,int nBins,FILE *opoo,float az,int start,int end,int azInd,char calibrate,calibration *cal,control *options,int band)
{
  int i=0,j=0;  /*loop variables*/
  int jS=0;
  int place=0;  /*array place*/
  //int maxPlace=0; /*for tracking the smoothed wave*/
  int nIn=0;    /*number of features*/
  int waveStart=0; /*outgoing pulse bin*/
  int width=0;      
  int findStart(int,int,char *,int);
  float findQstart(int,int,char *,int,control *);
  float findSstart(int,int,char *,int,float,float *);
  float pRange=0;         /*peak range*/
  float sRange=0;         /*mean range weighted by intensity*/
  float qPeak=0,qRange=0; /*quadratic fit parameters*/
  float qIntegral=0;      /*quadratic integral*/
  float sum=0.0;     /*return energy sum*/
  float rho=0;            /*target reflectance*/
  float res=0;
  float sStart=0;   /*signal starts found from quadrativ and sum methods*/
  float salcaRange(float,float,float,float);
  float *ranges=NULL,*refls=NULL;
  float tZen=0,tAz=0;     /*angles after matrix has been applied*/
  char max=0;             /*intensity threshold and max value*/
  char satTest=0;         /*saturation test indicator*/
  char ringTest=0;        /*test for ringing*/
  char checkRinging(int,int,char *,int);
  char brEak=0;           /*break indicator*/
  void markPoint(float,int,int,float,int,char,float,char,float,int,float,int,FILE *,control *,float,float,float,double *,float,float);
  void fitFunction(char *,int,float *,float *,float *,float,int,int,float *,control *);
  void fitQuad(char *,int,float *,float *,float,control *,float *);
  void writeBinAngles(control *,int,FILE *,float,float);
  void writeBinHits(int,float *,float *,FILE *);
  void recordHit(float **,float **,int,float,float);
  void applyMatrix(float *,float *,float,float,double **);
  void leadingEdgeRange(int,int,int,control *,char *,int,int,float *,float,float);
  float sCount=0;


  res=options->res;
  rho=-1.0;           /*nonesense value*/

  qRange=qPeak=qIntegral=-1.0; /*set to nonesense*/

  for(i=0;i<numb;i++){ /*zenith loop*/
    nIn=1;                      /*number of points per waveform*/

    /*apply collocation matrix*/
    applyMatrix(&tZen,&tAz,options->zen[band][i],options->azOff[band][i]+az,options->trans);
  
    /*if needed, write the angles*/
    if(options->writeBin)writeBinAngles(options,band,opoo,tZen,tAz);


    /*find the waveform start and check for saturation*/
    /*peak*/
    waveStart=findStart(start,end,data,i*nBins);
    if((band==0)&&(waveStart<0))waveStart=0; /*as the 1545 pulse is often lost off the end*/
    /*quadratic*/
    /*if(options->doQuad)qStart=findQstart(start,end,data,i*nBins,options);*/
    /*sum*/
    sStart=findSstart(start,end,data,i*nBins,options->res,&sCount);

    /*check for ringing*/
    ringTest=checkRinging(start,end,data,i*nBins);

    if(waveStart>=0){   /*if there is a return here*/
      for(j=start;j<end;j++){     /*loop along waveform*/
        place=i*nBins+j;            /*array place*/
        if(data[place]>thresh){   /*signal above noise level*/
          satTest=0;      /*reset variables*/
          sRange=0.0;
          width=0;
          max=thresh;
          sum=0.0;

          jS=j;                   /*mark data for function fitting*/
          for(;j<end;j++){        /*step to end of feature*/
            place=i*nBins+j;        /*array place*/
            sum+=(float)data[place]-backSig;
            sRange+=(((float)(j-start)+1.5)*options->res-sStart)*((float)data[place]-backSig);  /*weighted mean range*/
            if(data[place]>max){  /*record maximum intensity and position*/
              max=data[place];      /*set max*/
              //maxPlace=place;       /*record array point for later tracking*/
              pRange=((float)(j-waveStart)+1.5)*0.15;  /*use outgoing pulse to set start*/  /*why +1.5?*/
            }                     /*max test*/

            if(data[place]==127)satTest=1;

            if(data[place]<=thresh){ /*we have left feature*/
              if(sum>0.0)sRange/=sum;
              else       sRange=pRange;

              /*fit quadratic*/
              if(options->doQuad)fitQuad(&(data[i*nBins+jS-1]),j-jS+1,&qRange,&qPeak,(float)(jS-waveStart)*res,options,&qIntegral);

              if(calibrate){
                if(((int)sum>=0)&&((int)sum<cal->numb))rho=cal->LUT[(int)sum]/salcaRange(cal->a,cal->k2,cal->k3,sRange);
                else{
                  if((int)sum<=0)rho=0.0;
                  else if((int)sum>=cal->numb)rho=1.0;
                }
                if(rho>1.0)rho=1.0;
              }else rho=sum;

              /*if saturated, use leading edge to get range*/
              if(ringTest&&(max>=satThresh))leadingEdgeRange(start,jS,j,options,data,i,nBins,&sRange,backSig,sStart);

              if(options->writeBin==0){  /*write an ASCII file*/
                markPoint(sRange,i,azInd,sum,band,max,pRange,satTest,rho,nIn,sStart,width,opoo,options,qRange,qPeak,qIntegral,options->scanCent,tZen,tAz);
              }else{            /*save up hits to write binary*/
                recordHit(&ranges,&refls,nIn-1,sRange,rho);
              }

              width=0;
              if(ringTest&&(max>=satThresh))brEak=1;   /*if saturated only take the first return to avoid ringing*/
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

    /*write out hits if any*/
    if(options->writeBin){
      writeBinHits(nIn-1,ranges,refls,opoo);
      TIDY(ranges);
      TIDY(refls);
    }
  }/*zenith loop*/
  return;
}/*pointOut*/


/*##########################################*/
/*range from leading edge*/

void leadingEdgeRange(int start,int jS,int jEnd,control *options,char *data,int i,int nBins,float *sRange,float backSig,float sStart)
{
  int j=0,jE=0;
  int place=0;
  int nAbove=0;
  float hm=0,d=0;
  float sigma=0;
  char max=0;

  /*pulse length in matres*/
  sigma=(1.1/1000000000)*2.98*100000000;

  /*does it need doing?*/
  nAbove=0;
  for(j=jS;j<jEnd;j++){
    place=i*nBins+j;
    if(data[place]>=110)nAbove++;
  }

  /*did it need fixing?*/
  max=-125;
  if(nAbove>2){
    /*find area of interest*/
    for(j=jS;j<jEnd;j++){
      place=i*nBins+j;
      if(data[place]>=satThresh){
        jE=j;
        max=data[place];
        break;
      }
    }

    /*determine range from these points*/
    hm=((float)max-(float)backSig)/2.0;
    for(;j>=jS;j--){
      place=i*nBins+j;
      /*interpolate*/
      if((float)data[place]<=hm){
        d=((float)max-hm)/(float)(max-data[place]);
        *sRange=((d*(float)(jE-j)+(float)j-(float)start)*options->res-sStart)+sigma*sqrt(2.0*log(2.0));
        break;
      }
    }
  }/*need fixing?*/


  return;
}/*leadingEdgeRange*/


/*##########################################*/
/*record hiyts for binary file*/

void recordHit(float **ranges,float **refls,int nIn,float range,float rho)
{
  ranges[0]=markFloat(nIn,ranges[0],range);
  refls[0]=markFloat(nIn,refls[0],rho);

  return;
}/*recordHit*/


/*##########################################*/
/*apply collocation matrix*/

void applyMatrix(float *tZen,float *tAz,float zen,float az,double **trans)
{
  double inVect[3],*outVect=NULL;

  /*get vector from angles*/
  inVect[0]=sin(az)*sin(zen);
  inVect[1]=cos(az)*sin(zen);
  inVect[2]=cos(zen);

  /*apply rotation*/
  outVect=rotateVect(inVect,trans);

  /*copy results*/
  *tZen=acos(outVect[2]);
  *tAz=atan2(outVect[0],outVect[1]);

  TIDY(outVect);
  return;
}/*applyMatrix*/


/*#######################################################################*/
/*rotate a vector*/

double *rotateVect(double *vect,double **matrix)
{
  int i=0,j=0;
  double *rotated=NULL;

  rotated=dalloc(3,"rotated vector",0);

  for(i=0;i<3;i++){
    rotated[i]=0.0;
    for(j=0;j<3;j++)rotated[i]+=vect[j]*matrix[i][j];
  }

  return(rotated);
}/*rotateVect*/


/*##########################################*/
/*write binary hits*/

void writeBinHits(int nIn,float *ranges,float *refls,FILE *opoo)
{
  int i=0,len=0;
  uint16_t offset=0;
  char *buff=NULL;
  uint8_t nHits=0;

  len=1+nIn*8;
  buff=challoc(len,"hit buffer",0);
  if((nIn>(int)pow(2,8))||(nIn<0)){
    fprintf(stderr,"Too many hits to be recorded. File will be corrupted %d\n",nIn);
    exit(1);
  }
  nHits=(uint8_t)nIn;

  memcpy(&(buff[0]),&nHits,1);
  offset=1;
  for(i=0;i<nIn;i++){
    memcpy(&(buff[offset]),&ranges[i],4);
    offset+=4;
    memcpy(&(buff[offset]),&refls[i],4);
    offset+=4;
  }

  if(fwrite(&(buff[0]),sizeof(char),len,opoo)!=len){
    printf("Error writing to output\n");
    exit(1);
  }

  TIDY(buff);
  return;
}/*writeBinHits*/


/*##########################################*/
/*write binary ending*/

void writeBinEnding(uint32_t *shotN,double *scanCent,FILE **opoo)
{
  int band=0,len=0;
  char *buff=NULL;

  len=4+3*8;
  buff=challoc(len,"binary ending buffer",0);

  for(band=0;band<2;band++){
    memcpy(&(buff[0]),&scanCent[0],8);
    memcpy(&(buff[8]),&scanCent[1],8);
    memcpy(&(buff[16]),&scanCent[2],8);
    memcpy(&(buff[24]),&shotN[band],4);

    /*write to binary file*/
    if(fwrite(&(buff[0]),sizeof(char),len,opoo[band])!=len){
      printf("Error writing to output\n");
      exit(1);
    }
  }

  TIDY(buff);
  return;
}/*writeBinEnding*/


/*##########################################*/
/*write binary angles*/

void writeBinAngles(control *options,int band,FILE *opoo,float tZen,float tAz)
{
  int len=6*4;
  char buff[len];
  float cZen=0,cAz=0;
  float cent=0;

  cZen=tZen*180.0/M_PI;
  cAz=tAz*180.0/M_PI;
  cent=0.0;   /*SALCA's optical centre is 0,0,0*/

  /*load variables in to buffer*/
  memcpy(&(buff[0]),&cZen,4);
  memcpy(&(buff[4]),&cAz,4);
  memcpy(&(buff[8]),&cent,4);
  memcpy(&(buff[12]),&cent,4);
  memcpy(&(buff[16]),&cent,4);
  memcpy(&(buff[20]),&options->shotN[band],4);

  /*increment counter*/
  options->shotN[band]++;

  /*write to binary file*/
  if(fwrite(&(buff[0]),sizeof(char),len,opoo)!=len){
    printf("Error writing to output\n");
    exit(1);
  }

  return;
}/*writeBinAngles*/


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

void markPoint(float sRange,int zenInd,int azInd,float sum,int band,char max,float pRange,char satTest,float rho,int nIn,float sStart,int width,FILE *opoo,control *options,float qRange,float qPeak,float qIntegral,double *scanCent,float tZen,float tAz)
{
  double x=0,y=0,z=0;  /*cartesian coordinates*/

  x=sRange*sin(tZen)*sin(tAz)+scanCent[0];
  y=sRange*sin(tZen)*cos(tAz)+scanCent[1];
  z=sRange*cos(tZen)+scanCent[2];

  fprintf(opoo,"%.3f %.3f %.3f %g %g %g %g %d %d %d %g %d %g %d %g %g %g %d\n",x,y,z,sum,tZen*180.0/M_PI,tAz*180.0/M_PI,sRange,nIn,azInd,zenInd,sStart,satTest,rho,width,qRange,qPeak,qIntegral,max);

  return;
}/*markPoint*/


/*##########################################*/
/*open two output pointers and write header*/

FILE **openOutput(char *outRoot,char writeBin)
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
    if(writeBin==0){  /*ASCII files*/
      sprintf(outNamen,"%s.band.%d.data",outRoot,band);
      if((opoo[band]=fopen(outNamen,"w"))==NULL){
        fprintf(stderr,"Error opening output file %s\n",outNamen);
        exit(1);
      }
      fprintf(opoo[band],"# 1 x, 2 y, 3 z, 4 sum, 5 zen, 6 az, 7 range, 8 number , 9 azInd, 10 zenInd, 11 wavstart, 12 satTest, 13 rho, 14 width, 15 qRange, 16 qPeak, 17 qIntegral,18 peak\n");
    }else{   /*binary files*/
      sprintf(outNamen,"%s.band.%d.bin",outRoot,band);
      if((opoo[band]=fopen(outNamen,"wb"))==NULL){
        fprintf(stderr,"Error opening output file %s\n",outNamen);
        exit(1);
      }
    }
  }
  return(opoo);
}/*openOutput*/


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

    /*check file size via bins*/
    if(((*nBins)*(*numb))>(*length)){
      fprintf(stderr,"File size mismatch\n");
      exit(1);
    }
    data=challoc(*length,"data",0);

    /*now we know how long, read the file*/
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
/*close output files*/

void closeFiles(FILE **opoo,char *outRoot,char writeBin)
{
  int band=0;

  for(band=0;band<2;band++){
    if(opoo[band]){
      fclose(opoo[band]);
      opoo[band]=NULL;
    } 
    if(writeBin==0)fprintf(stdout,"Written to %s.band.%d.data\n",outRoot,band);
    else           fprintf(stdout,"Written to %s.band.%d.bin\n",outRoot,band);
  }
  TIDY(opoo);

  return;
}/*closeFiles*/


/*##########################################*/
/*read calibration LUT*/

calibration *readCalibration(char *calFile,float filt,float filt1,float filt2,control *options)
{
  int band=0,i=0;
  calibration *cal=NULL;
  char line[200];
  char temp1[100],temp2[100],temp3[100],temp4[100];
  void setTransmission(calibration *,float);
  FILE *ipoo=NULL;

  if(!(cal=(calibration *)calloc(2,sizeof(calibration)))){
    fprintf(stderr,"error in calibration allocation.\n");
    exit(1);
  }
  cal[0].minDN=cal[1].minDN=0;
  cal[0].maxDN=cal[1].maxDN=4000;
  cal[0].numb=cal[1].numb=cal[0].maxDN-cal[0].minDN+1;  /*hard wired for now*/

  if(options->calOldStyle){
    if(filt1<0.0)setTransmission(cal,filt);
    else{
      cal[0].tran=filt1;
      cal[1].tran=filt2;
    }
  }else{
    cal[0].tran=cal[1].tran=1.0;
  }

  /*read LUT file*/
  for(band=0;band<2;band++)cal[band].LUT=falloc(cal[band].numb,"calibration LUT",band);
  if((ipoo=fopen(calFile,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",calFile);
    exit(1);
  }
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
        i=atoi(temp1);
        if((i<0)||(i>=cal[0].numb)){
          fprintf(stderr,"Stepped too far\n");
          exit(1);
        }
        cal[0].LUT[i]=atof(temp2)/(cal[0].tran*cal[0].tran*options->I0[0]);  /*scale by filter transmission*/
        cal[1].LUT[i]=atof(temp3)/(cal[1].tran*cal[1].tran*options->I0[1]);
      }
    }else{
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
        if(!strncasecmp(temp2,"a",1)){ /*read the a*/
          cal[0].a=atof(temp3);
          cal[1].a=atof(temp4);
        }else if(!strncasecmp(temp2,"k2",2)){ /*read the a*/
          cal[0].k2=atof(temp3);
          cal[1].k2=atof(temp4);
        }else if(!strncasecmp(temp2,"k3",2)){ /*read the a*/
          cal[0].k3=atof(temp3);
          cal[1].k3=atof(temp4);
        }
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
  int i=0,band=0;
  float zen=0,az=0;
  float cZen=0,cAz=0;
  void squint(float *,float *,float,float,float,float,float);

  options->zen=fFalloc(2,"zenith squin",0);
  options->azOff=fFalloc(2,"azimuth squin",0);
  for(band=0;band<2;band++){
    options->zen[band]=falloc(numb,"zenith squint",0);
    options->azOff[band]=falloc(numb,"azimuth squint",0);

    az=0.0;
    for(i=0;i<numb;i++){
      zen=((float)(options->nZen/2)-(float)i)*options->zStep+options->mSquint;
      squint(&(cZen),&(cAz),zen,az,options->zenSquint[band],options->azSquint[band],options->omega);
      options->zen[band][i]=cZen;
      options->azOff[band][i]=cAz;
    }/*zenith loop*/
  }/*band loop*/

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
  int band=0;
  float sinAz=0,sinZen=0;

  for(band=0;band<2;band++){
    sinZen=sin(options->zenSquint[band]);
    sinAz=sin(options->azSquint[band]);

    options->azSquint[band]=atan2(sinAz,sinZen);
    options->zenSquint[band]=atan2(sqrt(sinAz*sinAz+sinZen*sinZen),1.0);
  }

  return;
}/*translateSquint*/


/*##############################################*/
/*salca intensity range function*/

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
/*read the command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *options=NULL;

  /*allocate space*/
  if(!(options=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }

  /*defaults*/
  strcpy(options->inRoot,"/home/server/users/bakgrp3/nsh103/data/SALCA/raw/hedge_plot34_f1_1936/hedge_plot34_f1_1936");
  strcpy(options->outRoot,"salcaTest");
  strcpy(options->calFile,"/home/sh563/src/salcaMake/salcaLUT.dat");
  strcpy(options->transName,"none");
  strcpy(options->globName,"none");

  /*switches*/
  options->coarsen=1;   /*use at native resolution*/
  options->calibrate=0;/*output DN rather than reflectance*/
  options->writeBin=0;  /*make an ASCII file by default*/
  options->maxR=30.0;
  /*set shot counter*/
  options->shotN[0]=options->shotN[1]=0;
  /*filters*/
  options->filt=0.0;
  options->filt1=options->filt2=-1.0;

  options->azStep=0.06;   /*in degrees*/
  options->maxZen=190.0;   /*full hemispheric scan*/
  options->azSquint[0]=-2.311151*M_PI/180.0;   /*these angles found by*/
  options->azSquint[1]=-2.272465*M_PI/180.0;   /*~/src/geomSALCA/geomSALCAline.c*/
  options->zenSquint[0]=-1.199219*M_PI/180.0;
  options->zenSquint[1]=-1.181007*M_PI/180.0;

  options->azSquint[0]=options->azSquint[1];
  options->zenSquint[0]=options->zenSquint[1];

  options->omega=M_PI/4.0;  /*45 degrees*/
  options->mSquint=1.490432*M_PI/180.0;
  options->azStart=0.0;
  options->nZen=3200;
  options->nAz=666;
  options->joy=1;        /*print out status*/
  options->res=0.15;
  options->smoother=NULL;
  options->sWidth=1.1;
  options->nSmoo=0;
  options->doQuad=0;  /*don't do Jupp's quadratic fitting*/
  options->I0[0]=options->I0[1]=1.0;
  options->calOldStyle=0;
  options->scanCent[0]=options->scanCent[1]=options->scanCent[2]=0.0;

  satThresh=125;      /*ringing thresholdf*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-inRoot",7)){
        checkArguments(1,i,argc,"-inRoot");
        strcpy(options->inRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(options->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-trans",6)){
        checkArguments(1,i,argc,"-trans");
        strcpy(options->transName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-globMat",8)){
        checkArguments(1,i,argc,"-globMat");
        strcpy(options->globName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-nAz",4)){
        checkArguments(1,i,argc,"-nAz");
        options->nAz=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azStep",7)){
        checkArguments(1,i,argc,"-azStep");
        options->azStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxR",5)){
        checkArguments(1,i,argc,"-maxR");
        options->maxR=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-calibrate",10)){
        options->calibrate=1;
      }else if(!strncasecmp(argv[i],"-calFile",8)){
        checkArguments(1,i,argc,"-calFile");
        strcpy(options->calFile,argv[++i]);
        options->calibrate=1;
      }else if(!strncasecmp(argv[i],"-filt",5)){
        checkArguments(1,i,argc,"-filt");
        options->filt=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-laser",6)){
        checkArguments(2,i,argc,"-laser");
        options->I0[0]=atof(argv[++i]);
        options->I0[1]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nFilt",6)){
        checkArguments(2,i,argc,"-nFilt");
        options->filt1=atof(argv[++i]);
        options->filt2=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-coarsen",8)){
        checkArguments(1,i,argc,"-coarsen");
        options->coarsen=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        options->maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-azSquint",9)){
        checkArguments(2,i,argc,"-azSquint");
        options->azSquint[0]=atof(argv[++i])*M_PI/180.0;
        options->azSquint[1]=atof(argv[++i])*M_PI/180.0;
      }else if(!strncasecmp(argv[i],"-zenSquint",10)){
        checkArguments(2,i,argc,"-zenSquint");
        options->zenSquint[0]=atof(argv[++i])*M_PI/180.0;
        options->zenSquint[1]=atof(argv[++i])*M_PI/180.0;
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
      }else if(!strncasecmp(argv[i],"-quad",5)){
        options->doQuad=1;
      }else if(!strncasecmp(argv[i],"-writeBin",9)){
        options->writeBin=1;
      }else if(!strncasecmp(argv[i],"-satThresh",11)){
        checkArguments(1,i,argc,"-satThresh");
        satThresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-inRoot root;   binary input filename root (less the _0.bin)\n-outRoot root;  output filename root\n-nAz n;         number of azimuth steps (number of binary files)\n-azStep ang;    azimuth step in degrees\n-maxR range;    maximum recorded range, metres\n-calibrate;     calibrate to reflectance, need calibration file\n-calFile name;  calibrate to reflectance, need calibration file\n-laser p0 p1;   laser power and filter scalar\n-rangePar file;      file holding calibration parameters\n-filt n;        filter optical depth (0, 0.6, 1 or 1.6)\n-nFilt t1 t2;   specify transmissions rather than use defaults\n-coarsen n;     coarsen by a factor. CAUTION, intensities will not scale properly yet\n-maxZen zen;    maximum zenith angle, degrees. 190 by default\n-azSquint ang0 ang1;    azimuth squint angle in degrees for each beam\n-zenSquint ang0 ang1;   zenith squint angle, degrees for each beam\n-mSquint angle;     mirror squint angle in degrees\n-omega angle;       mirror angle in degrees\n-azStart angle;     azimuth start, degrees\n-background DN;     background DN value\n-joyless;           don't print out status\n-sWidth width;      smoothing width for quadratic fitting\n-quad;              do quadratic fitting\n-satThresh DN;      DN to test for saturation. 125 by defaultz\n-writeBin;          write a binary file for use in voxelTLS\n\nGeolocation\n-trans name;     collocation matrix filename\n-globName name;  geolocation matrix filename\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  return(options);
}/*readCommands*/


/*the end*/
/*##########################################*/

