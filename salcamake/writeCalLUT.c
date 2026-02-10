#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"


/*##########################*/
/*# Uses Newton Raphson to #*/
/*# solve SALCA calibration#*/
/*# in terms of DN         #*/
/*##########################*/


typedef struct{
  float m;
  float c;
  float k;
  float a;
  float k2;
  float k3;
  float *LUT;
  float minDN;
  float maxDN;
  int numb;
}calibration;

int main(int argc,char **argv)
{
  int i=0;
  calibration *cal=NULL;
  calibration *setupCalibration(char *);
  void tidyCal(calibration *);
  void writeLUT(calibration *,char *);
  char calFile[200];
  char outNamen[200];

  strcpy(outNamen,"salcaLUT.dat");
  /*strcpy(calFile,"/mnt/geodesy38/nsh103/SALCA/calibration/salcaCal.dat");*/
  strcpy(calFile,"/mnt/urban-bess/shancock_work/data/teast/salca_align/cal_targ/cockelCal.dat");

  /*read command line*/
  for(i=1;i<argc;i++){
    if(*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(calFile,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;   salca calibration parameter file\n-outut name;       output filename\n\n");
        exit(0);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry echiMake -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*read the command line*/

  cal=setupCalibration(calFile); /*set up calibration bits*/

  /*write out data*/
  writeLUT(cal,outNamen);

  tidyCal(cal);
  return(0);
}/*main*/


/*##########################################*/
/*write output*/

void writeLUT(calibration *cal,char *outNamen)
{
  int band=0,i=0;
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }
  fprintf(opoo,"# 1 DN, 2 rho 1545nm, 3 rho 1064nm\n");
  fprintf(opoo,"# a %f %f\n",cal[0].a,cal[1].a);
  fprintf(opoo,"# k2 %f %f\n",cal[0].k2,cal[1].k2);
  fprintf(opoo,"# k3 %f %f\n",cal[0].k3,cal[1].k3);
  for(i=0;i<cal[0].numb;i++){
    fprintf(opoo,"%d",i);
    for(band=0;band<=1;band++)fprintf(opoo," %f",cal[band].LUT[i]);
    fprintf(opoo,"\n");
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);

  return;
}/*writeLUT*/


/*##########################################*/
/*fill up calibration arrays*/

void calibrateArray(calibration *cal)
{
  int i=0;
  int *contN=NULL;
  float rho=0;
  float DN=0,res=0;
  float lastDN=0;
  float expo=0;
  char increaseRes=0;

  cal->minDN=0;
  cal->maxDN=4000;
  cal->numb=cal->maxDN-cal->minDN+1;
  res=0.1;

  cal->LUT=falloc(cal->numb,"calibration",0);
  contN=ialloc(cal->numb,"contribution",0);

  do{  /*iterate until resolution is just high enough*/
    increaseRes=0;
    for(i=0;i<cal->numb;i++){
     cal->LUT[i]=0.0;
     contN[i]=0;
    }

    rho=lastDN=0.0;
    do{    /*loop over reflectances*/
      if((cal->k*rho)<700.0)expo=exp(-1.0*cal->k*rho);
      else                  expo=0.0;
      DN=(cal->m*rho+cal->c)*(1.0-expo);
      if((DN-lastDN)>1.0){
        res/=10.0;
        increaseRes=1;
        break;
      }else{
        lastDN=DN;
        if(DN>=cal->numb)break;
        cal->LUT[(int)DN]+=rho;
        contN[(int)DN]++;
        rho+=res;
      }
    }while((int)DN<cal->numb);  /*reflectance loop*/
  }while(increaseRes);  /*resolution refining*/

  for(i=0;i<cal->numb;i++)if(contN[i]>0)cal->LUT[i]/=(float)contN[i];

  TIDY(contN);
  return;
}/*calibrateArray*/


/*##########################################*/
/*calibrate to reflectance*/

calibration *setupCalibration(char *calFile)
{
  int band=0;
  void readCalFile(char *,calibration *);
  void calibrateArray(calibration *);
  calibration *cal=NULL;

  if(!(cal=(calibration *)calloc(2,sizeof(calibration)))){
    fprintf(stderr,"error in calibration structure.\n");
    exit(1);
  }

  readCalFile(calFile,cal);

  for(band=0;band<2;band++){
    calibrateArray(&(cal[band]));
  }/*band loop*/
  return(cal);
}/*setupCalibration*/


/*##########################################*/
/*read calibration file*/

void readCalFile(char *calFile,calibration *cal)
{
  char line[200];
  char temp1[100],temp2[100],temp3[100];
  FILE *ipoo=NULL;

  /*read file*/
  if((ipoo=fopen(calFile,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",calFile);
    exit(1);
  }
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
        if(!strncasecmp(temp1,"m",1)){
          cal[0].m=atof(temp2);
          cal[1].m=atof(temp3);
        }else if(!strncasecmp(temp1,"c",1)){
          cal[0].c=atof(temp2);
          cal[1].c=atof(temp3);
        }else if(!strncasecmp(temp1,"k2",2)){
          cal[0].k2=atof(temp2);
          cal[1].k2=atof(temp3);
        }else if(!strncasecmp(temp1,"k3",2)){
          cal[0].k3=atof(temp2);
          cal[1].k3=atof(temp3);
        }else if(!strncasecmp(temp1,"k",1)){
          cal[0].k=atof(temp2);
          cal[1].k=atof(temp3);
        }else if(!strncasecmp(temp1,"a",1)){
          cal[0].a=atof(temp2);
          cal[1].a=atof(temp3);
        }
      }
    }
  }/*file reading loop*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*readCalFile*/


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


/* the end */
/*##############################################################################*/

