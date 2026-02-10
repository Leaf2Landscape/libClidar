#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLidVoxel.h"



/*###################################*/
/*main*/

int main()
{
  uint32_t i=0;
  char namen[200];
  tlsScan *scan=NULL;
  void readPTX(char *,uint32_t,tlsScan **);

  sprintf(namen,"botanics1.ptx");

  /*initial read*/
  readPTX(namen,0,&scan);

  /*loop over beams*/
  for(i=0;i<scan->nBeams;i++){
    readPTX(namen,i,&scan);

    if(scan->beam[i-scan->pOffset].nHits==0)fprintf(stdout,"%f %f %d %u\n",scan->beam[i-scan->pOffset].zen,scan->beam[i-scan->pOffset].az,scan->beam[i-scan->pOffset].nHits,i);
  }

  return(0);
}/*main*/


/*###################################*/
/*read a PTX file*/

void readPTX(char *namen,uint32_t place,tlsScan **scan)
{
  int j=0;
  uint32_t i=0;
  uint32_t xInd=0,yInd=0;
  uint32_t contN=0;
  uint64_t fStart=0;
  double x=0,y=0,z=0;
  float zen=0,az=0,diff=0;
  float lastZen=-400.0,lastAz=-400.0;
  char line[100];
  char temp1[25],temp2[25];
  char temp3[25],temp4[25];
  void translateLeica(double *,double *,double *,float **);
  /*keep a list of done beams, only in here*/
  static uint32_t nRow;
  static double res;

  /*read header if first time*/
  if((*scan)==NULL){
    res=0.0;
    contN=0;
    /*allocate space*/
    if(!((*scan)=(tlsScan *)calloc(1,sizeof(tlsScan)))){
      fprintf(stderr,"error in tls structure allocation.\n");
      exit(1);
    }
    (*scan)->matrix=fFalloc(4,"translation matrix",0);
    for(j=0;j<4;j++)(*scan)->matrix[j]=falloc(4,"translation matrix",j+1);

    /*open file*/
    if(((*scan)->ipoo=fopen(namen,"r"))==NULL){
      fprintf(stderr,"Error opening input file %s\n",namen);
      exit(1);
    }

    /*read header and determine file length*/
    i=0;
    while(fgets(line,100,(*scan)->ipoo)!=NULL){
      /*read the header*/
      if(i<10){
        /*if(i==0)nCol=atoi(line);*/
        if(i==1)nRow=atoi(line);
        else{
          if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){  /*translation matrix*/
            (*scan)->matrix[0][j]=atof(temp1);
            (*scan)->matrix[1][j]=atof(temp2);
            (*scan)->matrix[2][j]=atof(temp3);
            (*scan)->matrix[3][j]=atof(temp4);
            j++;
          }else if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            if(i==2){
              (*scan)->xOff=atof(temp1);
              (*scan)->yOff=atof(temp2);
              (*scan)->zOff=atof(temp3);
              j=0;
            }
          }
        }
      }else if(i==10)fStart=ftell((*scan)->ipoo);
      else{  /*work out resolution*/
        if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
          x=atof(temp1);
          y=atof(temp2);
          z=atof(temp3);
          if((fabs(x)+fabs(y)+fabs(z))>0.0001){
            zen=atan2(sqrt(x*x+y*y),z)*180.0/M_PI;
            az=atan2(x,y)*180.0/M_PI;
            if(lastZen>=-180.0){
              diff=sqrt(pow(lastZen-zen,2)+pow(lastAz-az,2));
              if(diff<1.0){  /*avoid ends of scan lines*/
                res+=(double)diff;
                contN++;
              }
            }
            lastZen=zen;
            lastAz=az;
  
          }else lastZen=lastAz=-400.0;
        }
      }
      i++;
    }
    res/=(double)contN;
    (*scan)->nBeams=i-10; /*ignore header*/
    (*scan)->maxRead=500000000/((uint64_t)sizeof(tlsBeam)+8);
    (*scan)->nRead=0;
    (*scan)->totRead=0;
    (*scan)->pOffset=0;
    fprintf(stdout,"Scan contains %u beams\n",(*scan)->nBeams);

    /*allocate space*/
    if(!((*scan)->beam=(tlsBeam *)calloc((long)(*scan)->maxRead,sizeof(tlsBeam)))){
      fprintf(stderr,"error beam allocation.\n");
      exit(1);
    }
    /*rewind to start of data blocks*/
    if(fseek((*scan)->ipoo,(long)fStart,SEEK_SET)){
      fprintf(stderr,"fseek error to start\n");
      exit(1);
    }
  }else if((place-(*scan)->pOffset)<(*scan)->nRead){  /*if still within block, return*/
    return;
  }/*header reading*/


  /*update offset*/
  if(place>0)(*scan)->pOffset+=(*scan)->nRead;

  /*clear out old memory*/
  for(i=0;i<(*scan)->maxRead;i++){
    TIDY((*scan)->beam[i].refl);
    TIDY((*scan)->beam[i].r);
  }

  /*read block of data*/
  i=0;
  while(fgets(line,100,(*scan)->ipoo)!=NULL){
    if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
      x=atof(temp1);
      y=atof(temp2);
      z=atof(temp3);

      (*scan)->beam[i].x=0.0;  /*origin is always 0 for a Leica*/
      (*scan)->beam[i].y=0.0;
      (*scan)->beam[i].z=0.0;
      (*scan)->beam[i].shotN=i;

      /*contains hits?*/
      if((fabs(x)+fabs(y)+fabs(z))>0.0){
        /*get angles*/
        translateLeica(&x,&y,&z,(*scan)->matrix);
        (*scan)->beam[i].zen=atan2(sqrt(x*x+y*y),z)*180.0/M_PI;
        (*scan)->beam[i].az=atan2(x,y)*180.0/M_PI;
        /*mark hit*/
        (*scan)->beam[i].nHits=1;
        (*scan)->beam[i].r=falloc(1,"range",i+1);
        (*scan)->beam[i].refl=falloc(1,"refl",i+1);
        (*scan)->beam[i].r[0]=sqrt(x*x+y*y+z*z);
        (*scan)->beam[i].refl[0]=atof(temp4);

        /*mark which pixels have been done for gap tracking*/
        xInd=(int)(((*scan)->beam[i].az+180.0)/(float)res+0.5);
        yInd=(int)((*scan)->beam[i].zen/(float)res+0.5);
      }else{  /*otherwise gap. NOTE that this should be diaganol across x*/
        if(yInd>0)yInd--;
        else{
          yInd=nRow;
          xInd++;
        }
        (*scan)->beam[i].zen=(float)res*(float)yInd;
        (*scan)->beam[i].az=(float)res*(float)xInd-180.0;
        (*scan)->beam[i].nHits=0;
        (*scan)->beam[i].shotN=i;
        (*scan)->beam[i].r=NULL;
        (*scan)->beam[i].refl=NULL;
      }

      i++;
      /*if we have read the whole chunk, break*/
      if(i>=(*scan)->maxRead){
        (*scan)->totRead=ftell((*scan)->ipoo);
        break;
      }
    }
  }/*data reading loop*/
  (*scan)->nRead=i;

  /*if this is the last call, close file*/
  if(((*scan)->pOffset+(*scan)->nRead)>=(*scan)->nBeams){
    if((*scan)->ipoo){
      fclose((*scan)->ipoo);
      (*scan)->ipoo=NULL;
    }
    TTIDY((void **)(*scan)->matrix,4);
  }/*file closing check*/

  return;
}/*readPTX*/


/*###################################*/
/*translate leica coord*/

void translateLeica(double *x,double *y,double *z,float **matrix)
{
  double tempX=0,tempY=0,tempZ=0;

  tempX=*x;
  tempY=*y;
  tempZ=*z;

  *x=tempX*matrix[0][0]+tempY*matrix[1][0]+tempZ*matrix[2][0];
  *y=tempX*matrix[0][1]+tempY*matrix[1][1]+tempZ*matrix[2][1];
  *z=tempX*matrix[0][2]+tempY*matrix[1][2]+tempZ*matrix[2][2];

  return;
}/*translateLeica*/


/*the end*/
/*###################################*/

