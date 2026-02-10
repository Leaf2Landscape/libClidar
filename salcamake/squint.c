#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"


#define TIDY(arr) if((arr)){free((arr));(arr)=NULL;}  /*free an array*/
float *falloc(int,char *,int);

int main()
{
  float zen=0,az=0;
  float poynt[2];
  void squint(float *,float *,float,float,float,float,float);
  float zenM=0,azM=0; /*motor angles*/
  float zE=0,aE=0;    /*squint angles*/
  float sZen=0,eZen=0,res=0;
  float omega=0;

  aE=0.0*M_PI/180.0;
  zE=0.0001*M_PI/180.0;

  azM=0.0;
  sZen=-100.0*M_PI/180.0;
  eZen=100.0*M_PI/180.0;
  res=(eZen-sZen)/100.0;
  omega=45*M_PI/180.0;

  azM=0.0;
  while(azM<=0.0){//M_PI){
    zenM=sZen;
    while(zenM<=eZen){
      squint(&(poynt[0]),&(poynt[1]),zenM,azM,zE,aE,omega);
      zen=poynt[0];
      az=poynt[1]+azM;
      if(az>(2.0*M_PI))az-=2.0*M_PI;
  
      //fprintf(stdout,"%f %f %f %f\n",zenM*180.0/M_PI,azM*180.0/M_PI,zen*180.0/M_PI,az*180.0/M_PI);

      zenM+=res;
    }/*zenith loop*/
    azM+=0.001047197;
  }

  return(0);
}/*main*/


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
  float slope=0;
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


  slope=atan2(sqrt(rX*rX+rY*rY),fabs(rZ));
  if(fabs(slope)<(M_PI/4.0)){  /*rotate about z axis*/
    thetaX=-1.0*atan2(sqrt(rX*rX+rY*rY),rZ);
    thetaZ=-1.0*atan2(rX,rY);
    rotateZ(vect,thetaZ);
    rotateX(vect,thetaX);
    rotateZ(vect,-2.0*inc);
    rotateX(vect,-1.0*thetaX);
    rotateZ(vect,-1.0*thetaZ);
  }else{                        /*rotate about y axis*/
    thetaZ=-1.0*atan2(rX,rY);
    thetaX=atan2(rZ,sqrt(rX*rX+rY*rY));
    rotateZ(vect,thetaZ);
    rotateX(vect,thetaX);
    rotateY(vect,-2.0*inc);
    rotateX(vect,-1.0*thetaX);
    rotateZ(vect,-1.0*thetaZ);
  }

  *cZen=atan2(sqrt(vect[0]*vect[0]+vect[1]*vect[1]),vect[2]);
  if(vect[1]!=0.0)*cAz=atan2(vect[0],vect[1])+aM;
  else            *cAz=aM;

  //if(zM<0.0){  /*to keep consistent with SALCA input*/
  //  *cZen*=-1.0;
  //  *cAz+=M_PI;
  //}  /*to keep consistent with SALCA input*/

  //if(*cAz>M_PI){
  //  (*cAz)-=M_PI;
  //  (*cZen)*=-1.0;
 // }

  fprintf(stdout,"zen %f az %f inc %f thetaZ %f thetaX %f slope %f r %f %f %f vect %f %f %f %f\n",(*cZen)*180.0/M_PI,(*cAz)*180.0/M_PI,inc*180.0/M_PI,thetaZ*180.0/M_PI,thetaX*180.0/M_PI,slope*180.0/M_PI,rX,rY,rZ,vect[0],vect[1],vect[2],atan2(vect[0],vect[1]));

  /*test the r vector*/
  /*vect[0]=rX;
  vect[1]=rY;
  vect[2]=rZ;
  rotateZ(vect,thetaZ);
  rotateX(vect,thetaX);
  fprintf(stdout,"r %f %f %f\n",vect[0],vect[1],vect[2]);
  rotateX(vect,-1.0*thetaX);
  rotateZ(vect,-1.0*thetaZ);
  fprintf(stdout,"rr %f %f %f before %f %f %f\n",vect[0],vect[1],vect[2],rX,rY,rZ);*/

  TIDY(vect);
  return;
}/*squint*/


/*########################################################################*/
/*allocate a float array, checking for errors*/

float *falloc(int length,char *namen,int n)
{
  float *jimlad=NULL;
  if(!(jimlad=(float *)calloc(length,sizeof(float)))){
    fprintf(stderr,"error in %s array allocation %d\n",namen,n);
    fprintf(stderr,"allocating %d\n",length);
    exit(1);
  }
  return(jimlad);
}/*falloc*/


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


/* the end */
/*#####################################################################*/

