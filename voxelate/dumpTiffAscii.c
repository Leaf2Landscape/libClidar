#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdint.h"
#include "math.h"
#include "tools.h"
#include "tools.c"
#include "tiffRead.h"
#include "libLasRead.h"

/*####################*/
/*# Dumps a geotiff  #*/
/*# to an ascii file #*/
/*####################*/

#define TOL 0.000001

/*#####################################################################################*/
/*control structure*/

typedef struct{
  int nFiles;          /*number of input files*/
  char **inList;       /*input filename list*/
  char outNamen[200];  /*output filename*/
  int maxNx;           /*maximum number of x pixels*/
  int maxNy;           /*maximum number of y pixels*/

  /*output image*/
  double minX;
  double maxX;
  double minY;
  double maxY;
  int nX;
  int nY;
  double res;
  float *image;
  float nodata;

  /*denosing*/
  float maxVal;    /*maximum value to allow*/
  float filtAbove; /*value to filter above*/
  int median;      /*median filter above maxVal*/
}control;


/*#####################################################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  void findFileBounds(char **,int,double *,double *,double *,double *,double *);
  void readData(control *);
  void processData(control *);
  void writeResults(control *);


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*determine output file bounds*/
  findFileBounds(dimage->inList,dimage->nFiles,&dimage->minX,&dimage->maxX,&dimage->minY,&dimage->maxY,&dimage->res);

  /*fill up data*/
  readData(dimage);

  /*process if needed*/
  processData(dimage);

  /*output data*/
  writeResults(dimage);

  /*tidy up*/
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage->image);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*#####################################################################################*/
/*priocess data*/

void processData(control *dimage)
{
  int i=0,j=0,ii=0,jj=0;
  int place=0,iPlace=0;
  int nIn=0;
  float *processed=NULL;
  float *temp=NULL;
  float singleMedian(float *,int);

  /*allocate space*/
  processed=falloc((uint64_t)dimage->nX*(uint64_t)dimage->nY,"processed image",0);

  if(dimage->median>0){
    temp=falloc((uint64_t)dimage->median*(uint64_t)dimage->median,"temp median",0);
    for(i=0;i<dimage->nX;i++){
      for(j=0;j<dimage->nY;j++){
        place=j*dimage->nX+i;
        if((dimage->image[place]>dimage->filtAbove)||(dimage->image[place]>dimage->maxVal)||(fabs(dimage->image[place]-dimage->nodata)<TOL)){  /*then filter*/
          nIn=0;
          for(ii=i-dimage->median/2;ii<=i+dimage->median/2;ii++){
            if((ii<0)||(ii>=dimage->nX))continue;  /*if exceed bounds, skip*/
            for(jj=j-dimage->median/2;jj<=j+dimage->median/2;jj++){
              if((jj<0)||(jj>=dimage->nY))continue;  /*if exceed bounds, skip*/
              iPlace=jj*dimage->nX+ii;
              if((dimage->image[iPlace]>(dimage->nodata+TOL))&&(dimage->image[iPlace]<dimage->filtAbove)){
                temp[nIn]=dimage->image[iPlace];
                nIn++;
              }
            }/*adjacent y loop*/
          }/*adjacent x loop*/

          /*median filter*/
          if(nIn>0)processed[place]=singleMedian(temp,nIn);
          else     processed[place]=dimage->nodata;
          if(processed[place]>dimage->maxVal)processed[place]=dimage->nodata;
        }else processed[place]=dimage->image[place];  /*otherwise copy over*/
      }/*y loop*/
    }/*x loop*/
    TIDY(temp);
  }/*median filter if needed*/

  /*copy over data*/
  TIDY(dimage->image);
  dimage->image=processed;
  processed=NULL;
  return;
}/*processData*/


/*#####################################################################################*/
/*write out results*/

void writeResults(control *dimage)
{
  int i=0,j=0,place=0;
  int x=0,y=0;
  int nX=0,nY=0;
  int xOff=0,yOff=0;
  int globI=0,globJ=0;
  char namen[200],hasData=0;
  FILE *opoo=NULL;


  if(((dimage->maxNx<0)&&(dimage->maxNy<0))||((dimage->nX<=dimage->maxNx)&&(dimage->nY<=dimage->maxNy))){  /*then a single tile*/
    if((opoo=fopen(dimage->outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
      exit(1);
    }

    /*write header*/
    fprintf(opoo,"ncols %d\n",dimage->nX);
    fprintf(opoo,"nrows %d\n",dimage->nY);
    fprintf(opoo,"xllcorner %.2f\n",dimage->minX);
    fprintf(opoo,"yllcorner %.2f\n",dimage->minY);
    fprintf(opoo,"cellsize %f\n",dimage->res);
    fprintf(opoo,"NODATA_value %.1f\n",dimage->nodata);

    for(j=0;j<dimage->nY;j++){
      for(i=0;i<dimage->nX;i++){
        place=j*dimage->nX+i;
        fprintf(opoo,"%f ",dimage->image[place]);
      }
      fprintf(opoo,"\n");
    }/*y loop*/

    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",dimage->outNamen);
  }else{                       /*then multiple tiles*/
    nX=(int)(dimage->nX/dimage->maxNx)+1;
    nY=(int)(dimage->nX/dimage->maxNy)+1;

    /*loop over multiple tiles*/
    for(x=0;x<nX;x++){  /*x tile loop*/
      if(dimage->maxNx>=0)xOff=x*dimage->maxNx;
      else                xOff=0;
      for(y=0;y<nY;y++){  /*y tile loop*/
        if(dimage->maxNy>=0)yOff=y*dimage->maxNy;
        else                yOff=0;


        /*check contains data*/
        hasData=0;
        for(j=0;j<dimage->maxNy;j++){
          globJ=j+y*dimage->maxNy;
          if(globJ>=dimage->nY)continue;
          for(i=0;i<dimage->maxNy;i++){
            globI=i+x*dimage->maxNx;
            if(globI>=dimage->nX)continue;
            place=globJ*dimage->nX+globI;
            if(fabs(dimage->image[place]-dimage->nodata)>TOL){
              hasData=1;
              break;
            }
          }
          if(hasData)break;
        }/*y loop*/


        if(hasData){
          /*open file*/
          sprintf(namen,"%s.%d.%d.asc",dimage->outNamen,x,y);
          if((opoo=fopen(namen,"w"))==NULL){
            fprintf(stderr,"Error opening output file %s\n",namen);
            exit(1);
          }

          /*write header*/
          if(x<(nX-1))fprintf(opoo,"ncols %d\n",dimage->maxNx);
          else        fprintf(opoo,"ncols %d\n",dimage->nX%dimage->maxNx);
          if(y<(nY-1))fprintf(opoo,"nrows %d\n",dimage->maxNy);
          else        fprintf(opoo,"nrows %d\n",dimage->nY%dimage->maxNy);
          fprintf(opoo,"xllcorner %.2f\n",dimage->minX+(double)(x*dimage->maxNx)*dimage->res);
          if(y<(nY-1))fprintf(opoo,"yllcorner %.2f\n",dimage->maxY-(double)((y+1)*dimage->maxNy)*dimage->res);
          else        fprintf(opoo,"yllcorner %.2f\n",dimage->minX);
          fprintf(opoo,"cellsize %f\n",dimage->res);
          fprintf(opoo,"NODATA_value %.1f\n",dimage->nodata);

          for(j=0;j<dimage->maxNy;j++){
            globJ=j+y*dimage->maxNy;
            if(globJ>=dimage->nY)continue;
            for(i=0;i<dimage->maxNy;i++){
              globI=i+x*dimage->maxNx;
              if(globI>=dimage->nX)continue;
              place=globJ*dimage->nX+globI;
              fprintf(opoo,"%f ",dimage->image[place]);
            }
            fprintf(opoo,"\n");
          }/*y loop*/

          if(opoo){
            fclose(opoo);
            opoo=NULL;
          }
          fprintf(stdout,"Written to %s\n",namen);
        }/*has data check*/
      }/*y tile loop*/
    }/*x tile loop*/
  }/*number of tiles check*/
  return;
}/*writeResults*/


/*#####################################################################################*/
/*read data input output array*/

void readData(control *dimage)
{
  int i=0,j=0,numb=0;
  int xBin=0,yBin=0;
  int place=0,iPlace=0;
  double x=0,y=0;
  double tMinX=0,tMaxY=0;
  geot *data=NULL;

  /*work out boundaries*/
  dimage->nX=(int)((dimage->maxX-dimage->minX)/dimage->res);
  dimage->nY=(int)((dimage->maxY-dimage->minY)/dimage->res);
  dimage->image=falloc((uint64_t)dimage->nX*(uint64_t)dimage->nY,"image",0);


  /*set all to nodata*/
  for(i=dimage->nX*dimage->nY-1;i>=0;i--)dimage->image[i]=dimage->nodata;


  for(numb=0;numb<dimage->nFiles;numb++){  /*file loop*/
    /*allocate space*/
    if(!(data=(geot *)calloc(1,sizeof(geot)))){
      fprintf(stderr,"error control allocation.\n");
      exit(1);
    }

    /*read data*/
    readGeotiff(data,dimage->inList[numb],1);
    tMinX=data->tiepoints[3]-data->tiepoints[0]*data->scale[0];
    tMaxY=data->tiepoints[4]+data->tiepoints[1]*data->scale[1];


    /*store in output array*/
    for(i=0;i<data->nX;i++){
      x=(double)i*data->scale[0]+tMinX;
      xBin=(int)((x-dimage->minX)/dimage->res);
      for(j=0;j<data->nY;j++){
        y=tMaxY-(double)j*data->scale[1];
        yBin=(int)((y-dimage->minY)/dimage->res);
        place=(dimage->nY-yBin)*dimage->nX+xBin;
        iPlace=j*data->nX+i;

        if((xBin<0)||(xBin>=dimage->nX)||(yBin<0)||(yBin>dimage->nY)){
          fprintf(stderr,"bin error %d of %d\n",yBin,dimage->nY);
          exit(1);
        }

        if(data->image)dimage->image[place]=(float)data->image[iPlace]*(float)data->scale[2];
        else if(data->fImage)dimage->image[place]=data->fImage[iPlace]*(float)data->scale[2];
        else if(data->dImage)dimage->image[place]=(float)data->dImage[iPlace]*(float)data->scale[2];
        else{
          fprintf(stderr,"No data in this file\n");
          exit(1);
        }
      }/*y loop*/
    }/*x loop*/

    /*tidy as we go along*/
    data=tidyTiff(data);
  }/*loop over files*/

  return;
}/*readData*/


/*#####################################################################################*/
/*find file bounds*/

void findFileBounds(char **inList,int nFiles,double *minX,double *maxX,double *minY,double *maxY,double *res)
{
  int i=0;
  double tMinX=0,tMinY=0;
  double tMaxX=0,tMaxY=0;
  geot *data=NULL;

  /*work out file bounds*/
  *minX=*minY=1000000000000.0;
  *maxX=*maxY=-1000000000000.0;
  for(i=0;i<nFiles;i++){  /*file loop*/
    /*allocate space*/
    if(!(data=(geot *)calloc(1,sizeof(geot)))){
      fprintf(stderr,"error control allocation.\n");
      exit(1);
    } 

    /*read data*/
    readGeotiff(data,inList[i],0);

    /*check resolution*/
    if(*res<=0.0)*res=data->scale[0];
    if((fabs(*res-data->scale[0])>=TOL)||((fabs(*res-data->scale[1])>=TOL))){
      fprintf(stderr,"Resolution mismatch - %f %f %f\n",*res,data->scale[0],data->scale[1]);
      exit(1);
    }

    /*work out bounds*/
    tMinX=data->tiepoints[3]-data->tiepoints[0]*data->scale[0];
    tMaxX=tMinX+(double)data->nX*data->scale[0];
    tMaxY=data->tiepoints[4]+data->tiepoints[1]*data->scale[1];
    tMinY=tMaxY-(double)data->nY*data->scale[1];

    if(tMinX<(*minX))*minX=tMinX;
    if(tMaxX>(*maxX))*maxX=tMaxX;
    if(tMinY<(*minY))*minY=tMinY;
    if(tMaxY>(*maxY))*maxY=tMaxY;

    /*tidy as we go along*/
    data=tidyTiff(data);
  }/*loop over files*/
  return;
}/*findFileBounds*/


/*#####################################################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->res=-1.0;   /*set blank*/
  dimage->nodata=-9999.0;
  dimage->maxNx=-1;   /*outpiut all pixels*/
  dimage->maxNy=-1;   /*outpiut all pixels*/

  /*denoising*/
  dimage->maxVal=1000000000.0;     /*maximum value to allow*/
  dimage->filtAbove=1000000000.0;  /*filter all above*/
  dimage->median=0;                /*median filter above maxVal*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->nFiles=1;
        dimage->inList=chChalloc(dimage->nFiles,"input name list",0);
        dimage->inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxPix",7)){
        checkArguments(1,i,argc,"-maxPix");
        dimage->maxNx=dimage->maxNy=sqrt(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-maxXpix",8)){
        checkArguments(1,i,argc,"-maxXpix");
        dimage->maxNx=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxYpix",8)){
        checkArguments(1,i,argc,"-maxYpix");
        dimage->maxNy=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxVal",7)){
        checkArguments(1,i,argc,"-maxVal");
        dimage->maxVal=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-filtAbove",10)){
        checkArguments(1,i,argc,"-filtAbove");
        dimage->filtAbove=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-median",7)){
        checkArguments(1,i,argc,"-median");
        dimage->median=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to convert geotiffs into ascii\n#####\n\n-input name;     lasfile input filename\n-inList list;    input file list for multiple files\n-output name;    output filename\n-maxPix n;       maximum number of pixels per output file\n-maxXpix n;      maximum number of x pixels per output file\n-maxYpix n;      maximum number of y pixels per output file\n-maxVal max;     maximum value to allow\n-filtAbove max;  median filter values above\n-median n;       remove values above max with n*n median filter\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/


  return(dimage);
}/*readCommands*/

/*end*/
/*#####################################################################################*/

