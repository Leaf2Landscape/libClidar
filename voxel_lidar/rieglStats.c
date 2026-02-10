#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "lasHead.h"
#include "lasProcess.h"


/*#####################*/
/*# Reads binary dump #*/
/*# of Riegl rxp and  #*/
/*# makes a histogram #*/
/*# S Hancock 03 2015 #*/
/*#####################*/


/*####################################*/
/*control structure*/

typedef struct{
  char tlsNamen[200];   /*list of TLS files*/
  char output[200];     /*output filename*/
  char waveOut;         /*output waveforms*/

  /*TLS wave bits*/
  int nTLSwaves;        /*number of TLS waves per beam*/
  int binOff;           /*number of prebins to look for geolocation*/
  float appRefl1064;    /*apparent reflectance at 1064nm*/
  double beamRad;       /*beam radius*/
  uint64_t counter;     /*output number of waves output*/
  char pointOut;        /*output individual points*/
  float tlsTanDiv;      /*tan of half beam divergence*/
  float tlsRad;         /*beam start diameter*/

  /*voxel area*/
  double bounds[6];     /*area of interest*/
  float vRes;           /*voxel resolution*/

  /*denoising parameters*/
  denPar denoise;      /*structure of denoising variables*/
}control;


/*####################################*/
/*important las point data*/

typedef struct{
  int32_t x;           /*index to calculate coordinate*/
  int32_t y;           /*index to calculate coordinate*/
  int32_t z;           /*index to calculate coordinate*/
  uint16_t refl;       /*point intensity*/
  b_field field;        /*bit field. Defined in lasHead.h*/
  unsigned char scanAng;  /*scan angle*/
  uint16_t psID;          /*point source ID, used for Leica's AGC*/
  unsigned char packetDes;/*waveform packed description*/
  uint64_t waveMap;       /*pointer to waveform in file*/
  uint32_t waveLen;       /*length of waveform in bins*/
  float time;             /*time in picoseconds of this wave*/
  float grad[3];          /*waveform gradient: Vector*/
  unsigned char *wave;    /*waveform*/
}lasPoint;


/*####################################*/
/*Riegl beam, polar coords*/

typedef struct{
  double zen;     /*zenith*/
  double az;      /*azimuth*/
  float x;        /*beam origin*/
  float y;        /*beam origin*/
  float z;        /*beam origin*/
  uint32_t shotN; /*beam number*/  /*is this necessary? Could I save space without it*/
  uint8_t nHits;  /*number of hits of this beam*/
  float *r;       /*range*/
  float *refl;    /*reflectance*/
}rBeam;


/*####################################*/
/*Riegl point, Cartesian coords*/

typedef struct{
  float x;        /*point coordinate*/
  float y;        /*point coordinate*/
  float z;        /*point coordinate*/
  uint8_t hitN;   /*hit number*/
  uint8_t nHits;  /*number of hits of this beam*/
  float r;        /*range*/
  float refl;     /*reflectance*/
}rPoint;


/*####################################*/
/*Riegl scan file*/

typedef struct{
  rBeam *beam;       /*pointer to beam structures*/
  rPoint *point;     /*pointer to point structures*/
  double xOff;       /*offset to allow coords to be floats*/
  double yOff;       /*offset to allow coords to be floats*/
  double zOff;       /*offset to allow coords to be floats*/
  uint32_t nBeams;   /*number of beams in this scan*/
  uint32_t nPoints;  /*number of points in this scan*/
}rScan;


/*####################################*/
/*voxel gap structure*/

typedef struct{
  uint32_t **hits;   /*hits per scan location*/
  uint32_t **miss;   /*misses per scan location*/
  int nVox;          /*total number of voxels*/
  int nX;
  int nY;
  int nZ;
  float res;
  double bounds[6];

  /*TLS map*/
  int **mapFile;        /*file per voxel*/
  uint32_t **mapPoint;  /*point per voxel*/
  int *nIn;             /*number of points per voxel*/
}voxStruct;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  listLas *tlsList=NULL;
  rScan *tls=NULL;
  void readTLS(listLas *,control *,voxStruct *);
  rScan *tidyTLS(rScan *,int);
  voxStruct *vox=NULL;
  voxStruct *voxAllocate(int,control *);
  voxStruct *tidyVox(voxStruct *,int);
  void intersectVoxels(rScan **,listLas *,listLas *,control *);
  void readPulse(denPar *);


  /*read commands*/
  dimage=readCommands(argc,argv);

  /*read TLS input list*/
  tlsList=readLasList(dimage->tlsNamen);

  /*read TLS data*/
  readTLS(tlsList,dimage,vox);
 
  /*tidy up arrays*/
  tls=tidyTLS(tls,tlsList->nFiles);
  vox=tidyVox(vox,tlsList->nFiles);
  tidyListLas(tlsList);
  if(dimage){
    TTIDY((void **)dimage->denoise.pulse,2);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################*/
/*read TLS within plot*/

void readTLS(listLas *tlsList,control *dimage,voxStruct *vox)
{
  int i=0;
  uint32_t j=0;
  rScan *tls=NULL,*tempTLS=NULL;
  rScan *readScan(char *);
  rScan *tidyTLS(rScan *,int);
  uint32_t *histogram=NULL,total=0;
  int maxHits=0;
  FILE *opoo=NULL;

  maxHits=20;
  if(!(histogram=(uint32_t *)calloc(maxHits,sizeof(uint32_t)))){
    fprintf(stderr,"error in histogram allocation.\n");
    exit(1);
  }
  for(i=0;i<maxHits;i++)histogram[i]=0;
  total=0;


  if(!(tls=(rScan *)calloc(tlsList->nFiles,sizeof(rScan)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }

  /*loop over TLS files*/
  for(i=0;i<tlsList->nFiles;i++){
    /*read a scan*/
    tempTLS=readScan(tlsList->nameList[i]);

    /*loop through beams and see if they intersect the plot*/
    for(j=0;j<tempTLS->nBeams;j++){
      if((int)tempTLS->beam[j].nHits>maxHits){
        fprintf(stderr,"Not enough array %d\n",(int)tempTLS->beam[j].nHits);
        exit(1);
      }
      histogram[tempTLS->beam[j].nHits]++;
    }/*beam loop*/

    tempTLS=tidyTLS(tempTLS,1);
  }/*file loop*/

  if((opoo=fopen(dimage->output,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->output);
    exit(1);
  }
  for(i=1;i<maxHits;i++)total+=histogram[i];
  for(i=0;i<maxHits;i++)fprintf(opoo,"%d %d %f\n",i,histogram[i],(float)histogram[i]*100.0/(float)total);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->output);

  TIDY(histogram);
  return;
}/*readTLS*/


/*###################################################################*/
/*add one integer to the end of an array of integers more efficiently*/

uint32_t *markUint32(int length,uint32_t *jimlad,uint32_t new)
{
  if(length>0){
    if(!(jimlad=(uint32_t *)realloc(jimlad,(length+1)*sizeof(uint32_t)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
  }else{
    if(!(jimlad=(uint32_t *)calloc(length+1,sizeof(uint32_t)))){
      fprintf(stderr,"error tls allocation.\n");
      exit(1);
    }
  }
  jimlad[length]=new;
  return(jimlad);
}/*markUint32*/


/*####################################*/
/*tidy TLS structures*/

rScan *tidyTLS(rScan *tls,int nFiles)
{
  int i=0,j=0;

  if(tls){
    for(i=0;i<nFiles;i++){
      if(tls[i].beam){
        for(j=0;j<tls[i].nBeams;j++){
          TIDY(tls[i].beam[j].r);
          TIDY(tls[i].beam[j].refl);
        }
        TIDY(tls[i].beam);
      }
      TIDY(tls[i].point);
    }
    TIDY(tls);
  }
  return(NULL);
}/*tidyTLS*/


/*####################################*/
/*tidy voxel structure*/

voxStruct *tidyVox(voxStruct *vox,int nFiles)
{
  if(vox){
    TTIDY((void **)vox->hits,nFiles);
    TTIDY((void **)vox->miss,nFiles);
    TTIDY((void **)vox->mapFile,vox->nVox);
    TTIDY((void **)vox->mapPoint,vox->nVox);
    TIDY(vox->nIn);
    TIDY(vox);
  }
  return(NULL);
}/*tidyVox*/


/*####################################*/
/*read a binary file*/

rScan *readScan(char *namen)
{
  int i=0,offset=0;
  unsigned char j=0;
  uint64_t buffSize=0;   /*buffer size*/
  rScan *scan=NULL;
  double tempX=0,tempY=0,tempZ=0;
  char *buffer=NULL;
  FILE *ipoo=NULL;

  if(!(scan=(rScan *)calloc(1,sizeof(rScan)))){
    fprintf(stderr,"error scan allocation.\n");
    exit(1);
  }
  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*last 4 bytes state number of beams*/
  if(fseek(ipoo,(long)-4,SEEK_END)){ /*skip to 4 bytes from the end*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  if(fread(&scan->nBeams,sizeof(uint32_t),1,ipoo)!=1){
    fprintf(stderr,"Error reading number of points\n");
    exit(1);
  }
  fprintf(stdout,"There are %d TLS beams\n",scan->nBeams);

  /*determine file position to set size*/
  buffSize=ftell(ipoo);

  /*read data*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  buffer=challoc(buffSize,"buffer",0);   /*allocate spave*/
  if(fread(&buffer[0],sizeof(char),buffSize,ipoo)!=buffSize){  /*read beams*/
    fprintf(stderr,"Error reading point data\n");
    exit(1);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*allocate space*/
  if(!(scan->beam=(rBeam *)calloc((long)scan->nBeams,sizeof(rBeam)))){
    fprintf(stderr,"error beam allocation. Allocating %d\n",scan->nBeams);
    exit(1);
  }
  /*load into array*/
  offset=0;
  scan->xOff=scan->yOff=scan->zOff=-1.0;
  for(i=0;i<scan->nBeams;i++){
    memcpy(&(scan->beam[i].zen),&buffer[offset],8);
    scan->beam[i].zen*=M_PI/180.0; /*convert to radians*/
    offset+=8;
    memcpy(&(scan->beam[i].az),&buffer[offset],8);
    scan->beam[i].az*=M_PI/180.0; /*convert to radians*/
    offset+=8;
    memcpy(&tempX,&buffer[offset],8);
    offset+=8;
    memcpy(&tempY,&buffer[offset],8);
    offset+=8;
    memcpy(&tempZ,&buffer[offset],8);

    offset+=8;
    memcpy(&(scan->beam[i].shotN),&buffer[offset],4);
    offset+=4;
    memcpy(&(scan->beam[i].nHits),&buffer[offset],1);
    offset+=1;

    /*apply offset to coords*/
    if(scan->yOff<0.0){
      scan->xOff=tempX;
      scan->yOff=tempY;
      scan->zOff=tempZ;
    }
    scan->beam[i].x=(float)(tempX-scan->xOff);
    scan->beam[i].y=(float)(tempY-scan->yOff);
    scan->beam[i].z=(float)(tempZ-scan->zOff);

    if(scan->beam[i].nHits>0){
      scan->beam[i].r=falloc((int)scan->beam[i].nHits,"range",i);
      scan->beam[i].refl=falloc((int)scan->beam[i].nHits,"refl",i);
      for(j=0;j<scan->beam[i].nHits;j++){
        memcpy(&(scan->beam[i].r[j]),&buffer[offset],4);
        offset+=4;
        memcpy(&(scan->beam[i].refl[j]),&buffer[offset],4);
        offset+=4;
      }/*hit loop*/
    }/*hit check*/

  }/*load into array loop*/

  TIDY(buffer);
  return(scan);
}/*readScan*/


/*####################################*/
/*allocate voxel structure*/

voxStruct *voxAllocate(int nFiles,control *dimage)
{
  int i=0;
  voxStruct *vox=NULL;

  if(!(vox=(voxStruct *)calloc(1,sizeof(voxStruct)))){
    fprintf(stderr,"error voxel structure allocation.\n");
    exit(1);
  }

  /*note that findVoxels() needs minX maxX etc, different to dimage's minX minY etc*/
  vox->res=dimage->vRes;
  for(i=0;i<6;i++)vox->bounds[i]=dimage->bounds[i];
  vox->nX=(int)((vox->bounds[1]-vox->bounds[0])/vox->res);
  vox->nY=(int)((vox->bounds[3]-vox->bounds[2])/vox->res);
  vox->nZ=(int)((vox->bounds[5]-vox->bounds[4])/vox->res);
  vox->nVox=vox->nX*vox->nY*vox->nZ;
 
  if(!(vox->hits=(uint32_t **)calloc(nFiles,sizeof(uint32_t *)))){
    fprintf(stderr,"error voxel hit array allocation.\n");
    exit(1);
  }
  if(!(vox->miss=(uint32_t **)calloc(nFiles,sizeof(uint32_t *)))){
    fprintf(stderr,"error voxel miss array allocation.\n");
    exit(1);
  }

  for(i=0;i<nFiles;i++){
    if(!(vox->hits[i]=(uint32_t *)calloc(vox->nVox,sizeof(uint32_t)))){
      fprintf(stderr,"error voxel hit array allocation.\n");
      exit(1);
    }
    if(!(vox->miss[i]=(uint32_t *)calloc(vox->nVox,sizeof(uint32_t)))){
      fprintf(stderr,"error voxel miss array allocation.\n");
      exit(1);
    }
  }/*file loop*/

  vox->mapFile=iIalloc(vox->nVox,"voxel file map",0);
  vox->nIn=ialloc(vox->nVox,"",0);
  if(!(vox->mapPoint=(uint32_t **)calloc(vox->nVox,sizeof(uint32_t *)))){
    fprintf(stderr,"error voxel map allocation.\n");
    exit(1);
  }


  return(vox);
}/*voxAllocate*/


/*####################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control structure allocation.\n");
    exit(1);
  }

  strcpy(dimage->tlsNamen,"/home/sh563/data/teast/binary.list");
  strcpy(dimage->output,"teast.hist");
  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=-100000000.0;
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=100000000.0;

  /*TLS parameters for passing around*/
  dimage->nTLSwaves=4;  /*0 target profile, 1 visible projected area, 2 projected area, 3 point count*/
  dimage->beamRad=0.165;    /*33 cm footprints*/
  dimage->binOff=0;         /*start TLS wave from ALS wave*/
  dimage->appRefl1064=0.2;  /*apparent reflectance, at 1545 nm actually*/
  dimage->tlsTanDiv=sin(0.35/(2.0*1000.0))/cos(0.35/(2.0*1000.0));  /*tan half of 0.35 mrad*/
  dimage->tlsRad=0.035;     /*7 mm diameter beam start*/
  dimage->counter=0;        /*reset counter*/
  dimage->pointOut=0;       /*don't output points*/
  dimage->waveOut=1;


  /*centre on plot 1 for the tests*/
  dimage->bounds[0]=507201.6157387340-10.0;
  dimage->bounds[1]=225391.1395131220-10.0;
  dimage->bounds[2]=127.0679805710-10.0;
  dimage->bounds[3]=507201.6157387340+10.0;
  dimage->bounds[4]=225391.1395131220+10.0;
  dimage->bounds[5]=127.0679805710+20.0;
  dimage->vRes=1.0;


  /*denoising*/
  dimage->denoise.meanN=12.0;
  dimage->denoise.tailThresh=-1.0;
  dimage->denoise.thresh=15.0;
  dimage->denoise.minWidth=6;
  dimage->denoise.sWidth=0.0;
  dimage->denoise.psWidth=0.0;
  dimage->denoise.medLen=0;
  dimage->denoise.varNoise=0;

  /*deconvolution*/
  dimage->denoise.deconMeth=-1;     /*do not deconvolve*/
  dimage->denoise.pScale=1.0;      /*scale pulse length by*/
  dimage->denoise.maxIter=2000;     /*maximum number of iterations*/
  dimage->denoise.deChang=0.00001;  /*change between decon iterations to stop*/
  strcpy(dimage->denoise.pNamen,"/home/sh563/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  dimage->denoise.pulse=NULL;       /*pulse to deconvolve by*/
  dimage->denoise.pBins=0;          /*number of pulse bins*/


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-tls",4)){
        checkArguments(1,i,argc,"-tls");
        strcpy(dimage->tlsNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->output,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-tls name;        input binary filename list\n-output name;    output filename\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry snow_depletion -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*arg loop*/

  /*set to meanN if left blank*/
  if(dimage->denoise.tailThresh<0.0)dimage->denoise.tailThresh=dimage->denoise.meanN;


  return(dimage);
}/*readCommands*/

/*the end*/
/*####################################*/

