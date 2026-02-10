#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "tiffio.h"


typedef struct{
  int nX;
  int nY;
  int nH;
  int *contN;
  int *contH;
  float maxWhite;
  float maxRefl;
  float res;
  float *imageCyl;
  float *imageHem;
  int intCol;   /*intensity column*/
}salcaImage;


/*#######################*/
/*# Reads an ascii file #*/
/*# outputs two tiffs   #*/
/*#######################*/

int main(int argc,char **argv)
{
  int i=0;
  char namen[200];
  char outRoot[100];
  char outNamen[110];
  salcaImage *dimage=NULL;
  void arrangeImage(char *,salcaImage *);
  void setupArr(salcaImage *);
  void drawPicture(char *,int,int,float *,salcaImage *);


  strcpy(namen,"/mnt/geodesy38/nsh103/SALCA/raw/brisbane/os_geom_0_1mrad_2926/salcaTest.band.0.coords");  
  strcpy(outRoot,"quickView");
  if(!(dimage=(salcaImage *)calloc(1,sizeof(salcaImage)))){
    fprintf(stderr,"error in structure allocation.\n");
    exit(1);
  }

  dimage->res=0.06;
  dimage->maxWhite=800.0;
  dimage->intCol=12-1;   /*12th column is intensity*/


  /*read command line*/
  for(i=1;i<argc;i++){
    if(*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(namen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxWhite",9)){   
        checkArguments(1,i,argc,"-maxWhite");
        dimage->maxWhite=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-hedgehog",9)){
        dimage->intCol=4-1;   /*4th column is intensity*/
      }else if(!strncasecmp(argv[i],"-13th",5)){
        dimage->intCol=13-1;   /*13th column is intensity*/
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;       input filename .coords\n-outRoot name;     output filename root\n-res res;          pixel size in degrees\n-maxWhite refl;    maximum DN to use as white, to avoid bright spots\n-hedgehog;         4th column is intensity rather than 12th\n-13th;    use 13th column, rho for hedgehog\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry echiMake -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }


  /*allocate arrays*/
  setupArr(dimage);

  /*convert point cloud to image*/
  arrangeImage(namen,dimage);

  /*draw image*/
  sprintf(outNamen,"%s.cyl.tif",outRoot);
  drawPicture(outNamen,dimage->nX,dimage->nY,dimage->imageCyl,dimage);
  sprintf(outNamen,"%s.hem.tif",outRoot);
  drawPicture(outNamen,dimage->nH,dimage->nH,dimage->imageHem,dimage);


  if(dimage){
    TIDY(dimage->imageCyl);
    TIDY(dimage->contN);
    TIDY(dimage);
  }

  return(0);
}/*main*/


/*#######################################################*/
/*draw picture*/

void drawPicture(char *namen,int nX,int nY,float *image,salcaImage *dimage)
{
  int i=0,j=0,place=0;
  unsigned char *buff=NULL;
  TIFF *file=NULL;


  file=TIFFOpen(namen,"w");
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, nX);                /* set the width of the image*/
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, nY);               /* set the height of the image*/
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 1);   /* set number of channels per pixel*/
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);                  /* set the size of the channels*/
  TIFFSetField(file, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);  /* set the origin of the image.*/
  /*Some other essential fields to set that you do not have to understand for now.*/
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC,PHOTOMETRIC_MINISBLACK);   /*mode 1 is greyscale with 0 as black*/

  /*apparently the simplest method is to write a row at a time*/
  buff=(unsigned char *)_TIFFmalloc(nX);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(file,nX));
  for(j=0;j<nY;j++){
    for(i=0;i<nX;i++){
      place=j*nX+i;
      if(dimage->imageCyl[place]>=dimage->maxWhite)buff[i]=255;
      else buff[i]=(unsigned char)(image[j*nX+i]*255.0/dimage->maxWhite);
    }

    if(TIFFWriteScanline(file,buff,j,0)<0){
      fprintf(stderr,"Error writing\n");
      exit(1);
    }
  }

  if(buff){
    _TIFFfree(buff);
    buff=NULL;
  }
  TIFFClose(file);

  fprintf(stdout,"Drawn to %s\n",namen);
  return;
}/*drawPicture*/


/*#######################################################*/
/*set up arrays*/

void setupArr(salcaImage *dimage)
{
  int i=0;

  /*cylindrical*/
  dimage->nX=(int)(360.0/dimage->res+1);
  dimage->nY=(int)(192.0/dimage->res+1);   /*-96 to 96*/
  dimage->imageCyl=falloc(dimage->nX*dimage->nY,"cylindrical image",0);
  dimage->contN=ialloc(dimage->nX*dimage->nY,"counter",0);
  for(i=dimage->nX*dimage->nY-1;i>=0;i--){
    dimage->imageCyl[i]=0.0;
    dimage->contN[i]=0;
  }

  /*hemispherical*/
  dimage->nH=(int)(180.0/dimage->res+1);
  dimage->imageHem=falloc(dimage->nH*dimage->nH,"hemispherical image",0);
  dimage->contH=ialloc(dimage->nH*dimage->nH,"counter",0);
  for(i=dimage->nH*dimage->nH-1;i>=0;i--){
    dimage->imageHem[i]=0.0;
    dimage->contH[i]=0;
  }

  return;
}/*setupArr*/


/*#######################################################*/
/*read data and arrange image*/

void arrangeImage(char *namen,salcaImage *dimage)
{
  int i=0,j=0;
  int con=0;
  int place=0;
  int number=0;
  float zen=0,az=0;
  float refl=0;
  float r=0;    /*distance from image centre*/
  char line[600];
  char temp[15][100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*get brightness bounds*/
  dimage->maxRefl=-100.0;
  con=0;
  while(fgets(line,600,ipoo)!=NULL){
    if(sscanf(line,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",temp[0],temp[1],temp[2],temp[3],\
                temp[4],temp[5],temp[6],temp[7],temp[8],temp[9],temp[10],temp[11],temp[12],temp[13],temp[14])==15){
      if(strncasecmp(temp[0],"#",1)){  /*not a comment*/
        number=atoi(temp[7]);

        if(number==1){
          az=atof(temp[5]);
          zen=atof(temp[4]);
          refl=atof(temp[dimage->intCol]);

          if(az<-180.0)az+=360.0;
          else if(az>180.0)az-=360.0;

          /*cylindrical image*/
          i=(int)((az+180.0)/dimage->res);
          j=(int)((zen+95.0)/dimage->res);

          if((i>=0)&&(i<dimage->nX)&&(j>=0)&&(j<dimage->nY)){
            place=j*dimage->nX+i;
            if(refl>dimage->maxRefl)dimage->maxRefl=refl;

            dimage->imageCyl[place]+=refl;
            dimage->contN[place]++;
          }

          /*hemispherical image*/
          if(zen<=90.0){
            r=zen/90.0;
            i=(int)((r*sin(az*M_PI/180.0)+1.0)*(float)dimage->nH/2.0);
            j=(int)((r*cos(az*M_PI/180.0)+1.0)*(float)dimage->nH/2.0);
            if((i>=0)&&(i<dimage->nH)&&(j>=0)&&(j<dimage->nH)){
              place=j*dimage->nH+i;
              dimage->imageHem[place]+=refl;
              dimage->contH[place]++;
            }

          }/*horizon check*/
        }/*return number check*/
        con++;
      }/*comment check*/
    }/*line splitting*/
  }/*file reading loop*/


  if(dimage->maxWhite>dimage->maxRefl)dimage->maxWhite=dimage->maxRefl;

  /*normalise*/
  for(i=dimage->nX*dimage->nY-1;i>=0;i--)if(dimage->contN[i])dimage->imageCyl[i]/=(float)dimage->contN[i];
  for(i=dimage->nH*dimage->nH-1;i>=0;i--)if(dimage->contH[i])dimage->imageHem[i]/=(float)dimage->contH[i];

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*arrangeImage*/


/*the end*/
/*#######################################################*/

