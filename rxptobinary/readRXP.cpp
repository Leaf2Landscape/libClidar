#include <riegl/scanlib.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <cstring>
#include <cstdlib>
#include <hdf5.h>
//#include <H5Cpp.h>
//#include <libLidarHDF.h>

using namespace scanlib;
using namespace std;

/*thresholds for recording points/beams*/
float maxZen,maxDev;


/*#########################################################*/
/* A class that reads the .rxp files*/

class reader : public pointcloud  /*we have inherited pointcloud public members which in turn inherits "compressed_packets"*/
{
private:
  ostream& o;
  ostream& oAsc;
  ostream& oSSCI;
  unsigned long line;
  int doneThis;
  float zen,az;      /*beam angle in degrees*/
  float range[20],refl[20];
  float dev;         /*deviation of pulse shape*/
  float scanCent[3];
  double trans[4][4];      /*4*4 translation matrix for geolocation*/
  void transform();
  uint8_t nHits;
  /*for the SSCI*/
  void calcSSCI();
  void writeSSCI();
  int nSSCIshot;      /*number of shots in this segment*/
  int nInSSCI;        /*number of shots in this segment*/
  float ssciFrac;     /*FRAC parameter for SSCI*/
  float ssciA;        /*area of the polygon*/
  float ssciP;        /*perimeter of the polygon*/
  float lastZen;      /*to keep track of where we are in scan*/
  float lastAz;       /*to keep track of where we are in scan*/
  double lastX;       /*to keep track of where we are in scan*/
  double lastY;       /*to keep track of where we are in scan*/
  double lastZ;       /*to keep track of where we are in scan*/
  double lastR;       /*to keep track of where we are in scan*/
  int lastHit;        /*to keep track of where we are in scan*/
  float ssciMinZen;   /*minimum zenith for this segment*/
public:
  uint32_t shotN,totShot;         // total number of shots
  double coordOffset[3];
  double bounds[6];
  int ascOut,binOut,hdfOut,findSSCI;  // switches copied from control structure
  int nDec;
  /*arrays to hold data for HDF5 file*/
  float *zenArr,*azArr;
  float *x,*y,*z;
  uint8_t *nHitsArr;
  uint32_t *shotNArr,totHits;
  float *rangeArr,*reflArr;

  // the constructor function
  reader(ostream& o_,ostream& oAsc_,ostream& oSSCI_)
        : pointcloud(false) // set this to true if you need gps aligned timing. What happens in the poincoud constructor?
        , o(o_)             // initialises o as o_, which is what is passed to the class constructor.
        , oAsc(oAsc_)
        , oSSCI(oSSCI_)
  {
    if(o)o.precision(10);    /*sets the number of decimal places of output*/
    line=0;
    shotN=0;
    totShot=0;
    totHits=0;
    zenArr=NULL;  /*force all arrays to NULL*/
    azArr=NULL;
    x=NULL;
    y=NULL;
    z=NULL;
    nHitsArr=NULL;
    shotNArr=NULL;
    rangeArr=NULL;
    reflArr=NULL;
    ssciFrac=0.0;
    ssciA=0.0;
    ssciP=0.0;
    ssciMinZen=10000000000.0;
    nInSSCI=0;
    nSSCIshot=0;
  }/*creator function*/

  /*This call is invoked for every pulse, even if there is no retur, by importer, before on_echo_transformed()*/
  void on_shot()   // this seems to be an overloading of a pointcloud member function.
  {
    /*transform vector and origin to geolocate*/
    transform();

    if(shotN==0){
      coordOffset[0]=beam_origin[0];
      coordOffset[1]=beam_origin[1];
      coordOffset[2]=beam_origin[2];
    }

    zen=(float)(atan2(sqrt(beam_direction[0]*beam_direction[0]+beam_direction[1]*beam_direction[1]),beam_direction[2])*180.0/M_PI);
    az=(float)(atan2(beam_direction[0],beam_direction[1])*180.0/M_PI);
    if((zen>361.0)||(zen<-361.0)||(az>361.0)||(az<-361.0))fprintf(stderr,"Angle error %f %f\n",zen,az);
    scanCent[0]=(float)(beam_origin[0]-coordOffset[0]);
    scanCent[1]=(float)(beam_origin[1]-coordOffset[1]);
    scanCent[2]=(float)(beam_origin[2]-coordOffset[2]);

    nHits=0;

    return;
  }/*on_shot*/

  /*Invoked for every return by importer, after on_shot()*/
  void on_echo_transformed(echo_type echo)   // this seems to be an overloading of a pointcloud member function.
  {
    target& t(targets[target_count-1]);

    /*check deviation before recording*/
    if(t.deviation<maxDev){
      if(nHits>=20){
        fprintf(stderr,"Too many hits in this beam\n");
        exit(1);
      }

      range[nHits]=(float)t.echo_range;
      refl[nHits]=pow(10.0,t.reflectance/10.0);

      /*write point cloud*/
      if(ascOut&&((shotN%nDec)==0)){
        dev=t.deviation;
        writePoint();
      }

      nHits++;
    }/*deviation checker*/

    return;
  }/*on_echo_transformed*/


  /*Check if we're done with a pulse*/
  void on_shot_end()   /*called at end of pulse*/
  {
    /*this is where we will call the voxelisation program*/
    if((zen<maxZen)&&((totShot%nDec)==0)){
      if(binOut&&!hdfOut)writeBeam();   /*write data to binary file*/
      else               packHDF();     /*store data in HDF5 arrays*/

      if(findSSCI)calcSSCI();           /*caclulate SSCI metric*/
      shotN++;
    }
    totShot++;

    return;
  }/*on_shot_end*/

  /*reset markers*/
  void reset(){doneThis=1;}

  /*read translation matrix*/
  void readTrans(char *,char *);

  /*write a beam to binary output*/
  void writeBeam();

  /*pack data in to HDF5 array*/
  void packHDF();

  /*write data to HDF5*/
  void writeHDFfile(char *);

  /*write point to ASCII file*/
  void writePoint();

  /*end of class*/
};/*reader class*/


/*#######################################################################*/
/*write HDF file*/

void reader::writeHDFfile(char *outNamen)
{
  hid_t file;         /* Handles */
  char varNamen[100];
  void writeComp1dFloatHDF5(hid_t,char *,float *,int);
  void writeComp1dUint8HDF5(hid_t,char *,uint8_t *,int);
  void writeComp1dUint32HDF5(hid_t,char *,uint32_t *,int);
  void write1dDoubleHDF5(hid_t,char *,double *,int);
  void write1dUint32HDF5(hid_t,char *,uint32_t *,int);


  /*open file*/
  file=H5Fcreate(outNamen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  /*write header parts*/
  strcpy(&varNamen[0],"XCENT");
  write1dDoubleHDF5(file,&varNamen[0],&coordOffset[0],1);
  strcpy(&varNamen[0],"YCENT");
  write1dDoubleHDF5(file,&varNamen[0],&coordOffset[1],1);
  strcpy(&varNamen[0],"ZCENT");
  write1dDoubleHDF5(file,&varNamen[0],&coordOffset[2],1);
  strcpy(&varNamen[0],"TOTSHOT");
  write1dUint32HDF5(file,&varNamen[0],&shotN,1);
  strcpy(&varNamen[0],"TOTHITS");
  write1dUint32HDF5(file,&varNamen[0],&totHits,1);


  /*write data in to file*/
  strcpy(&varNamen[0],"X0");
  writeComp1dFloatHDF5(file,&varNamen[0],x,shotN);
  strcpy(&varNamen[0],"Y0");
  writeComp1dFloatHDF5(file,&varNamen[0],y,shotN);
  strcpy(&varNamen[0],"Z0");
  writeComp1dFloatHDF5(file,&varNamen[0],z,shotN);
  strcpy(&varNamen[0],"ZEN");
  writeComp1dFloatHDF5(file,&varNamen[0],zenArr,shotN);
  strcpy(&varNamen[0],"AZ");
  writeComp1dFloatHDF5(file,&varNamen[0],azArr,shotN);
  strcpy(&varNamen[0],"NHITS");
  writeComp1dUint8HDF5(file,&varNamen[0],nHitsArr,shotN);
  strcpy(&varNamen[0],"SHOTN");
  writeComp1dUint32HDF5(file,&varNamen[0],shotNArr,shotN);
  strcpy(&varNamen[0],"RANGE");
  writeComp1dFloatHDF5(file,&varNamen[0],rangeArr,totHits);
  strcpy(&varNamen[0],"REFL");
  writeComp1dFloatHDF5(file,&varNamen[0],reflArr,totHits);


  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }
  fprintf(stdout,"HDF5 data written to %s\n",outNamen);

  /*tidy up*/
  if(x){free(x);x=NULL;}
  if(y){free(y);y=NULL;}
  if(z){free(z);z=NULL;}
  if(zenArr){free(zenArr);zenArr=NULL;}
  if(azArr){free(azArr);azArr=NULL;}
  if(rangeArr){free(rangeArr);rangeArr=NULL;}
  if(reflArr){free(reflArr);reflArr=NULL;}
  if(nHitsArr){free(nHitsArr);nHitsArr=NULL;}
  if(shotNArr){free(shotNArr);shotNArr=NULL;}

  return;
}/*writeHDF*/


/*#######################################################################*/
/*write a binary beam*/

void reader::writeBeam()
{
  int i=0,len=25+8*nHits;
  int offset=0;
  char data[len];

  /*out data into char array and write it*/
  memcpy(&(data[0]),&zen,4);
  memcpy(&(data[4]),&az,4);
  memcpy(&(data[8]),&(scanCent[0]),4);
  memcpy(&(data[12]),&(scanCent[1]),4);
  memcpy(&(data[16]),&(scanCent[2]),4);
  memcpy(&(data[20]),&shotN,4);
  memcpy(&(data[24]),&nHits,1);

  /*if((nHits==0)&&((zen>90.0)||(zen<=-90.0))){
    cout << zen<<" "<<az<< " " << scanCent[0] << " " << scanCent[1] << " " << scanCent[2] << endl;
  }*/

  offset=25;
  for(i=0;i<nHits;i++){
    memcpy(&(data[offset]),&(range[i]),4);
    offset+=4;
    memcpy(&(data[offset]),&(refl[i]),4);
    offset+=4;
  }
  o.write(data,len);

  return;
}/*writeBeam*/


/*#######################################################################*/
/*pack data in to HDF5 array*/

void reader::packHDF()
{
  uint8_t i;
  uint8_t *markUint8(uint32_t,uint8_t *,uint8_t);
  uint32_t *markUint32(uint32_t,uint32_t *,uint32_t);
  float *markFloat(uint32_t,float *,float);
  double *markDo(uint32_t,double *,double);

  /*increase array sizes and add values*/
  zenArr=markFloat(shotN,zenArr,zen);
  azArr=markFloat(shotN,azArr,az);
  x=markFloat(shotN,x,scanCent[0]);
  y=markFloat(shotN,y,scanCent[1]);
  z=markFloat(shotN,z,scanCent[2]);
  nHitsArr=markUint8(shotN,nHitsArr,nHits);
  shotNArr=markUint32(shotN,shotNArr,shotN);

  /*loop over variable hits*/
  for(i=0;i<nHits;i++){
    rangeArr=markFloat((uint32_t)totHits,rangeArr,range[i]);
    reflArr=markFloat((uint32_t)totHits,reflArr,refl[i]);
    totHits++;
  }

  return;
}/*packHDF*/


/*#######################################################################*/
/*write an ASCII point*/

void reader::writePoint()
{
  double x=0,y=0,z=0;
  double zenRad=0,azRad=0;

  zenRad=(double)zen*M_PI/180.0;
  azRad=(double)az*M_PI/180.0;

  /*is this correct?*/
  x=sin(zenRad)*sin(azRad)*(double)range[nHits]+(double)scanCent[0]+coordOffset[0];
  y=sin(zenRad)*cos(azRad)*(double)range[nHits]+(double)scanCent[1]+coordOffset[1];
  z=cos(zenRad)*(double)range[nHits]+(double)scanCent[2]+coordOffset[2];

  if((x>=bounds[0])&&(y>=bounds[1])&&(z>=bounds[2])&&(x<=bounds[3])&&(y<=bounds[4])&&(z<=bounds[5])){
    oAsc << std::setprecision(10) << x << " " << y << " " << z << " " << refl[nHits] << " " << zen << " " << az << " " << range[nHits] << " hits " << (float)nHits << " " << dev << endl;
  }

  return;
}/*writePoint*/


/*#######################################################################*/
/*read transformation matrix*/

void reader::readTrans(char *namen,char *globNamen)
{
  int i=0,j=0,k=0;
  ifstream ipoo(namen); /*file is opened by constructor*/
  double tempTrans[4][4],newTrans[4][4];

  /*do we need to apply a local transformation?*/
  if(strncasecmp(namen,"none",4)){
    /*check that file is opened*/
    if(!ipoo.is_open()){
      cerr << "Input file not open " << namen << endl;
      exit(1);
    }

    /*read file*/
    i=0;
    while((!ipoo.eof())&&(i<4)){
      ipoo >> trans[i][0] >> trans[i][1] >> trans[i][2] >> trans[i][3];
      i++;
    }
    ipoo.close();
  }else{    /*unity matrix*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        trans[i][j]=0.0;
      }
    }
    for(i=0;i<4;i++)trans[i][i]=1.0;
  }

  /*do we need to apply a global transformation too?*/
  if(strncasecmp(globNamen,"none",4)){
    ifstream gpoo(globNamen); /*file is opened by constructor*/
    /*read global ,matrix*/
    if(!gpoo.is_open()){
      cerr << "Input file not open " << globNamen << endl;
      exit(1);
    }
    i=0;
    while((!gpoo.eof())&&(i<4)){
      gpoo >> tempTrans[i][0] >> tempTrans[i][1] >> tempTrans[i][2] >> tempTrans[i][3];
      i++;
    }
    gpoo.close();

    /*apply global matrix*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        newTrans[i][j]=0.0;
        for(k=0;k<4;k++){
          newTrans[i][j]+=tempTrans[i][k]*trans[k][j];
        }
      }
    }

    /*copy over*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        trans[i][j]=newTrans[i][j];
      }
    }
  }/*global matrix if needed*/

  return;
}/*readTrans*/


/*#######################################################################*/
/*calculate SSCI if at end of track*/

void reader::calcSSCI()
{
  float sep=0;
  double r=0,a=0,b=0,c=0,s=0;
  double x=0,y=0,z=0,p=0;
  double zenRad=0,azRad=0;

  /*record highest zenith of segment*/
  if(zen<ssciMinZen)ssciMinZen=zen;

  /*add up area and last hit*/
  zenRad=(double)zen*M_PI/180.0;
  azRad=(double)az*M_PI/180.0;

    /*is this correct?*/
  if(nHits>0)r=(double)range[nHits-1];
  else       r=100.0;

  x=sin(zenRad)*sin(azRad)*r+(double)scanCent[0]+coordOffset[0];
  y=sin(zenRad)*cos(azRad)*r+(double)scanCent[1]+coordOffset[1];
  z=cos(zenRad)*r+(double)scanCent[2]+coordOffset[2];

  /*only if this is the second point*/
  if(nInSSCI>=1){
    /*add up perimeter*/
    p=sqrt((x-lastX)*(x-lastX)+(y-lastY)*(y-lastY)+(z-lastZ)*(z-lastZ));
    ssciP+=p;
    /*segment area from Heron's forumla*/
    a=lastR;
    b=r;
    c=p;
    s=(a+b+c)/2.0;
    ssciA+=sqrt(s*(s-1)*(s-b)*(s-c));
  }
  nInSSCI++;

  /*if we have got to the end of a segment, record*/
  if(shotN>0)sep=sqrt((lastZen-zen)*(lastZen-zen)+(lastAz-az)*(lastAz-az));
  else       sep=0.0;
  if(sep>1.0*nDec){  /*then we are in a new segment calculate SSCI*/
    ssciFrac=2.0*log(0.25*ssciP)/log(ssciA);

    /*write results if more than 1 beam*/
    if(nInSSCI>1)writeSSCI();

    /*reset counters*/
    ssciP=ssciA=0.0;
    nInSSCI=0;
    ssciMinZen=100000000.0;
    nSSCIshot++;
  }

  lastZen=zen;
  lastAz=az;
  lastX=x;
  lastY=y;
  lastZ=z;
  lastR=r;

  return;
}/*reader::calcSSCI*/


/*#######################################################################*/
/*write SSCI results*/

void reader::writeSSCI()
{
  oSSCI << std::setprecision(10) << ssciFrac << "," << ssciP << "," << ssciA << "," << ssciMinZen << "," <<  nInSSCI << "," << nSSCIshot << endl;

  return;
}/*reader::writeSSCI*/


/*#######################################################################*/
/*apply transformation*/

void reader::transform()
{
  double *tempDir=NULL,*tempCent=NULL;
  double *setBits(double *,double **);
  double *rotateVect(double[3],double[4][4]);

  /*rotate*/
  tempDir=rotateVect(beam_direction,trans);
  tempCent=rotateVect(beam_origin,trans);

  /*translate and put into original vectors*/
  beam_origin[0]=tempCent[0]+trans[0][3];
  beam_origin[1]=tempCent[1]+trans[1][3];
  beam_origin[2]=tempCent[2]+trans[2][3];
  beam_direction[0]=tempDir[0];
  beam_direction[1]=tempDir[1];
  beam_direction[2]=tempDir[2];

  /*free arrays*/
  if(tempDir)delete[] tempDir;
  if(tempCent)delete[] tempCent;

  return;
}/*transform*/


/*#######################################################################*/
/*rotate a vector*/

double *rotateVect(double vect[3],double matrix[4][4])
{
  int i=0,j=0;
  double *rotated=NULL;

  rotated=new double[3];

  for(i=0;i<3;i++){
    rotated[i]=0.0;
    for(j=0;j<3;j++)rotated[i]+=vect[j]*matrix[i][j];
  }

  return(rotated);
}/*rotateVect*/


/*#######################################################################*/
/*control structure*/
 
typedef struct {
  char inName[200];    // I'll need a list of these eventually.
  char transName[200]; // transformation file
  char globName[200];  // global transformation file
  char outNamen[200];  // output filename
  char ascNamen[200];  // ascii output filename
  char ssciNamen[200]; // ssci output filename
  char ascOut;         // ascii output switch
  char binOut;         // binary output switch
  char hdfOut;         // hdf output switch
  char findSSCI;       // Calculate SSCI complexity index
  int nDec;            // amount to decimate by
  double bounds[6];    // bounds for ASCII output
}control;


/*#######################################################################*/
/*main*/

int main(int argc, char** argv)
{
  int i=0;
  shared_ptr<basic_rconnection> rc;   /*the file object? Exciting bits are in here*/
  control *dimage=NULL;
  control *readCommands(int,char **);
  ofstream oAsc;
  ofstream opoo;
  ofstream oSSCI;

  /*read the command line*/
  dimage=readCommands(argc,argv);

  if(dimage->binOut&&(!dimage->hdfOut))opoo.open(dimage->outNamen,ios::binary);
  if(dimage->ascOut)oAsc.open(dimage->ascNamen);
  if(dimage->findSSCI){
    oSSCI.open(dimage->ssciNamen);
    oSSCI << "ssciFrac,ssciP,ssciA,ssciMinZen,nInSSCI,nSSCIshot" << endl;
  }


  /*open the file*/
  cout << "# Reading " << dimage->inName << "\n";
  rc=basic_rconnection::create(dimage->inName);
  rc->open();

  /*decode the file*/
  decoder_rxpmarker dec(rc);     /*needs to be declared after rc is opened for constructor*/
  reader            imp(opoo,oAsc,oSSCI);   /*this is a declaration*/
  buffer            buf;

  /*copy switches to class for later*/
  imp.ascOut=dimage->ascOut;
  imp.binOut=dimage->binOut;
  imp.hdfOut=dimage->hdfOut;
  imp.nDec=dimage->nDec;
  imp.findSSCI=dimage->findSSCI;
  for(i=0;i<6;i++)imp.bounds[i]=dimage->bounds[i];

  /*read the transformation file*/
  imp.readTrans(dimage->transName,dimage->globName);

  /*loop over data packets. Functions are called within the inherited classes*/
  for(dec.get(buf);!dec.eoi();dec.get(buf)){ /*loops through the data packets*/
    imp.reset();   /*reset markers*/
    /*many shots are called in the packets, so calls need to be in there*/
    imp.dispatch(buf.begin(), buf.end());   /*dispatch is a member of basic_packets*/
  }

  /*close file*/
  if(dimage->binOut&&!dimage->hdfOut){
    /*write out total number of beams written*/
    opoo.write((char *)(&imp.coordOffset[0]),8);
    opoo.write((char *)(&imp.coordOffset[1]),8);
    opoo.write((char *)(&imp.coordOffset[2]),8);
    opoo.write((char *)(&imp.shotN),4);
    opoo.close();
    cout << "Binary to " << dimage->outNamen << endl;
  }else if(dimage->hdfOut)imp.writeHDFfile(dimage->outNamen);

  if(dimage->findSSCI){
    oSSCI.close();
    cout << "SSCI to " << dimage->ssciNamen << endl;
  }

  if(dimage->ascOut){
    oAsc.close();
    cout << "ASCII to " << dimage->ascNamen << endl;
  }

  /*close file and tidy arrays*/
  rc->close();
  if(dimage){
    free(dimage);
    dimage=NULL;
  }
  return(0);
}/*main*/


/*#######################################################################*/
/*read the command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  void checkArguments(int,int,int,string);

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  /*defaults*/
  strcpy(dimage->inName,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/2014-08-05.001.riproject/ScanPos002/140805_094556.rxp");
  strcpy(dimage->transName,"none");
  strcpy(dimage->outNamen,"teast.bin");
  strcpy(dimage->ascNamen,"teast.pts");
  strcpy(dimage->ssciNamen,"ssci.txt");
  strcpy(dimage->globName,"none");
  dimage->ascOut=0;
  dimage->binOut=1;
  dimage->hdfOut=0;
  dimage->findSSCI=0;
  dimage->nDec=1;
  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=-100000000.0;
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=100000000.0;
  maxZen=10000000.0;
  maxDev=10000000.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-trans",6)){
        checkArguments(1,i,argc,"-trans");
        strcpy(dimage->transName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-globMat",8)){
        checkArguments(1,i,argc,"-globMat");
        strcpy(dimage->globName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-ascii",6)){
        checkArguments(1,i,argc,"-ascii");
        dimage->ascOut=1;
        dimage->binOut=0;
        strcpy(dimage->ascNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-dec",4)){
        checkArguments(1,i,argc,"-dec");
        dimage->nDec=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxDev",7)){
        checkArguments(1,i,argc,"-maxDev");
        maxDev=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noBin",6)){
        dimage->binOut=0;
      }else if(!strncasecmp(argv[i],"-hdf",4)){
        dimage->hdfOut=1;
        dimage->binOut=1;
      }else if(!strncasecmp(argv[i],"-ssci",9)){
        checkArguments(1,i,argc,"-ssci");
        strcpy(dimage->ssciNamen,argv[++i]);
        dimage->findSSCI=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;       input rxp filename\n-trans name;       transformation filename\n-globMat name;     global coordinate matrix if needed. Unity otherwise\n-output name;      output filename\n-hdf;              write output in HDF format\n-ascii name;       ascii output and filename\n-ssci name;       calculate SSCI complexity metric\n-noBin;            don't output binary\n-maxZen zen;       max zen to trust. Degrees\n-maxDev dev;       maximum return shape deviation to accept\n-dec n;            rate to decimate pointcloud by\n-bounds minX minY minZ maxX maxY maxZ;    bounds for ASCII output\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry readRXP -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/


/*#########################################################################*/
/*check the number of command line arguments*/

void checkArguments(int numberOfArguments,int thisarg,int argc,string option)
{
  int i=0;

  for(i=0;i<numberOfArguments;i++){
    if(thisarg+1+i>=argc){
      cerr << "error in number of arguments for " << option << " option: " << numberOfArguments << "required\n";
      exit(1);
    }
  }
  return;
}/*checkArguments*/


/*###################################################################*/
/*add one integer to the end of an array of floats*/

float *markFloat(uint32_t length,float *jimlad,float newVal)
{

  if(length>0){
    if(!(jimlad=(float *)realloc(jimlad,(size_t)((uint64_t)(length+1)*(uint64_t)sizeof(float))))){
      fprintf(stderr,"Error in float reallocation, %lu\n",(length+1)*sizeof(float));
      exit(1);
    }
  }else{
    if(!(jimlad=(float *)calloc((size_t)((uint64_t)length+1),sizeof(float)))){
      fprintf(stderr,"error in float buffer allocation.\n");
      exit(1);
    }
  }
  jimlad[length]=newVal;
  return(jimlad);
}/*markFloat*/


/*#############################################################*/
/*add one integer to the end of an array of doubles*/

double *markDo(uint32_t length,double *jimlad,double newVal)
{
  if(length>0){
    if(!(jimlad=(double *)realloc(jimlad,(length+1)*sizeof(double)))){
      fprintf(stderr,"Error in double reallocation %lu\n",(length+1)*sizeof(double));
      exit(1);
    }
  }else{
    if(!(jimlad=(double *)calloc(length+1,sizeof(double)))){
      fprintf(stderr,"error in float buffer allocation.\n");
      exit(1);
    }
  }
  jimlad[length]=newVal;
  return(jimlad);
}/*markDo*/


/*#############################################################*/
/*add one integer to the end of an array of uint8*/

uint8_t *markUint8(uint32_t length,uint8_t *jimlad,uint8_t newVal)
{
  if(length>0){
    if(!(jimlad=(uint8_t *)realloc(jimlad,(length+1)*sizeof(uint8_t)))){
      fprintf(stderr,"Error in double reallocation %lu\n",(length+1)*sizeof(double));
      exit(1);
    }
  }else{
    if(!(jimlad=(uint8_t *)calloc(length+1,sizeof(uint8_t)))){
      fprintf(stderr,"error in float buffer allocation.\n");
      exit(1);
    }
  }
  jimlad[length]=newVal;
  return(jimlad);
}/*uint8_t*/


/*#############################################################*/
/*add one integer to the end of an array of uint32*/

uint32_t *markUint32(uint32_t length,uint32_t *jimlad,uint32_t newVal)
{
  if(length>0){
    if(!(jimlad=(uint32_t *)realloc(jimlad,(length+1)*sizeof(uint32_t)))){
      fprintf(stderr,"Error in double reallocation %lu\n",(length+1)*sizeof(double));
      exit(1);
    }
  }else{
    if(!(jimlad=(uint32_t *)calloc(length+1,sizeof(uint32_t)))){
      fprintf(stderr,"error in float buffer allocation.\n");
      exit(1);
    }
  }
  jimlad[length]=newVal;
  return(jimlad);
}/*uint32_t*/


/*####################################################*/
/*write a compressed 1D float array*/

void writeComp1dFloatHDF5(hid_t file,char *varName,float *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_FLOAT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return;
}/*writeComp1dFloatHDF5*/


/*####################################################*/
/*write a 1D uint8 array*/

void writeComp1dUint8HDF5(hid_t file,char *varName,uint8_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_UCHAR);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return;
}/*writeComp1dUint8HDF5*/


/*####################################################*/
/*write a 1D uint32 array*/

void writeComp1dUint32HDF5(hid_t file,char *varName,uint32_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_UINT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return;
}/*writeComp1dUint32HDF5*/


/*####################################################*/
/*write a 1D double array*/

void write1dDoubleHDF5(hid_t file,char *varName,double *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_DOUBLE);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return;
}/*write1dDoubleHDF5*/


/*####################################################*/
/*write a 1D uint32 array*/

void write1dUint32HDF5(hid_t file,char *varName,uint32_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_UINT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);

  return;
}/*write1dUint32HDF5*/


/*the end*/
/*#######################################################################*/

