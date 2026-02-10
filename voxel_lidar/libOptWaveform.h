/*#######################################*/
/*# Copyright 2014-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of the lidar voxel library, compareWaves.          #*/
/*#                                                                      #*/
/*# compareWaves is free software: you can redistribute it and/or modify #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# compareWaves is distributed in the hope that it will be useful,      #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#  You should have received a copy of the GNU General Public License   #*/
/*#  along with compareWaves.  If not, see <http://www.gnu.org/licenses/>#*/
/*########################################################################*/


/*##########################################*/
/*options and control variables*/

typedef struct{
  /*universal parameters*/
  int nFiles;     /*number of input files*/
  float res;
  /*options*/
  char weight;         /*weight to avoid subterranean*/
  char optSmoo;        /*optimise smoothing*/
  char optPsmoo;       /*optimise pre-smoothing*/
  char fixPS;          /*fix pulse length scaling*/
  float denWidth;      /*smoothing width when testing denoising accuracy*/
  char testAtten;      /*optimise to the attenuation corrected waveform*/
  float flWeight;      /*optimisation weight*/
  char doGBIC;         /*apply GBIC*/
  char fixGap;         /*fix min gap correction*/
  /*input/output*/
  FILE *errOut;        /*error output file*/
}optControl;


/*##########################################*/
/*data structure*/

typedef struct{
  /*raw data*/
  float *als;             /*als waveform*/
  unsigned char *rawALS;  /*raw ALS waveform*/
  float gbic;             /*intensity correction*/
  int nALS;               /*number of ALS bins*/
  int nTLS;               /*number of TLS points*/
  tlsPoint *tls;          /*pointer to tls structure*/
  float aE;               /*ALS waveform energy*/
  double x0;              /*point coordinate*/
  double y0;              /*point coordinate*/
  double z0;              /*point coordinate*/

  /*reconstructed waveforms*/
  float **tlsWave;   /*0 ALS, 1 visible, 2 area*/
  float *smoothTLS;  /*smoothed TLS to optimise to*/
  float tE;          /*visible TLS wave energy*/
  float *denoised;   /*denoised ALS wave*/
  int tStart;        /*start bin of TLS wave*/
  int tEnd;          /*end bin of TLS wave*/

  /*TLS parameters*/
  lidVoxPar tlsPar;
}dataStruct;


/*##########################################*/
/*overarching data structure*/

typedef struct{
  dataStruct *dArr;   /*pointer to array of data structures*/
  int thisTLS;        /*tls index for keeping track when optimising TLS parameters*/

  /*ALS pulse*/
  float *pulse;
  float res;
  int pBins;
  int maxBin;

  /*ALS denoising*/
  denPar den;        /*ALS denoising parameters*/

  /*pointer to control to pass options*/
  optControl *dimage;
}overStruct;

/*the end*/
/*##########################################*/

