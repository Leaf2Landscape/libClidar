#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"
#include "scanifc.h"


/*####################*/
/*# Learn to use the #*/
/*# RiVLib libraries #*/
/*####################*/


int main()
{
  char inName[200];
  scanifc_uint16_t major=0,minor=0,build=0;
  scanifc_bool sync_to_pps=0;
  point3dstream_handle h3ds = 0; /* This is the handle to the data stream. */



  strcpy(inName,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/2014-08-05.004.riproject/ScanPos001/140805_134208.rxp");

  scanifc_get_library_version(&major,&minor,&build);
  fprintf(stdout,"Versions %d %d %d\n",major,minor,build);

  if(scanifc_point3dstream_open_with_logging(inName,sync_to_pps,"log.rxp", &h3ds)){
    fprintf(stderr, "Cannot open: %s\n",inName);
    //print_last_error();
    return(1);
  }


  return(0);
}/*main*/

