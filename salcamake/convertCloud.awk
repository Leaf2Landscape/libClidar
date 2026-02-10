BEGIN{
  found=0;
  p=0;
}


($0&&($1!="#")){
  n=$8;

  if(n==1){
    if(found){ # output points
      printf("%g,%g,%g,%g",x[0],y[0],z[0],I[0]);        # first return
      if(p>1)printf("%g,%g,%g,%g",x[1],y[1],z[1],I[1]); # second return
      else   printf("0,0,9999,9999"); # second return
      if(p>2)printf("%g,%g,%g,%g",x[p-1],y[p-1],z[p-1],I[p-1]); # last return
      else   printf("0,0,9999,9999"); # last return
      printf("\n");
      p=0;
    }else found=1;
    x[p]=$1;
    y[p]=$2;
    z[p]=$3;
    I[p]=$12;
    p++;
  }
}

