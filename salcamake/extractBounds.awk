
{
  if($1!="#"){
    x=$1;
    y=$2;
    z=$3;
    r=sqrt(x*x+y*y+z*z);
    azInd=$9;
    zenInd=$10;
    if(NF>11){
      refl=$12;
      satFlag=$13;
    }else{
      refl=10.0;
      satFlag=0.0;
    }

    if((x>=minX)&&(x<=maxX)&&(y<=maxY)&&(y>=minY)&&(z<=maxZ)&&(z>=minZ)&&(satFlag<maxSat)&&(refl>=minRefl)&&(r>=minR)&&(r<=maxR)&&(azInd>=minAzI)&&(azInd<=maxAzI)&&(zenInd>=minZenI)&&(zenInd<=maxZenI))print $0;
  }else print $0;
}

