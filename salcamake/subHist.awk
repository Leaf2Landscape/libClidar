BEGIN{
  numb=0;
  min=100000;
  max=-100000;
}


($0&&($1!="#")){
  x=$1;
  y=$2;
  z=$3;
  r=sqrt(x*x+y*y+z*z);
  azInd=$9;
  zenInd=$10;

  if((x>=minX)&&(x<=maxX)&&(y<=maxY)&&(y>=minY)&&(z<=maxZ)&&(z>=minZ)&&(azInd>=minAz)&&(azInd<=maxAz)&&(zenInd>=minZen)&&(zenInd<=maxZen)&&(r>=minR)&&(r<=maxR)){
    refl[numb]=$col;
    if(refl[numb]>max)max=refl[numb];
    if(refl[numb]<min)min=refl[numb];
    numb++;
  }
}

END{
  res=(max-min)/nBins;
  for(bin=0;bin<=nBins;bin++)n[bin]=0;
  for(i=0;i<numb;i++){
    bin=int((refl[i]-min)/(max-min)*nBins);
    n[bin]++;
  }

  for(bin=0;bin<=nBins;bin++)print (bin+0.5)*res+min,n[bin];
}

