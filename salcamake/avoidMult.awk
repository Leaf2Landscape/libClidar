BEGIN{
  meanDN=meanRho=0.0;
  numb=0;
  for(az=minAz;az<=maxAz;az++){
    for(zen=minZen;zen<=maxZen;zen++)nBack[az,zen]=0;
  }

  output=sprintf("%s.hard",root);
  print "# 1 azInd, 2 zenInd, 3 DN, 4 rho, 5 range, 6 zen" > output
  output=sprintf("%s.soft",root);
  print "# 1 azInd, 2 zenInd, 3 DN, 4 rho, 5 range, 6 zen" > output
}

($0&&($1!="#")){
  r=sqrt($1*$1+$2*$2+$3*$3);
  az=$9;
  zen=$10;

  if((r>=minR)&&(r<=maxR)){  # a leaf hit
    z[numb]=$5;
    range[numb]=r;
    rho[numb]=$14;
    DN[numb]=$12;
    azI[numb]=az;
    zenI[numb]=zen;
    numb++;
  }else{
    DNback[az,zen]=$12;
    rhoBack[az,zen]=$14;
    nBack[az,zen]++;
  }
}

END{
  nIn=0;
  for(i=0;i<numb;i++){
    az=azI[i];
    zen=zenI[i];

    if(nBack[az,zen]==0)output=sprintf("%s.hard",root);
    else                output=sprintf("%s.soft",root);
    print az,zen,DN[i],rho[i],range[i],z[i] > output
  }
}

