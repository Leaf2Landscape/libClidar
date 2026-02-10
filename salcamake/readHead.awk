
{
  if($1=="Range:")                   maxR=$NF;
  else if(($1=="AZ")&&($2=="step"))  azStep=$NF;
  else if(($1=="AZ")&&($2=="steps:"))nAz=$NF;
  else if($1=="Start")               azStart=$NF;
  else if($1=="Finish")              maxZen=$NF;
}

END{
  printf("%s\n%s\n%s\n%s\n%s\n",maxR,azStep,nAz,azStart,maxZen);
}

