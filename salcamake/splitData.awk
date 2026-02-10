BEGIN{
  numb=0;
  n=0;
  splitter=50000;
  name=sprintf("%s.%d.dat",root,n);
}

($0&&($1!="#")){
  print $0 >> name
  numb++;
  if(numb>=splitter){
    n++;
    name=sprintf("%s.%d.dat",root,n);
    numb=0;
  }
}

