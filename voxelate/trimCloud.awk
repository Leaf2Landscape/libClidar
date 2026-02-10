
{
  if($1!="#"){
    x=$1;
    y=$2;
    z=$3;

    if((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)&&(z>=minZ)&&(z<=maxZ))print $0;
  }else print $0;
}

