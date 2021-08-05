GCM_wrap<-function(x,y,z,suffStat){

  x1=suffStat$data[,x];
  y1=suffStat$data[,y];
  z1=suffStat$data[,z];

  if (length(z)>0){
    out = gcm.test(x1,y1,z1)
  } else{
    out = gcm.test(x1,y1)
  }

  return(out$p.value)

}

