RCoT_wrap<-function(x,y,z,suffStat){
  
  x1=suffStat$data[,x];
  y1=suffStat$data[,y];
  z1=suffStat$data[,z];
  
  out = RCIT:::RCoT(x1,y1,z1,num_f=100);
  
  return(out$p)
  
}

