extract_from_bed = function(wig,bed,win_len=2000,fold=10){
  wig=as.numeric(wig[,1])
  signal=matrix(ncol=win_len/fold,nrow=nrow(bed))
  offset=win_len/(2*fold)
  for(i in 1:nrow(bed)){
    center0=(as.numeric(bed[i,2])+as.numeric(bed[i,3]))/(2*fold)
    center=round(center0)
    flag=center0-center
    if(flag>=0){
      start=center-offset+2
      end=center+offset+1
      signal[i,]=wig[start:end]
    }else{
      start=center-offset+1
      end=center+offset
      signal[i,]=wig[start:end]
    }  
  }
  return(signal)
}