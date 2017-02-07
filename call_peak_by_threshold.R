call_peak_by_threshold=function(wig,threshold,chr){  
#  result=data.frame(matrix(ncol=3))
  result=matrix(ncol=3)
  wig=as.numeric(wig[,1])
  flag=0 #### 0:out of the peak  1:in the peak
  for(i in 1:(length(wig)-2)){
    if(wig[i]>threshold & flag==0){
      start=i
      flag=1
    }else if(wig[i]<threshold & flag==1){
      end=i-1
      flag=0
      result=rbind(result,c(chr,start*100+1,end*100+2))
    }
  }
  return(result[2:nrow(result),])
}