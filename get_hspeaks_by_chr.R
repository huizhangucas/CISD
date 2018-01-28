Args=commandArgs()
dir_code=Args[6]
source(paste0(dir_code,"/write.table0.R"))
#source("~/project1/code/call_peak_by_threshold.R")

threshold = Args[7]
dir_p_wig=Args[8]
dir_output=Args[9]
chr=Args[10]

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

cat("calling high score peaks\n")
cat("p wig files are from this dir:",dir_p_wig,"\n")
  
  cat(chr,"...... ") 
  inputfile=paste0(dir_p_wig,"/",chr,"_p")
  outputfile=paste0(dir_output,"/hspeaks",threshold,chr,".bed")
  
  p=read.table(inputfile)
  
  hspeaks=call_peak_by_threshold(wig = p,threshold,chr)
  write.table0(x=hspeaks,outputfile)
cat(" Finished!\n")
