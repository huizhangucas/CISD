
Args=commandArgs()
dir_code=Args[6]
source(paste0(dir_code,"/write.table0.R"))

dir_fftwig=Args[7]
model_file=Args[8]
dir_p=Args[9]
chr=Args[10]
load(model_file)
cat("fft wig file is from this directory: ",dir_fftwig,"/\n")


  
  cat(chr)  
  cat("...... ")
  cat(paste0(dir_fftwig,"/",chr,"_4.normalized.fft\n"))
  inputfile=paste0(dir_fftwig,"/",chr,"_4.normalized.fft")
  
  outputfile=paste0(dir_p,"/",chr,"_p")
  fft_profile=read.table(inputfile)[,c(1,6,7)]
  
  len=nrow(fft_profile)+1
  
  
  f1=c(rep(0,10),as.numeric(fft_profile[,1]),rep(0,9))
  f6=c(rep(0,10),as.numeric(fft_profile[,2]),rep(0,9))
  f7=c(rep(0,10),as.numeric(fft_profile[,3]),rep(0,9))
  
  
  feature=matrix(ncol=3,nrow=len)
  
  for (k in 1:len){
    j=(k-1)*100+1
    
    feature[k,1]=mean(f1[k:(k+19)])
    feature[k,2]=mean(f6[k:(k+19)])
    feature[k,3]=mean(f7[k:(k+19)])
    
  }  
  feature=data.frame(feature)
  colnames(feature)=c("V1","V2","V3")
  p=predict(model0,feature,type="response")
  
  write.table0(cbind(p),outputfile)
cat("Finished!\n")
cat(paste0("the p of logistic regression genomewide is saved to ",dir_p,"/\n\n"))

