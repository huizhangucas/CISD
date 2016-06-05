Args=commandArgs()
dir_code=Args[6]
dir_wig=Args[7]
dir_out=Args[8]
chr=Args[9]


source(paste0(dir_code,"/write.table0.R"))
source(paste0(dir_code,"/extract_from_bed.R"))


FFT4 <- function(d1){ 
  length=ncol(d1)
  q=floor((length-1)/2)
  return(t(abs(apply(as.matrix(d1),1,fft)[c(1:(q+1)),])))
}


cat(paste0("Reading MNase-seq signal from ",dir_wig,"\n"))
cat("Calculating FFT for sites in\n")



  cat(chr)
  cat("...... ")
  inputfile=paste0(dir_wig,"/",chr,"_4.normalized")
  wig=as.numeric(read.table(inputfile)[,1])
  wig_len=length(wig)
  endpoint=floor((wig_len-100)/10)
  data=matrix(ncol=100,nrow=endpoint+1)
  for (k in 0:endpoint){
    start=k*10+1
    end=start+99
    data[k+1,]=wig[start:end]
  }

  fft_profile=round(FFT4(data))
  zero50=matrix(0,ncol=50,nrow=5)
  fft_profile=rbind(zero50,fft_profile,zero50)
  outputfile=paste0(dir_out,"/",chr,"_4.normalized.fft")
  write.table0(fft_profile,outputfile)

cat("Finished!\n")
cat(paste0("FFT results genomewide are saved to ",dir_out,"\n\n"))
