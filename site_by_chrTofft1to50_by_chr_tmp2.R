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

write.table0 <- function(x,file){
  write.table(x, file = file, append = FALSE, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = F, qmethod = c("escape", "double"),
              fileEncoding = "")
}


extract_feature_from_FftAndBed_f_1to50 = function(fft_profile,bed){
  
  f6=rowMeans(extract_from_bed(data.frame(fft_profile[,6]),bed,2000,100))  
  f7=rowMeans(extract_from_bed(data.frame(fft_profile[,7]),bed,2000,100))
  
  f1=rowMeans(extract_from_bed(data.frame(fft_profile[,1]),bed,2000,100))
  f2=rowMeans(extract_from_bed(data.frame(fft_profile[,2]),bed,2000,100))
  f3=rowMeans(extract_from_bed(data.frame(fft_profile[,3]),bed,2000,100))
  f4=rowMeans(extract_from_bed(data.frame(fft_profile[,4]),bed,2000,100))
  f5=rowMeans(extract_from_bed(data.frame(fft_profile[,5]),bed,2000,100))
  f8=rowMeans(extract_from_bed(data.frame(fft_profile[,8]),bed,2000,100))
  f9=rowMeans(extract_from_bed(data.frame(fft_profile[,9]),bed,2000,100))
  f10=rowMeans(extract_from_bed(data.frame(fft_profile[,10]),bed,2000,100))
  
  f11=rowMeans(extract_from_bed(data.frame(fft_profile[,11]),bed,2000,100))
  f12=rowMeans(extract_from_bed(data.frame(fft_profile[,12]),bed,2000,100))
  f13=rowMeans(extract_from_bed(data.frame(fft_profile[,13]),bed,2000,100))
  f14=rowMeans(extract_from_bed(data.frame(fft_profile[,14]),bed,2000,100))
  f15=rowMeans(extract_from_bed(data.frame(fft_profile[,15]),bed,2000,100))
  f16=rowMeans(extract_from_bed(data.frame(fft_profile[,16]),bed,2000,100))
  f17=rowMeans(extract_from_bed(data.frame(fft_profile[,17]),bed,2000,100))
  f18=rowMeans(extract_from_bed(data.frame(fft_profile[,18]),bed,2000,100))
  f19=rowMeans(extract_from_bed(data.frame(fft_profile[,19]),bed,2000,100))
  f20=rowMeans(extract_from_bed(data.frame(fft_profile[,20]),bed,2000,100))
  
  f21=rowMeans(extract_from_bed(data.frame(fft_profile[,21]),bed,2000,100))
  f22=rowMeans(extract_from_bed(data.frame(fft_profile[,22]),bed,2000,100))
  f23=rowMeans(extract_from_bed(data.frame(fft_profile[,23]),bed,2000,100))
  f24=rowMeans(extract_from_bed(data.frame(fft_profile[,24]),bed,2000,100))
  f25=rowMeans(extract_from_bed(data.frame(fft_profile[,25]),bed,2000,100))
  f26=rowMeans(extract_from_bed(data.frame(fft_profile[,26]),bed,2000,100))
  f27=rowMeans(extract_from_bed(data.frame(fft_profile[,27]),bed,2000,100))
  f28=rowMeans(extract_from_bed(data.frame(fft_profile[,28]),bed,2000,100))
  f29=rowMeans(extract_from_bed(data.frame(fft_profile[,29]),bed,2000,100))
  f30=rowMeans(extract_from_bed(data.frame(fft_profile[,30]),bed,2000,100))
  
  f31=rowMeans(extract_from_bed(data.frame(fft_profile[,31]),bed,2000,100))
  f32=rowMeans(extract_from_bed(data.frame(fft_profile[,32]),bed,2000,100))
  f33=rowMeans(extract_from_bed(data.frame(fft_profile[,33]),bed,2000,100))
  f34=rowMeans(extract_from_bed(data.frame(fft_profile[,34]),bed,2000,100))
  f35=rowMeans(extract_from_bed(data.frame(fft_profile[,35]),bed,2000,100))
  f36=rowMeans(extract_from_bed(data.frame(fft_profile[,36]),bed,2000,100))
  f37=rowMeans(extract_from_bed(data.frame(fft_profile[,37]),bed,2000,100))
  f38=rowMeans(extract_from_bed(data.frame(fft_profile[,38]),bed,2000,100))
  f39=rowMeans(extract_from_bed(data.frame(fft_profile[,39]),bed,2000,100))
  f40=rowMeans(extract_from_bed(data.frame(fft_profile[,40]),bed,2000,100))
  
  f41=rowMeans(extract_from_bed(data.frame(fft_profile[,41]),bed,2000,100))
  f42=rowMeans(extract_from_bed(data.frame(fft_profile[,42]),bed,2000,100))
  f43=rowMeans(extract_from_bed(data.frame(fft_profile[,43]),bed,2000,100))
  f44=rowMeans(extract_from_bed(data.frame(fft_profile[,44]),bed,2000,100))
  f45=rowMeans(extract_from_bed(data.frame(fft_profile[,45]),bed,2000,100))
  f46=rowMeans(extract_from_bed(data.frame(fft_profile[,46]),bed,2000,100))
  f47=rowMeans(extract_from_bed(data.frame(fft_profile[,47]),bed,2000,100))
  f48=rowMeans(extract_from_bed(data.frame(fft_profile[,48]),bed,2000,100))
  f49=rowMeans(extract_from_bed(data.frame(fft_profile[,49]),bed,2000,100))
  f50=rowMeans(extract_from_bed(data.frame(fft_profile[,50]),bed,2000,100))
  return(cbind(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26,f27,f28,f29,f30,f31,f32,f33,f34,f35,f36,f37,f38,f39,f40,f41,f42,f43,f44,f45,f46,f47,f48,f49,f50))
}



Args=commandArgs()


fftwig=Args[6]
candidatesites=Args[7]
outbedwithfeature=Args[8]



  fft_file0=paste0(fftwig)
  fft_profile0=read.table(fft_file0)

  c0=read.table(paste0(candidatesites))
  c0f1to50=extract_feature_from_FftAndBed_f_1to50(fft_profile0,c0)

  write.table0(cbind(c0,c0f1to50),outbedwithfeature)

