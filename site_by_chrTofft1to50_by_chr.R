

Args=commandArgs()
dir_code=Args[6]
source(paste0(dir_code,"/write.table0.R"))
source(paste0(dir_code,"/extract_feature_from_FftAndBed.R"))


dir_fftwig=Args[7]
dir_candidatesites=Args[8]
prefix=Args[9]
chr=Args[10]

cat("fft wig file is from this directory: ",dir_fftwig,"\n")
cat("extracting fft 1-50 for all sites in ",dir_candidatesites,"\n")
  
  cat(chr,"...... ")
  fft_file0=paste0(dir_fftwig,"/",chr,"_4.normalized.fft")
  cat(paste0("reading fft file:",fft_file0,"\n")) 
  fft_profile0=read.table(fft_file0)
  
  c0=read.table(paste0(dir_candidatesites,"/",prefix,"_",chr,".bed"))
  c0f1to50=extract_feature_from_FftAndBed_f_1to50(fft_profile0,c0)

  write.table0(cbind(c0,c0f1to50),paste0(dir_candidatesites,"/",prefix,"_",chr,".bed.f1to50"))  
cat("Finished!\n\n")
