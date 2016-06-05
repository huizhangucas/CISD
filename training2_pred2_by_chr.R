library("e1071")



Args=commandArgs()
dir_code=Args[6]
model_file=Args[7]
source(paste0(dir_code,"/write.table0.R"))
dir=Args[8]
chr=Args[9]
load(model_file)

candidate_bed=read.table(paste0(dir,"/candidate_",chr,".bed.f1to50"))[,1:3]
candidate_fft=read.table(paste0(dir,"/candidate_",chr,".bed.f1to50"))[,4:53]
candidate_feature=data.frame(candidate_fft)
colnames(candidate_feature)=paste0("V",c(1:50))
pred <- predict(model,candidate_feature)

result=candidate_bed[pred=="i",]
write.table0(result,paste0(dir,"/pred2_",chr))


