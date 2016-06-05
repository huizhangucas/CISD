Args=commandArgs()
model_file=Args[6]
candidate_file=Args[7]
output_file=Args[8]

library("e1071")
write.table0 <- function(x,file){
  write.table(x, file = file, append = FALSE, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = F, qmethod = c("escape", "double"),
              fileEncoding = "")
}



candidate_bed=read.table(candidate_file)

candidate_f2=rowMeans(candidate_bed[,5:6])-rowMeans(candidate_bed[,2:3])
candidate_f1=candidate_bed[,7]
candidate_feature=data.frame(cbind(candidate_f1,candidate_f2))
#Tr=read.table(trainingset2)
#colnames(Tr)=paste0("V",c(1:(ncol(Tr)-1),9999))
#model <- svm( V9999 ~ ., data = Tr)

#save(model,file="~/model2_SVM")
load(model_file)
colnames(candidate_feature)=paste0("V",c(1:2))
pred <- predict(model,candidate_feature)
result=candidate_bed[pred=="i",1:6]

write.table0(result,output_file)
