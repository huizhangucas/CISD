Args=commandArgs()

inputfile=Args[6] # krnormalized file, dict format, 2 fields: anchor1_anchor2 krnormed_value
outputfile=Args[7] # output file, dict format, 2 fields: anchor1_anchor2 zscore

cat(paste0("calculating zscore of ",inputfile))
cat("\n")
a=read.table(inputfile)

zscore=scale(log(a[,2]))


write.table0 <- function(x,file){
  write.table(x, file = file, append = FALSE, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = F, qmethod = c("escape", "double"),
              fileEncoding = "")
}
write.table0(cbind(a[,1],zscore),outputfile)
