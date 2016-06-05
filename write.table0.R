write.table0 <- function(x,file){
  write.table(x, file = file, append = FALSE, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = F, qmethod = c("escape", "double"),
              fileEncoding = "")
}