Args=commandArgs()
file=Args[6]

a=as.numeric(read.table(file)[,1])
cat(paste0(sd(a),"\n"))
