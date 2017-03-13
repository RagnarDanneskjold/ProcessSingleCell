args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("usage: myout fcout file_num")

}

file_num <- args[3]

myout<-read.table(args[1], header=FALSE)
fcout<-read.table(args[2], header=FALSE)

diff<-myout$V2 - fcout$V2
diff<-cbind.data.frame(myout$V1, diff)
write.table(file=paste0("diff_", file_num, ".out"),diff[order(diff$diff, decreasing=TRUE),],col.names=FALSE,row.names=FALSE,quote=FALSE, sep='\t')
