library(ggplot2)
library(reshape2)
diff1<-read.table("diff_1.out", header=F)
diff2<-read.table("diff_2.out", header=F)
diff3<-read.table("diff_3.out", header=F)
diff4<-read.table("diff_4.out", header=F)
diff5<-read.table("diff_5.out", header=F)

diff1<-diff1[order(diff1$V1),]
diff2<-diff2[order(diff2$V1),]
diff3<-diff3[order(diff3$V1),]
diff4<-diff4[order(diff4$V1),]
diff5<-diff5[order(diff5$V1),]

diff1<-head(diff1,100)
diff2<-head(diff2,100)
diff3<-head(diff3,100)
diff4<-head(diff4,100)
diff5<-head(diff5,100)

diffs<-cbind.data.frame(diff1$V1, diff1$V2, diff2$V2, diff3$V2, diff4$V2, diff5$V2)


colnames(diffs)<-c("genes", "v1", "v2", "v3", "v4", "v5")

diffs.m<-melt(diffs)

pdf("compare.pdf")
ggplot(diffs.m, aes(variable, value, group=genes, color=genes)) + geom_line() + geom_point() +  theme(legend.position="none")
dev.off()
