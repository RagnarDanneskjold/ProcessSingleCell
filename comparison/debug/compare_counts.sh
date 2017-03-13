#!/bin/bash

if [ $# -ne 4 ]
then
    echo "usage: <htseq output> <featureCounts> <my output> <run number>"
    exit 1
fi

htseq=$1
fcounts=$2
mine=$3
run=$4

echo "grabbing correct lines"

cut -f 1,7 $fcounts | sed 1,2d | sort -k1,1 > fcounts_$run
cut -f 5,6 $mine | sort -k1,1 > mine_$run
head -n -5 $htseq | sort -k1,1 > htseq_$run

echo "testing wc"

fc_lines=$(wc -l < fcounts_$run)
m_lines=$(wc -l < mine_$run)
ht_lines=$(wc -l < htseq_$run)

if [ $fc_lines -ne $m_lines ] || [ $fc_lines -ne $ht_lines ] || [ $ht_lines -ne $m_lines ] 
then
    echo "Lines don't match: fc: $fc_lines mine:$m_lines htseq:$ht_lines"
    exit 1
fi

echo "printing R code"

r_code="library(ggplot2)
library(reshape2)

"


cbind_fc="diffs_fc<-cbind.data.frame(diff_fc_1\$genes"
cbind_ht="diffs_ht<-cbind.data.frame(diff_ht_1\$genes"
cbind_cntrl="diffs_cntrl<-cbind.data.frame(diff_cntrl_1\$genes"
colnames=")<-c(\"genes\""

for i in {1..$run}
do
    r_code+="
ht_$run<-read.table(\"htseq_$run\",header=F)
fc_$run<-read.table(\"fcounts_$run\",header=F)
m_$run<-read.table(\"mine_$run\",header=F)

diff_fc_$run<-cbind.data.frame(m_$run\$V1, m_$run\$V2 - fc_$run\$V2)
diff_ht_$run<-cbind.data.frame(m_$run\$V1, m_$run\$V2 - ht_$run\$V2)
diff_cntrl_$run<-cbind.data.frame(fc_$run\$V1, fc_$run\$V2 - ht_$run\$V2)

colnames(diff_fc_$run)<-c(\"genes\", \"vals\")
colnames(diff_ht_$run)<-c(\"genes\", \"vals\")
colnames(diff_cntrl_$run)<-c(\"genes\", \"vals\")

diff_fc_$run<-head(diff_fc_$run[order(diff_fc_$run\$genes),],100)
diff_ht_$run<-head(diff_ht_$run[order(diff_ht_$run\$genes),],100)
diff_cntrl_$run<-head(diff_cntrl_$run[order(diff_cntrl_$run\$genes),],100)

    "

    colnames+=",\"v$run\""

    cbind_fc+=",diff_fc_$run\$vals"
    cbind_ht+=",diff_ht_$run\$vals"
    cbind_cntrl+=",diff_cntrl_$run\$vals"
done

r_code+="
$cbind_fc)
$cbind_ht)
$cbind_cntrl)

colnames(diffs_fc$colnames)
colnames(diffs_ht$colnames)
colnames(diffs_cntrl$colnames)

diffs_fc.m<-melt(diffs_fc)

pdf(\"compare_fc.pdf\")

ggplot(diffs_fc.m, aes(variable, value, group=genes)) + geom_line() + geom_point() + theme(legend.position=\"none\")

dev.off()

diffs_ht.m<-melt(diffs_ht)

pdf(\"compare_ht.pdf\")

ggplot(diffs_ht.m, aes(variable, value, group=genes)) + geom_line() + geom_point() + theme(legend.position=\"none\")

dev.off()

diffs_cntrl.m<-melt(diffs_cntrl)

pdf(\"compare_cntrl.pdf\")

ggplot(diffs_cntrl.m, aes(variable, value, group=genes)) + geom_line() + geom_point() + theme(legend.position=\"none\")

dev.off()
"

echo "running R code"

echo "$r_code" > rscript_$run.R
Rscript rscript_$run.R

echo "done"
