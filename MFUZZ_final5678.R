#清除环境
rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)
#设置路径
setwd("G:/Project/NTC_CD8/Monocle处理/DNB_1_jie/")

DNB_YIJIE <- read.table("DNB_background_gene_network.txt")
TOP50 <- read.table("TOP50gene.txt")

index1 <- which(DNB_YIJIE$V1 %in% TOP50$V1)

Connect1JIE <- DNB_YIJIE[index1,]

#读取表达矩阵
setwd("G:/Project/NTC_CD8/Monocle处理/DNB_1_jie/")

genecount <- read.csv("TTCAll_DNBcount_raw(labelled_pseudotime).csv",header = T, row.names = 1)

index2 <- which(row.names(genecount) %in% Connect1JIE$V2)

genecountMFUZZ <- genecount[index2,]

genecountMFUZZ <- rbind(genecount[1,],genecountMFUZZ)

names(genecountMFUZZ) <- genecountMFUZZ[1,]

genecountMFUZZ <- genecountMFUZZ[-1,]

write.csv(genecountMFUZZ,"1jie_countexp.csv")

#取平均
index5 <- names(genecountMFUZZ) ==5
index6 <- names(genecountMFUZZ) ==6
index7 <- names(genecountMFUZZ) ==7
index8 <- names(genecountMFUZZ) ==8
genecount_5 <- genecountMFUZZ[,index5]
genecount_6 <- genecountMFUZZ[,index6]
genecount_7 <- genecountMFUZZ[,index7]
genecount_8 <- genecountMFUZZ[,index8]

Cluster5 <- apply(genecount_5,1,mean)
Cluster6 <- apply(genecount_6,1,mean)
Cluster7 <- apply(genecount_7,1,mean)
Cluster8 <- apply(genecount_8,1,mean)

MFUZZ_mean_exp <- data.frame(Cluster5,Cluster6,Cluster7,Cluster8)
write.csv(MFUZZ_mean_exp,"MFUZZ_5678_mean_exp.csv")

#Mfuzz软聚类 
library(Mfuzz)

set.seed(2020)
MFUZZ_mean_exp <- read.csv("MFUZZ_5678_mean_exp.csv", header = T, row.names = 1)
MFUZZ_mean_exp <- as.matrix(MFUZZ_mean_exp)


dat <- new('ExpressionSet',exprs = MFUZZ_mean_exp)
yeast.r <- filter.NA(dat, thres=0.25)
yeast.f <- fill.NA(dat,mode="mean")
tmp <- filter.std(yeast.f,min.std=0)
yeast.s <- standardise(tmp)

cl <- mfuzz(yeast.s,c=20,m=1.25)
mfuzz.plot(yeast.s,cl=cl,mfrow=c(4,5))

#保存MFUZZ对象
saveRDS(cl,"cl_v1.rds")
cl <- readRDS("cl_v1.rds")
#保存标准化的文件
a <- yeast.s@assayData$exprs
write.csv(a,"exp_normalized_MFUZZ.csv")


#每个簇下基因数量

cl$size

#每个基因所属簇

head(cl$cluster)

#基因和 cluster 之间的 membership，用于判断基因所属簇，对应最大值的那个簇

head(cl$membership)

#整合关系输出

gene_cluster <- cbind(cl$cluster, cl$membership)

colnames(gene_cluster)[1] <- 'cluster'

write.csv(gene_cluster, 'MFUZZgene_20cluster.csv')



