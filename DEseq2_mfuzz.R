#清除环境
rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)

##设置路径


setwd("G:/Project/NTC_CD8/Monocle处理")

setwd("G:/Project/NTC_CD8/Monocle处理/DEseq2/")

count <- read.csv("TTCAll_DNBcount_raw(labelled_pseudotime).csv", row.names = 1, header = T)
count <- as.data.frame(t(count))
index5 <- count$pseudotime == 5
index6 <- count$pseudotime == 6
index7 <- count$pseudotime == 7
index8 <- count$pseudotime == 8

count_5 <- count[index5,]
count_6 <- count[index6,]
count_7 <- count[index7,]
count_8 <- count[index8,]

newcount <- rbind(count_5,count_6,count_7,count_8)
newcount <- as.data.frame(t(newcount))

index1 <- newcount["pseudotime",]==5
index2 <- newcount["pseudotime",]==6
index3 <- newcount["pseudotime",]==7
index4 <- newcount["pseudotime",]==8

count5 <- newcount[,index1]
count6 <- newcount[,index2]
count7 <- newcount[,index3]
count8 <- newcount[,index4]
#构建dds对象,开始DESeq流程
library(DESeq2)

#构建dds对象需要的condition,添加分组信息
condition <- factor(c(rep("control",160),rep("treat",168)), 
                    levels = c("control","treat"))

Deseq6VS8 <- cbind(count6,count8)

coldata<-data.frame(row.names=colnames(Deseq6VS8), condition)
#构建DDS对象
dds <- DESeqDataSetFromMatrix(Deseq6VS8, coldata, design= ~ condition)

dds<-DESeq(dds)

res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
write.csv(res,file="6VS8_result.csv")
diff_gene <-subset(res, padj < 0.05 & abs(log2FoldChange) > 2)
write.csv(diff_gene,file = "6vs8_diff_gene.csv")


Deseq6VS7 <- cbind(count6,count7)
condition <- factor(c(rep("control",160),rep("treat",168)), 
                    levels = c("control","treat"))
coldata<-data.frame(row.names=colnames(Deseq6VS7), condition)
#构建DDS对象

dds <- DESeqDataSetFromMatrix(Deseq6VS7, coldata, design= ~ condition)

dds<-DESeq(dds)

res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
write.csv(res,file="6VS7_result.csv")
diff_gene <-subset(res, padj < 0.05 & abs(log2FoldChange) > 2)
write.csv(diff_gene,file = "6vs7_diff_gene.csv")

#合并6与7和8的差异基因，去除重复
diff678 <- read.table("6vs78diffgene.txt")

index1 <- duplicated(diff678$V1)
diff678clean <- diff678[!index1,]

write.csv(diff678clean,"diff678clean.csv")
