#清除环境
rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)
#设置路径
setwd("G:/Project/NTC_CD8/Monocle处理/")
genecount <- read.csv("TTCAll_DNBcount_raw(labelled_pseudotime).csv",header = T ,row.names = 1)
names(genecount) <- genecount[1,]
genecount <- genecount[-1,]


index_deg <- which(rownames(genecount) %in% Deg_all)

genecount_deg <- genecount[index_deg,]

genecount_deg = genecount

index1 <- names(genecount_deg) ==1
index2 <- names(genecount_deg) ==2
index3 <- names(genecount_deg) ==3
index4 <- names(genecount_deg) ==4
index5 <- names(genecount_deg) ==5
index6 <- names(genecount_deg) ==6
index7 <- names(genecount_deg) ==7
index8 <- names(genecount_deg) ==8

genecount_1 <- genecount_deg[,index1]
genecount_2 <- genecount_deg[,index2]
genecount_3 <- genecount_deg[,index3]
genecount_4 <- genecount_deg[,index4]
genecount_5 <- genecount_deg[,index5]
genecount_6 <- genecount_deg[,index6]
genecount_7 <- genecount_deg[,index7]
genecount_8 <- genecount_deg[,index8]

Cluster1 <- apply(genecount_1,1,mean)
Cluster2 <- apply(genecount_2,1,mean)
Cluster3 <- apply(genecount_3,1,mean)
Cluster4 <- apply(genecount_4,1,mean)
Cluster5 <- apply(genecount_5,1,mean)
Cluster6 <- apply(genecount_6,1,mean)
Cluster7 <- apply(genecount_7,1,mean)
Cluster8 <- apply(genecount_8,1,mean)

deg_mathix <-data.frame(Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster7,Cluster8)
deg_mathix <- as.data.frame(t(deg_mathix))


#need package GSVA
library(GSVA)

library(GSEABase)


library(limma)
#设置工作路径
setwd("G:/Project/NTC_CD8/GSVA/")
getwd()
#构建GMT文件，我彩玉的是免疫的
geneSets <- getGmt("immune_test1.gmt")
#geneSets <- getGmt("Data/kegg.v7.1.symbols.gmt")

mydata <- read.csv(file = "G:/Project/NTC_CD8/Monocle处理/TTCAll_DNBcount_raw(labelled_pseudotime).csv",header = T)

mydata = deg_mathix
name = mydata[,1]
index <- duplicated(mydata[,1])
fildup = mydata[!index,]
exp = fildup[,-1]
row.names(exp) = fildup[,1]

exp <- exp[-1,]

mydata= as.matrix(t(deg_mathix))
mydata=round(mydata)
#这一步很耗时
Es <- gsva(mydata, geneSets, parallel.sz=1,kcdf="Poisson")

pheatmap::pheatmap(Es,angle_col = 0,cluster_cols = F,cluster_rows = T)


library(pheatmap)
Es_1 <-  t(Es)
pheatmap(Es, cluster_cols = F, cluster_rows = T, border=FALSE,
         #color = colorRampPalette(c("navy","white","red"))(50),
         color = colorRampPalette(c("black","#191970","#00008B","#1E90FF",
                                    "#F0F8FF","#FF4500","red","#8B0000","black"))(100),
         scale = "row",show_colnames =F,show_rownames = T,
         main = "Expression in 8 different Pseudotime GSVA", angle_col = 0,
         treeheight_col = 9,
         fontsize = 14,
         fontsize_row = 13,
         fontsize_col = 13)








write.csv(Es,"GSVA_ES.csv")

#limma 宸紓鍩哄洜闆嗗垎鏋?
adjPvalueCutoff <- 0.05
logFCcutoff <- log2(2)
#璁剧疆鍒嗙粍
class<-c("treat","treat","control","control")
design<-model.matrix(~factor(class))
colnames(design)<-c("control","treat")

fit <- lmFit(Es, design)
fit <- eBayes(fit)

allgeneSets<-topTable(fit,coef = 2,adjust="fdr",number = 200000)
DEgeneSets <- topTable(fit, coef=2, number=Inf,p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)


#鎻愬彇宸紓鍩哄洜闆嗙储寮曞彲瑙嗗寲
#DEgeneSetsgenes<-row.names(DEgeneSets)
#index<-which(row.names(leukemia_es@assayData$exprs)%in%DEgeneSetsgenes)
#DEgeneSetsexp<-leukemia_es@assayData$exprs[index,]
#pheatmap::pheatmap(DEgeneSetsexp)


mydata<-read.csv("Data/Bcell_GSVA.csv",row.names = 1)
#name=mydata[,1]
#index <- duplicated(mydata[,1])
#fildup=mydata[!index,]
#mydata<-fildup
#write.csv(mydata,"Data/Bcell_GSVA.csv")
pheatmap::pheatmap(mydata[20:49,],cluster_cols = F,angle_col = 0)





