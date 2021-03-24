#加载相关的包
library(Matrix)
library(monocle)
library(patchwork)
library(Cairo)
library(ggplot2)
library(Seurat)
library(dplyr)

#清除环境
rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)
#设置路径
setwd("G:/Project/NTC_CD8/Monocle处理/")


#清洗只要特定组织(只要癌症组织)
genecount <- read.csv("Genecount_TTC.csv",header = T ,row.names = 1)
genemeta <- read.csv("TTC_monocleMeta.csv",header = T ,row.names = 1)

index1 <- genemeta$sampleType == "TTC"
genemeta1 <- genemeta[index1,]

index2 <- which(names(genecount) %in% rownames(genemeta1))

genecount <- genecount[,index2]

#write.csv(genecount,"Genecount_TTC.csv")
#write.csv(genemeta1,"Genemeta_TTC.csv")
rm(index2)

#进行monocle对象的构建
RA_matrix<-as(as.matrix(genecount), 'sparseMatrix')

feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))

rownames(feature_ann)<-rownames(RA_matrix)

RA_fd<-new("AnnotatedDataFrame", data = feature_ann)

sample_ann<- genemeta

rownames(sample_ann)<-colnames(RA_matrix)

RA_pd<-new("AnnotatedDataFrame", data =sample_ann)

#create monocle objection
RA.cds<-newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())

#查看phenodata、featuredata
head(pData(HSMM))

HSMM <- RA.cds

## 归一化 
HSMM <- estimateSizeFactors(HSMM)

HSMM <- estimateDispersions(HSMM)

#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 )

print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))

length(expressed_genes)

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 10, 
                        reduction_method = 'tSNE', verbose = T) 

HSMM <- clusterCells(HSMM, num_clusters = 9)

plot_cell_clusters(HSMM, color="Cluster")

## 挑选变异度大的基因，如图所示

HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM)
#saveRDS(HSMM,file = "TTC_monocle1021.rds")

HSMM <- readRDS("TTC_monocle1021.rds")

pData(HSMM) <- read.csv("TTC_monocleMeta.csv", header = T, row.names = 1)

#绘制拟时序图
pData(HSMM)$pseudotime <- as.factor(pData(HSMM)$pseudotime)
pData(HSMM)$seurat_clusters <- as.factor(pData(HSMM)$seurat_clusters)

plot_cell_trajectory(HSMM, color_by="pseudotime", show_tree=F, 
                     show_branch_points = F,show_cell_names =F,
                     show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#ADFF2F", "#32CD32", "#48D1CC", "red",
                                 "#4682B4", "#00008B", "#BA55D3","#8B008B"))


#保存meta信息和拟时序对象
write.csv(pData(HSMM),"TTC_monocleMeta.csv")


#SingleR注释看看效果
setwd("G:/Project/")
library(SingleR)


#提取Seurat相关的表达矩阵
HSMM_for_SingleR <- GetAssayData(HSMM,slot="data",generic = "CellDataSet")
?GetAssayData()
clusters <- HSMM@phenoData@data$Pseudotime


#refdata <- load("ref_Monaco_114s.RData")




pred.hesc <- SingleR(test = genecount, ref = ref_Monaco, 
                     labels = ref_Monaco$label.fine,
                     #因为样本主要为免疫细胞（而不是全部细胞），因此设置为label.fine
                     method = "cluster", clusters = clusters,
                     #这里我们为上一步分的9个cluster注释celltype
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels) 
plotScoreHeatmap(pred.hesc)


DimPlot(sce, reduction = "tsne",group.by = "seurat_clusters",label = T,pt.size = 1.2)


#tSNE可视化
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
#如下为sce对象注释细胞cluster鉴定结果。

sce@meta.data$celltype = "NA"

#先新增列celltype，值均为NA，然后利用下一行代码循环填充
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

DimPlot(sce, group.by="celltype", label=T, label.size=4, reduction='tsne',pt.size = 1.3)  
DimPlot(sce, group.by="seurat_clusters", label=T, label.size=4, reduction='tsne',pt.size = 1.3)


write.csv(sce@meta.data,"meta_singleR_seurat_v2.csv")


cg <- c("SRSF1","CCT6A","EIF4A3","IFI44L")
plot_genes_jitter(HSMM[cg,], grouping = "pseudotime",color_by = "pseudotime", plot_trend = T)+
scale_colour_manual(values = c("black", "black", "black", "black",
                                "black", "red", "black","black"))

?plot_genes_jitter()
plot_genes_in_pseudotime(HSMM[cg,], color_by = "pseudotime", min_expr = 2)


#制作DNB对象
test1 <- read.csv("TTC_monocleMeta.csv",header = T,row.names = 1)

test1 <- test1["pseudotime"]
test1 <- as.data.frame(t(test1))
havezero <- rbind(test1,genecount)
write.csv(havezero,"TTCAll_DNBcount_raw(labelled_pseudotime).csv")

