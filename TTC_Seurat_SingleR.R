#清除环境
rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)

#Seurat流程

#调用相关包
library(dplyr)
library(Seurat)
library(patchwork)
library(Cairo)
library(ggplot2)

setwd("G:/Project/NTC_CD8/肿瘤组织样本CD8/")

genecount <- read.csv("Genecount_TTC.csv",header = T ,row.names = 1)
genemeta <- read.csv("Genemeta_TTC.csv",header = T ,row.names = 1)
TTC <- CreateSeuratObject(genecount,
                            meta.data = genemeta,
                            project = "TTC",
                            min.cells = 5,
                            min.features  = 200)

#数据清洗，tsne降维对meta信息进行添加
#Normalizing the data

TTC <- NormalizeData(TTC, normalization.method = "LogNormalize", scale.factor = 10000)
TTC <- NormalizeData(TTC)

#Identification of highly variable features (feature selection)
TTC <- FindVariableFeatures(TTC, selection.method = "vst", nfeatures = 2000)
?FindVariableFeatures()
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(TTC), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(TTC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#Scaling the data
all.genes <- rownames(TTC)
TTC <- ScaleData(TTC, features = all.genes)

#Perform linear dimensional reduction

TTC <- RunPCA(TTC, features = VariableFeatures(object = TTC))

VizDimLoadings(TTC, dims = 1:2, reduction = "pca")

DimPlot(TTC, reduction = "pca")

DimHeatmap(TTC, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(TTC, dims = 1:6, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset

TTC <- JackStraw(TTC, num.replicate = 300)

TTC <- ScoreJackStraw(TTC, dims = 1:20)

JackStrawPlot(TTC, dims = 1:15)

ElbowPlot(TTC)

##Cluster the cells


TTC <- FindNeighbors(TTC, dims =1:15)
TTC <- FindClusters(TTC, resolution = 0.9)


# Look at cluster IDs of the first 5 cells
head(Idents(TTC), 5)

#Run non-linear dimensional reduction (umap)

TTC <- RunTSNE(TTC, dims = 1:10)

write.csv(TTC@meta.data,"TTC_seurat_meta.csv")
TTC@meta.data <- read.csv("TTC_seurat_meta.csv", header = T, row.names = 1)
DimPlot(TTC, reduction = "tsne",group.by = "Pseudotime",label = T,pt.size = 1.2)

TTC@meta.data


cluster1.markers <- FindMarkers(TTC, ident.1 = 1, min.pct = 0.25,group.by = "pseudotime")
cluster2.markers <- FindMarkers(TTC, ident.1 = 2, min.pct = 0.25,group.by = "pseudotime")
cluster3.markers <- FindMarkers(TTC, ident.1 = 3, min.pct = 0.25,group.by = "pseudotime")
cluster4.markers <- FindMarkers(TTC, ident.1 = 4, min.pct = 0.25,group.by = "pseudotime")
cluster5.markers <- FindMarkers(TTC, ident.1 = 5, min.pct = 0.25,group.by = "pseudotime")
cluster6.markers <- FindMarkers(TTC, ident.1 = 6, min.pct = 0.25,group.by = "pseudotime")
cluster7.markers <- FindMarkers(TTC, ident.1 = 7, min.pct = 0.25,group.by = "pseudotime")
cluster8.markers <- FindMarkers(TTC, ident.1 = 8, min.pct = 0.25,group.by = "pseudotime")

setwd("G:/Project/NTC_CD8/肿瘤组织样本CD8/耗竭性Marker差异分析/")
write.csv(cluster1.markers,"deg1.csv")
write.csv(cluster2.markers,"deg2.csv")
write.csv(cluster3.markers,"deg3.csv")
write.csv(cluster4.markers,"deg4.csv")
write.csv(cluster5.markers,"deg5.csv")
write.csv(cluster6.markers,"deg6.csv")
write.csv(cluster7.markers,"deg7.csv")
write.csv(cluster8.markers,"deg8.csv")
